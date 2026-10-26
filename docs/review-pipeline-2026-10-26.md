# Pipeline review against current single-cell practice — 2026-10-26

Scope: `pipeline/` as of commit 4759b4d, read-only. Reference guidance taken from the
scanpy, anndata, scvi-tools, scvelo and pathway-enrichment skills (Scientific Agent Skills
v2.72.0) and translated to the Seurat/R equivalents the pipeline actually uses. Every claim
cites a `file:line` that was read; anything not read is marked [unverified]. Decisions
already made (frozen reference, Harmony, cm `min_features` 1000, H1D20_1 exclusion) are not
re-litigated.

## 1. Verdict

The pipeline is a sound, reproducible Seurat 5 clustering-and-annotation workflow with
unusually good provenance (cell-fate accounting, version tables, frozen-reference drift
benchmark). Its weakest point is statistics: every between-sample comparison is per-cell
Wilcoxon or a within-sample bootstrap, which treats cells as replicates and cannot support
a condition-level claim. Second, QC thresholds are fixed constants chosen per species rather
than derived from each sample's distribution, and ambient-RNA correction is a manual
side-step rather than a pipeline stage.

## 2. Stage-by-stage

| Stage | Current approach | Gap or risk | Recommendation | Effort |
|---|---|---|---|---|
| Load / gene filter | `CreateSeuratObject(min.cells = 3, min.features = 200)` `01_load_qc.R:39-40` | `min.features = 200` is hard-coded in the loader, independent of `QC$min_features` (`config.R:129`). cm sets `QC$min_features` 1000, so cells with 200–999 genes are loaded, counted in `cell_fate.csv` as "After_load", then dropped at QC. Harmless but it is what produced the H1D20_2 "low QC" false flag. | Pass `min.features = QC$min_features` to the loader, or document that the two are intentionally different. | S |
| QC thresholds | Fixed per species: 200/5000/500/25000/20 % `config.R:129-133`; cm and mouse overlays override by hand | No data-derived option. The scanpy skill's rule is "inspect the distribution before choosing thresholds rather than copying defaults". Nothing in `01_load_qc.R` computes MAD or quantile cut-offs (grep for `mad(`/`quantile(` returned nothing). The mouse mito median was 0.4–1.6 % against a 30 % ceiling, so that cut is inert. | Add an optional `QC$adaptive = TRUE` branch: per-sample median ± 3 MAD on log counts, log genes and %mt, with the fixed values as floor/ceiling. Write the realised thresholds to `qc/qc_thresholds.csv` so the report can show them. | M |
| Ambient RNA | `soupx_adjust.R:90-94` runs `autoEstCont(forceAccept = TRUE)` and writes corrected matrices to `Samples_soupx/<s>_sx/`; report reads `qc/soupx_rho.csv` `08b_html_report.R:670` | Not a pipeline step: `run_pipeline.sh` has no reference to soupx. The corrected cohort is a second set of sample dirs, so cache keys, run names (`_sx`) and the cell-fate table all fork. `forceAccept = TRUE` suppresses SoupX's own sanity check on rho; H1D20_1 at rho 0.795 would otherwise have errored. | Make SoupX an optional step `00` in `run_pipeline.sh` that writes into the run dir, keyed into `cache_hash()`. Record whether `forceAccept` fired. Decide a rho ceiling (e.g. 0.3) above which the sample is flagged, not corrected. | M |
| Doublets | scDblFinder per sample, `dbr` auto unless listed, `nfeatures = 1000` `02_doublets.R:40-45`; seed set `02_doublets.R:39` | Sound. `clusters=` is not passed, so scDblFinder uses its random-pairs mode; the cluster-based mode is more sensitive for heterogeneous tissue (cm, mouse) [unverified for this data]. `DOUBLET$doublet_rate` lists only H1/H2 `config.R:148`. | Try `clusters = TRUE` on one cm sample and compare called rate. Drop the H1/H2 hard-coding or move it to a per-project overlay. | S |
| Normalisation | LogNormalize, 10 000, vst 2000 HVG `config.R:150-158`, `03_individual.R:42-47`; `ScaleData(vars.to.regress = "percent.mt")` `03_individual.R:56`, repeated on the merged object `04_integrate.R:59-62` | Regressing %mt is fine for blood; for cardiomyocytes, mito load is biology (the cm overlay comment at `config_species_mouse.R:31-34` says so) and regressing it removes the CM vs progenitor axis. HVGs are recomputed on the merged object without `selection.method = "vst"` or per-sample pooling (`04_integrate.R:59`), so highly expressed sample-specific genes dominate. | Make `vars.to.regress` a `NORM` field and set it to `NULL` in the cm overlay. Select HVGs per layer (`FindVariableFeatures` on split layers, Seurat 5 default) before `JoinLayers` at `04_integrate.R:107`. | S |
| Cell cycle | None. grep for `CellCycleScoring`/`S.Score` returned nothing | For a differentiation time-course (cm) and a proliferating cluster (mouse cluster 8, Mki67 77 %), phase is a confounder the clustering will pick up. | Add `CellCycleScoring` in step 03, store S/G2M scores, show phase on the UMAP in 06, and offer regression as a `NORM` option. | S |
| Dimensionality | `npcs = 30`, `dims_use = 1:20` fixed `config.R:164,192`; elbow plot drawn `03_individual.R:59-63` | The red line is drawn at a constant, not read from the curve. 20 PCs is low for tissue cohorts with 15+ cell types. | Compute the elbow (first PC where cumulative variance > 90 % and step < 0.1 %) and log it; keep the fixed value as default but warn when they disagree by > 5. | S |
| Integration | Harmony on `sample`, theta 2, lambda 1 `config.R:187-189`, `04_integrate.R:142-146`; over-correction warning `04_integrate.R:210` | Batch and biology are confounded by design in two cohorts: flowcell = timepoint in mouse, sample = condition in the 2-sample DE path. Harmony on `sample` removes the condition signal before DE. The scvi-tools skill: "do not add condition as a nuisance covariate automatically". No integration metric is computed (LISI, silhouette per batch, kBET) [grep for `silhouette` returned nothing]. | Add `HARMONY$group_by_vars` override per project (e.g. `replicate` or `flowcell`, never the contrast variable). Compute iLISI/cLISI or a per-cell-type batch silhouette and put the numbers in the report. | M |
| Clustering | Louvain (`algorithm = 1`) at 0.3–0.8, default 0.5 `config.R:171-174`; comparison UMAP `04_integrate.R:234-250` | Leiden (`algorithm = 4`) is the current default everywhere else (scanpy deprecates Louvain). Default resolution is a constant, not chosen per dataset. | Switch to `algorithm = 4` (needs `leidenbase`); add a clustree or silhouette-vs-resolution panel to the comparison PDF so the choice is evidence-based. | S |
| Annotation | SingleR per cell, fine-tune + prune `05_annotate.R:67-68`; delta stored as annotation only `05_annotate.R:282-283`; scType marker route; manual cluster map; frozen reference via 05r | Strong design. Two gaps: (a) SingleR runs per cell, not per cluster, so unmapped clusters shatter (`CLAUDE.md` gotcha) — a cluster-level SingleR pass (`clusters =`) would give a stable fallback; (b) cm marker panels are also used for mouse (`config_species_mouse.R:67-90`), so scType scores an immune dataset against Pluripotent/Cardiomyocyte panels. | Add a cluster-level SingleR call alongside the per-cell one and use it as the default for unmapped clusters. Give mouse immune panels (Monocyte, Neutrophil, DC, B, T, NK, Macrophage). | S/M |
| Marker genes | `FindAllMarkers` Wilcoxon, `min.pct 0.1`, `logfc 0.25` `03_individual.R:103-104`, `06_visualize.R:269-270` | Fine as exploratory markers. The scanpy skill's caveat applies: clusters chosen from the same data yield exploratory markers, p-values are not inferential. No `test.use = "MAST"` or pseudobulk option. | Label the marker tables "exploratory" in the report glossary. No code change needed. | S |
| Condition DE | Per-cell-type `FindMarkers(ident.1 = s1, ident.2 = s2, test.use = "wilcox")` `06b_differential.R:103-108`; supports exactly 2 samples `06b_differential.R:40-41` | Cells are treated as replicates; p-values are anti-conservative by orders of magnitude. Cannot use the 7-sample cm or 8-sample mouse designs at all. `08f_customer_export.R:80-89` already writes `pseudobulk_cpm.csv.gz`, so the aggregation exists but is not used for testing. | Rewrite 06b as pseudobulk: `AggregateExpression(group.by = c("sample","cell_type"))` → edgeR/DESeq2 with `~ condition` (or `~ timepoint + arm` for mouse). Require ≥ 2 replicates per arm, else emit a "descriptive only" notice. Use `SCRNA_CONDITION` for the design. | L |
| Proportions | Bootstrap within sample to `n_min` cells `09_bootstrap_proportions.R:51-60` | Bootstrap CI reflects sampling noise within one library, not replicate-to-replicate variance. No test between conditions. | Add propeller (speckle) or scCODA-style logit-linear test on sample-level proportions when ≥ 2 replicates per arm exist. Keep the bootstrap as a within-sample precision band. | M |
| Reproducibility | Seeds for UMAP (`config.R:165`), doublets, bootstrap; versions in report `08b_html_report.R:568-573` and export `08f_customer_export.R:206-208`; cache key `config.R:356` | Good. Gaps: `FindClusters` uses Seurat's default seed (0) implicitly; `RunHarmony` has no seed passed (`04_integrate.R:142-146`) and is stochastic across BLAS thread counts; no `renv.lock` or conda `environment.yml` export (`setup_env.sh` installs by name) [unverified whether a lock exists outside `pipeline/`]. | Pass `seed.use`/`set.seed` before Harmony; export `conda env export --no-builds` into the run dir at step 01. | S |
| Reporting | Cell fate, rho, doublet rate, versions, drift | Missing: realised QC thresholds per sample; integration mixing metric; cell-cycle phase; a "what this report cannot tell you" box (no condition inference without replicates). | Add the three numbers above once computed; add the caveat box to the customer template. | S |
| Tests | `tests/test_regressions.R`, 21 grep-based guards | No functional test. A 500-cell fixture run of 01–05 on `Samples/H1` would catch most breakage in minutes. | Add `tests/smoke.sh` that runs `run_pipeline.sh Samples/H1 01 02 03 04 05 --no-report` with `--max-cells 500` [flag unverified] and asserts the rds files exist. | M |

## 3. Missed entirely

Analyses the data supports that the pipeline never performs:

1. **Pseudobulk differential expression with biological replicates.** The 7-sample cm and
   8-sample mouse cohorts have replicates; the only DE path is 2-sample per-cell Wilcoxon.
2. **Compositional testing between conditions** (propeller / scCODA). Only within-sample
   bootstrap exists.
3. **Cell-cycle scoring.** Absent, despite a proliferating cluster in mouse and a
   differentiation time-course in cm.
4. **Pathway / gene-set enrichment on DE output.** Only `AddModuleScore` on hand-picked
   inflammatory sets `06b_differential.R:72`. No GSEA/ORA (fgsea, clusterProfiler) and no
   species-aware gene-set library. For cm the Yap 2023 module scores would be the natural
   first set.
5. **Integration quality metrics.** No LISI, kBET or batch silhouette; the only signal is a
   heuristic warning at `04_integrate.R:210`.
6. **QC drift across runs.** 08c benchmarks composition against a frozen baseline but not QC
   metrics (median genes, %mt, doublet rate, rho) per anchor sample, which is where a kit or
   flowcell change shows first.
7. **Trajectory on the cm time-course** exists only as a hand-run script
   (`projects/cm/trajectory.R`, monocle3); RNA velocity is not possible because the DNBelab
   quantifier output has no spliced/unspliced layers [unverified: check the C4 pipeline
   output for `velocity` matrices].
8. **Ligand–receptor analysis** exists only for bat wing (`projects/bat_wing/13_cellchat.R`);
   the mouse immune dataset is a textbook CellChat use case.
9. **Doublet rate sanity against loading.** Called doublet rate is reported but never
   compared with the expected rate from cells loaded, so a 2× discrepancy goes unnoticed.

## 4. Top 5 changes by impact / effort

| # | Change | Touches | Why first |
|---|---|---|---|
| 1 | Pseudobulk DE replacing per-cell Wilcoxon, replicate-aware design from `SCRNA_CONDITION` | `06b_differential.R`, new `DE` block in `config.R` | Every customer "what changed" question currently gets an invalid p-value. |
| 2 | Compositional test on sample-level proportions (propeller) | `09_bootstrap_proportions.R`, `08b_html_report.R` | Proportions are the headline figure in every report; they have no inferential test. |
| 3 | Stop regressing %mt for cm; per-layer HVGs; Leiden | `config.R` NORM/CLUSTER, `03_individual.R:56`, `04_integrate.R:59-62,75` | Three one-line changes that move the clustering to current defaults. |
| 4 | Adaptive QC option with realised thresholds written to the run dir and shown in the report | `01_load_qc.R:118-124`, `config.R:129-133`, `08b_html_report.R` | Removes the per-species guesswork that produced inert cut-offs for mouse. |
| 5 | SoupX as an optional step 00 inside `run_pipeline.sh`, no `forceAccept`, rho ceiling flag | `run_pipeline.sh`, `soupx_adjust.R:92`, `config.R` cache key | Ends the `_sx` parallel-cohort fork and makes ambient correction auditable. |

Cell-cycle scoring (S) and the Harmony seed (S) are worth doing in the same pass as #3.

## Citing

Methodological guidance was drawn from Scientific Agent Skills (Kassis et al., 2026,
arXiv:2609.00065); cite it if this review feeds a manuscript.
