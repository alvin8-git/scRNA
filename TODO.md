# TODO

## Pending

- [ ] **validate_config.R: warning-only mode for SAMPLE_PATHS** — Check 3 in `pipeline/validate_config.R` calls `dir.exists()` for every sample path. If data lives on a NAS or external drive that isn't mounted yet, the validator fails with a confusing error before any pipeline work. Add a `--skip-paths` flag or convert Check 3 to a warning (not an error) so validation can still run on config-only checks before drives are attached. *Surfaced during 2026-06-10 plan-eng-review of T6.*

- [ ] **bat_wing species documentation** — `config.R` advertises `bat_wing` as a supported species alongside `human` and `bat`, but none of the docs (howto-bat-whole-blood.md, reference-config.md) cover it. Add a section to `docs/howto-bat-whole-blood.md` or a stub noting "experimental" if the mode is not yet ready for external use.

### Office-hours architecture findings (2026-06-11) — ranked for implementation

Severity = blast radius if left unfixed, not effort. P0 (cache collision), P1, and the two P2 fixes are done (see Done).

**All resolved 2026-06-11 (see Done):** the BLAS item was dropped on inspection (every heavy step fans out `PARALLEL$workers` and each worker uses BLAS, so the global 1-thread pin is correct); it was replaced by per-phase worker tuning, now implemented. P3 (droplevels, doc drift, extract 11–14, decompose config + SINGLER_NORM, add.cell.ids) and the annotation-loop enhancement are all done and verified.

Residual (optional, deferred): a deeper split of `05_annotate.R` into annotate-core vs annotate-refine modules (SINGLER_NORM already extracted; the rest is logic-heavy and lower value); and the human-PBMC Azimuth/`MapQuery` fast path for the annotation loop (tracked under Optional Extensions).

---

## Done

- [x] **Annotation-loop automation (office-hours enhancement, 2026-06-11)** — new consensus block in `05_annotate.R` fuses the two per-cluster calls already computed (`cluster_singler$majority_singler` + scType `.cl_sctype$sctype`) into an AUTO/REVIEW decision gated by mean `singler_delta` (≥0.10). Writes `consensus_annotation.csv` and prints a pre-filled `CLUSTER_CELLTYPE_MAP` with `REVIEW` markers on ambiguous clusters — collapses "hand-build the whole map" into "resolve the few flagged clusters." Advisory only (does not change the assignment logic). Decision logic unit-tested.
- [x] **Per-phase worker tuning (office-hours P2, 2026-06-11)** — `config.R` `PARALLEL` now has `workers` (per-sample, 8 GB/worker) and `merge_workers` (merged-object phase 04/05/06/06b, 16 GB/worker → fewer workers), both RAM-governed. Verified: 8 per-sample / 5 merged on this box. Steps 04/05/06/06b switched to `merge_workers`/`merge_mem_gb`.
- [x] **droplevels() after T-cell subset (office-hours P3, 2026-06-11)** — `05_annotate.R` drops empty factor levels on the subcluster object so DotPlot legends/axes don't show stale non-T-cell categories. (Step 06 heatmap subset needs none — every cell_type group is sampled.)
- [x] **Doc drift "10 steps" → reality (office-hours P3, 2026-06-11)** — `reference-pipeline-steps.md` + README now describe the core steps (01–10, 06b) and point to the bat-wing extension (11–14) under `pipeline/projects/bat_wing/`.
- [x] **Extracted steps 11–14 to project extension (office-hours P3, 2026-06-11)** — `git mv` to `pipeline/projects/bat_wing/`; each script's `config.R` source repointed two dirs up; `run_pipeline.sh` STEP_SCRIPT paths + logfile basenaming updated. Off the shared trunk; not run in standard PBMC/whole-blood mode. Parse + path-resolution verified (not yet run end-to-end on a bat_wing dataset).
- [x] **Decomposed config.R + SINGLER_NORM (office-hours P3, 2026-06-11)** — `config.R` 783→609 lines: bat/bat_wing override block extracted to `config_species_bat.R` (181 lines), sourced conditionally. Proven behavior-preserving (all 17 key objects byte-identical for human/bat/bat_wing). `05_annotate.R` 756→696: `SINGLER_NORM` (49 entries) extracted to `05_annotate_singler_norm.R`, verified identical.
- [x] **`merge(add.cell.ids = SAMPLE_NAMES)` (office-hours P3, 2026-06-11)** — `04_integrate.R` prefixes merged barcodes with the sample name (self-documenting; provenance still in `$sample`). NOTE: changes post-merge barcode format — smoke-test on a real 2-sample run before relying on it.
- [x] **SingleR all-NA guard (office-hours P1, 2026-06-11)** — `05_annotate.R` now hard-stops with a diagnostic count if the barcode-indexed `singler_result$labels[colnames(merged)]` reindex yields any NA. `$labels` is never legitimately NA (unlike `$pruned.labels`), so NA means the row/colname alignment broke — e.g. a SingleR/Seurat upgrade returning positionally-named labels would otherwise silently turn every cell into "Unassigned". Tested (clean labels pass, misaligned caught).
- [x] **Per-step cache scoping (office-hours P2, 2026-06-11)** — `cache_hash(nm, step)` with cumulative `.cache_params`: 01 = QC+species, 02 += DOUBLET, 03 += NORM/DIM/CLUSTER. A clustering-resolution change no longer busts the QC (01) or doublet (02) caches; a QC change still cascades to all three. Verified by isolation test. Call sites in 01–03 pass their step.
- [x] **RAM governor (office-hours P2, 2026-06-11)** — `config.R` caps `PARALLEL$workers` so `workers × future_mem_gb` leaves ~20% of `/proc/meminfo` MemAvailable as headroom; prevents OOM in the fan-out phases (every step spawns workers, each holding an object copy). Degrades to the CPU-only cap when meminfo is unreadable; logs a line when it reduces. (Supersedes the stale "raised to 16 GB" note — per-worker budget is 8 GB, now bounded by available RAM.)
- [x] **Cache key collision fixed (office-hours P0, 2026-06-11)** — `cache_hash()` was keyed on sample name + config only, so two experiments sharing a folder name (e.g. `PBMC/`) under identical config silently served each other's cached `.rds` objects. Now keyed on the resolved absolute input path + a fingerprint (size+mtime) of the 10x matrix files (`config.R`); this also busts the cache when the source matrix changes under unchanged config (closes the content-blindness gap too). Signature changed to `cache_hash(nm)`; call sites updated in steps 01–03. Parse-checked + collision smoke-tested (same name / different path → distinct hashes). First run after this change recomputes each sample once (key format changed).
- [x] **Samba share** — installed samba 4.15.13; map from Windows: `\\192.168.1.168\alvin`
- [x] **Pipeline: dynamic sample folders** — `run_pipeline.sh` accepts any number of folder paths; names auto-derived; `SINGLE_SAMPLE` mode skips Harmony
- [x] **Pipeline: raw_matrix support** — `raw_matrix` / `raw_feature_bc_matrix` recognised alongside `filter_matrix`; results folder appends `_filtered` or `_raw` to differentiate runs
- [x] **Pipeline: combined PDF reports** — 5 numbered reports (`01-QC_report.pdf` … `05-Integrated_report.pdf`) + `Overall_report.pdf` (curated cross-stage summary); each figure has bold title banner + Good/Bad caption
- [x] **Pipeline: report naming** — reports renamed `01-QC_report.pdf` … `05-Integrated_report.pdf`; `Overall_report.pdf` combines QC + doublets + SingleR heatmap + individual HVG/markers + all integrated plots
- [x] **Pipeline: remove hardcoded H1/H2** — all scripts loop over `SAMPLE_NAMES` dynamically
- [x] **Step 07** — `07_finalize_reports.R` merges partial PDFs into category reports + `Overall_report.pdf`
- [x] **conda auto-activation** — `run_pipeline.sh` activates `scrna_seurat` automatically if not already active
- [x] **Relative path support** — `run_pipeline.sh` accepts relative paths (e.g. `H1/raw_matrix/`) as sample inputs
- [x] **Per-plot captions** — all 5 report PDFs include a description + Good/Bad interpretation guide under each figure
- [x] **ReportGuide.md** — detailed guide with plot descriptions and embedded images from H1+H2 run
- [x] **Expanded README** — QC threshold table, cell type detection explanation, marker contamination guide, annotation workflow, species adaptation table
- [x] **Extended marker list** — added Neutrophil (FCGR3B, CSF3R, CXCR2, CEACAM8), Treg (FOXP3, IL2RA, CTLA4), Plasma (MZB1, JCHAIN, SDC1), RBC (HBB, HBA1, HBA2, GYPA), HSPC (CD34, GATA2, AVP)
- [x] **Feature plot pagination** — all feature plots paginated 3 cols × 2 rows per page; dedicated Neutrophil and RBC plots added; "All Canonical Markers" replaced by "Remaining Markers (Treg/Plasma/HSPC)" to avoid duplicating genes already shown elsewhere
- [x] **Unassigned → grey** — Unassigned/Unknown forced to `#AAAAAA` in SingleR UMAP, cell type UMAP, and heatmap group bar
- [x] **Per-cell SingleR labels** — when no `CLUSTER_CELLTYPE_MAP` is set, per-cell labels used (not majority-vote per cluster) so rare types (Neutrophils, Unassigned) appear in the cell type UMAP
- [x] **Colour fixes** — Neutrophil → `#E377C2` (pink, distinct from NK teal); HSPC → `#8C564B` (brown, distinct from Monocyte salmon)
- [x] **Dot plot rotation** — gene labels rotated 90° CCW (`angle=90, hjust=1, vjust=0.5`) in integrated dot plot and individual Top 5 Markers plot (size 7)
- [x] **UMAP triptych labels** — cell type panel uses `label.size=2.5` to reduce clutter
- [x] **Cell composition plot** — stacked % bar (one bar per sample) with integer % labels; cell type name shown in large segments (≥7%); compact horizontal legend strip between panels; overall title removed from combined figure
- [x] **QC violin x-axis** — labels at 0° (horizontal)
- [x] **Doublet UMAP title** — includes count and % e.g. `H1 - Doublet UMAP (31 in 400 cells (7.8%))`
- [x] **Violin key lineage markers** — each sub-plot titled with cell type association e.g. `PPBP (Platelet)`, `CD3D (T cells)`
- [x] **Cell fate tracking** — `cell_fate.csv` saved to `results/qc/`; accounts for every GEM barcode: load losses (min.features < 200) → QC losses → doublet removal → final retained cells + %
- [x] **SingleR label normalisation** — `SINGLER_NORM` lookup in `05_annotate.R` maps raw SingleR labels to canonical pipeline names (Platelets→Platelet, Neutrophils→Neutrophil, NK_cell→NK, B_cell→B cell, Endothelial_cells→Endothelial, etc.); prevents split colours in UMAP
- [x] **Sub-type refinement** — `SUBTYPE_MARKERS` in `config.R`; after SingleR majority vote, generic labels (CD4 T, B cell, Monocyte) are refined to sub-types (naive/memory/effector, CD14+/FCGR3A+) using average expression of lineage markers; resolves cluster collapse problem
- [x] **Non-immune stromal labels** — Endothelial, Epithelial, Fibroblast, Smooth Muscle added to `CELLTYPE_COLORS` (grey scale); appear as named contamination labels instead of "Unknown"
- [x] **Auto-annotation implemented** — `05_annotate.R` uses SingleR majority per cluster when `CLUSTER_CELLTYPE_MAP = NULL`; prints copy-pasteable suggested map to `logs/05_annotate.log`
- [x] **Fallback logic** — partial `CLUSTER_CELLTYPE_MAP` entries override specific clusters; unmapped clusters fall back to SingleR automatically
- [x] **DemoScRNA independently verified** — all 13 clusters annotated by marker evidence; 3 errors in H1/H2 map caught; 4 SingleR mislabels corrected (clusters 7→Neutrophil, 11→DC, 12→Platelet, 9→γδ T)
- [x] **New cell type labels & colors** — added `CD4 T (effector)`, `γδ T`, `B cell (naive)`, `B cell (memory)` to `CELLTYPE_COLORS`
- [x] **Step 06b — Differential Expression** — `06b_differential.R` runs `FindMarkers` (Wilcoxon) per cell type between 2 samples; saves per-cell-type CSVs, volcano plots, inflammatory/interferon/T-activation module score violin plots, `DE_summary.csv`, `differential_report.pdf`; auto-skipped for single-sample runs
- [x] **Per-sample RDS cache** — steps 01–03 save results to `sample_cache/<name>/`; subsequent integration runs with the same sample skip reprocessing; cache invalidation by deleting `sample_cache/<name>/`
- [x] **future.globals.maxSize raised to 16 GB** — `config.R` `future_mem_gb = 16L`; prevents future worker memory limit errors with large merged objects (>40 K cells)
- [x] **Contamination cell type override** — `CONTAMINATION_TYPES` vector in `config.R`; after majority-vote labelling, per-cell SingleR labels for Neutrophil/RBC/HSPC/Platelet/Basophil/Eosinophil/Mast cell override cluster label so all contamination cells appear on UMAP
- [x] **SINGLER_NORM extended to 30 entries** — covers full granulocyte/erythroid lineage (Megakaryocyte, Neutrophil_-G-CSF, Pro-Myelocyte, Myelocyte, Basophils, Eosinophils, Mast_cells, Erythrocyte, BFU-E, CFU-E)
- [x] **Split-by-sample UMAP pagination** — at most 2 panels per page; 3+ samples produce multiple pages; `pdftools::pdf_length()` loop in `07_finalize_reports.R` includes all pages in reports
- [x] **Shared legend on split UMAP** — `cowplot::get_legend()` from full merged object; `NoLegend()` on each per-sample panel; eliminates duplicate legends caused by differing factor levels between subsets
- [x] **Step 08 — Comparison Report** — `08_comparison_report.R` generates `Comparison_report.pdf` with cover, sample quality, doublet rates, composition, integrated UMAP, and DE sections; runs via `bash run_pipeline.sh <samples> 08`

---

## RBC / Cell Type Detection Fixes (v0.7.0)

- [x] **RBC detection fixed** — `MARKERS$RBC` added to scType `.gs_pos`; cluster 10 (1,180 cells) correctly relabelled from "B cell" to "RBC" via Monaco-blind propagation
- [x] **Platelet detection fixed** — cluster 14 (587 cells) correctly relabelled from "CD14+ Mono" to "Platelet" via Monaco-blind propagation
- [x] **Monaco-blind propagation** — generalised override for RBC, Platelet, Eosinophil, Mast cell: scType label wins when MonacoImmune cannot detect the type; `.sctype_monaco_blind` saved before scType cleanup `rm()` call (bug fix)
- [x] **Eosinophil & Mast cell markers added** — to base `config.R` MARKERS and scType `.gs_pos`; all marker genes confirmed present in bat 10x feature matrix
- [x] **HSPC added to scType** — `MARKERS$HSPC` now in scType gene sets; dual-detection via scType + MonacoImmune "Progenitor cells"
- [x] **Composition plots paginated** — `06_visualize.R`: ≤3 samples per page; "(N/total)" title suffix; `celltype_composition_combined.pdf` is multi-page
- [x] **Overall_report composition pagination** — `07_finalize_reports.R` loops over all composition pages (same pattern as split UMAP)
- [x] **`limitsize = FALSE` fix** — `03_individual.R`: prevents ggsave crash when `ALL_MARKERS` exceeds ~50 unique genes
- [x] **Cross-species marker accuracy audit** — Tier 1 (CD3, CD14, CD79A): reliable; Tier 2 (NKG7, HBB): moderate; Tier 3 (CCR7, SIGLEC8): less reliable; NCAM1/CEACAM8/CST3/FCGR3A confirmed non-functional in bat

---

## 4-Sample Bat Analysis (Sample6, Sample7, 10, ES03_newkit)

- [x] **scType annotation** — marker-based scoring added to `05_annotate.R` Part 1b alongside SingleR; inverse marker-frequency weighting; HVG-safe scale.data filtering; outputs `sctype_labels_umap.pdf` + `singler_vs_sctype_comparison.csv`
- [x] **Bat marker additions** — IDO1 added to Neutrophil panel (bat-specific granulocyte marker confirmed at 43% expression); CSF1R added to CD14_mono panel
- [x] **NK → CD8 T cluster override** — clusters 3, 11, 15 confirmed CD3E 90.7%, NCAM1 0.5% → relabelled CD8 T; true NK is 0–0.7% as expected
- [x] **DC cluster restoration** — clusters 13, 18 confirmed FLT3 72.8%/48%, HLA-DRA ~98%/80% → relabelled DC; SingleR incorrectly called Monocyte
- [x] **Neutrophil cluster fix** — cluster 17 confirmed S100A8/A9 100%, CSF3R 39.6% → relabelled Neutrophil; SingleR called HSPC
- [x] **Bootstrap proportions** — `09_bootstrap_proportions.R` added; 1000× downsampling to n_min=3604 (Sample6); all 6 pairwise chi-squared p << 0.001; Sample10 outlier (CD4 T 43.6%, CD8 T 20.1%)
- [x] **Rarefaction analysis** — `10_rarefaction.R` added; ES03_newkit (n=6520) as ground truth; recommendation ≥5000 cells per sample for ±1% CI
- [x] **Cell proportion comparison vs reference image** — CellProportions_4samples.png compared; NK→CD8 T correction resolves CTL/NK discrepancy; neutrophil gap is bat-biology (human references lack bat granulocyte signatures); DC now restored

---

## Bat (*Eonycteris spelaea*) Support

- [x] **Bat species parameter** — `run_pipeline.sh` accepts `bat` or `human` keyword; exports `SCRNA_SPECIES`; `config.R` reads it and applies whole-blood overrides automatically
- [x] **Bat whole-blood overrides in config.R** — substitutes 3 absent marker genes (CST3→drop, FCGR3A→FCGR2A, CEACAM8→CEACAM6); removes RBC/Neutrophil from `CONTAMINATION_TYPES` (expected in whole blood); bat-specific B cell and Monocyte `SUBTYPE_MARKERS`
- [x] **γδ T cell markers** — `MARKERS$gamma_delta_T = c("TRDC", "TRGC1", "TRGC2")` added to bat override block; TRDC/TRGC1/TRGC2 present in merged bat GTF
- [x] **MonacoImmuneData reference for bat** — `SINGLER_REF <- "MonacoImmune"` set in bat block; resolves CD4 T, CD8 T, Treg, γδ T, classical/non-classical monocytes, pDC/mDC in a single pass (29 blood cell types vs 4 from HumanPrimaryCell)
- [x] **Higher clustering resolution for bat** — `CLUSTER$resolutions`, `CLUSTER$default_res = 1.0`, `CLUSTER$compare_res` all overridden in bat block; ensures sufficient clusters for whole-blood diversity
- [x] **Bat reference GTF** — merged old `finalsort_gtf4MT_recalc.gtf` (30,011 genes) with new TCR GTF `ESpe_Peaks2UTRed_genome_tcellgenes_ZF_Dec2024_v2.gtf`; result: `ESpe_merged_fullref.gtf` (30,180 genes = 30,011 + 169 new TCR loci); 38 TCR genes already present in old reference
- [x] **GTF gene count clarification** — new GTF has only 207 `gene`-level feature lines (TCR genes) but 30,156 unique `gene_name` values across transcript/exon/CDS features; gene count does not affect GEM (STAR counts via exon features)
- [x] **Bat pipeline run** — ES03_newkit + ES12_newkit (7,725 cells) analysed successfully; 4 initial cell types at res=0.5 (CD14+ Mono, CD4 T (memory), NK, B cell (memory)); re-running at res=1.0 with Monaco reference for improved type resolution

## Analysis Improvements

- [x] **Split CD4/CD8 T cells** — `05_annotate.R` PART 5: `FindSubCluster()` on T cell clusters (matching `SUBCLUSTER$t_patterns`); outputs `tcell_subclusters_umap.pdf`, `tcell_subclusters_dotplot.pdf`, `tcell_subcluster_summary.csv`. At res=1.2 gave 18 sub-clusters (over-split); suggest lowering `SUBCLUSTER$resolution` to 0.6–0.8
- [x] **Monocyte sub-typing** — `05_annotate.R` PART 4: per-cluster CD14 vs FCGR3A expression %; found FCGR3A monocytes dominant (~17.9%); CD14 monocytes at ~7%; both subtypes resolved
- [x] **Dendritic cells** — `05_annotate.R` PART 4: FCER1A/CLEC9A check; DC essentially absent (<0.3%); normal for small/healthy PBMC samples
- [x] **Harmony mixing quality** — `04_integrate.R`: per-cluster max sample % calculated; WARNING emitted if any cluster >80% from one sample with theta tuning suggestion
- [x] **Adjust clustering resolution** — `04_integrate.R`: `CLUSTER$compare_res = c(0.5, 0.6, 0.8)` produces side-by-side UMAP comparison saved as `cluster_resolution_comparison.pdf`
- [x] **Gene set scoring** — implemented in `06b_differential.R`: `AddModuleScore` for Inflammatory, Interferon, and T_Activation gene sets; violin+boxplot per cell type per sample

---

## QC Review

- [x] **H1 QC thresholds confirmed OK** — H1 started with 386 GEMs; 0 cells failed QC; MT% median 0.04%; H1 is a genuinely small library, not over-filtered
- [x] **Doublet rates reviewed** — H1: 7.0% (27/386), H2: 10.8% (103/951); both within expected range for 10x Genomics
- [x] **cell_fate.csv reviewed** — all removal counts confirmed: H1 retained 93%, H2 87.7%, DemoScRNA 88.2%

---

## Differential Expression

- [x] **Step 06b implemented** — `FindMarkers` per cell type (≥20 cells per group); volcano plots; DE CSVs; module scores
- [x] **Inflammatory signatures** — Inflammatory, Interferon, T_Activation gene sets scored and plotted per cell type × sample

---

## Figures for Publication / Report

- [ ] Finalise cell type labels (depends on manual annotation above)
- [ ] Export `umap_triptych.pdf` as main figure
- [ ] Export `integrated_dotplot.pdf` as supplementary marker validation figure
- [ ] Export `celltype_composition_combined.pdf` for sample comparison
- [ ] Consider PNG export for PowerPoint:
  ```bash
  conda activate scrna_seurat
  pdftoppm -r 150 results_*/integrated/umap_triptych.pdf umap_triptych
  ```

---

## Optional Extensions

- [ ] **Trajectory analysis** — run PAGA or Monocle3 on monocyte/DC lineage
- [ ] **Gene set scoring** — score cells for interferon response, exhaustion, activation signatures
- [ ] **Cell-cell communication** — run CellChat or NicheNet on integrated object
- [ ] **Reference mapping** — map to Azimuth PBMC reference for higher-resolution annotation
- [ ] **Batch effect validation** — compute LISI scores to quantify integration quality

---

## Environment

- [x] Install `qpdf` R package for PDF merging: `install.packages("qpdf")`
- [x] Add `scDblFinder` to `setup_env.sh` for fresh environment installs:
  ```bash
  mamba install -n scrna_seurat -c conda-forge -c bioconductor bioconductor-scdblfinder -y
  ```
- [x] Remove DoubletFinder GitHub install from `setup_env.sh` (replaced by scDblFinder)
- [x] Install `r-pdftools` — required for A4 Overall_report rasterization: `mamba install -n scrna_seurat -c conda-forge r-pdftools -y`
- [x] Install `r-magick` — required by `cowplot::draw_image()` for A4 page composition: `mamba install -n scrna_seurat -c conda-forge r-magick -y`
