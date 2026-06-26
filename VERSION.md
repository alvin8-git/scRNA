# Version History

## v0.9.0 — 2026-06-26

### Run-independent labels (frozen reference) + cross-run benchmark + report/PDF wiring

#### New features

- **Frozen bat reference** (`pipeline/build_reference.R`) — trains a SingleR classifier ONCE from a canonical annotated run (ES03 8-sample), marker-validated and population-complete (12 labels including CD8 T, RBC, Platelet, Neutrophil). Saved as a model bundle (`model` + `meta` + per-anchor `baseline`) under `Results/frozen_reference/` (gitignored, ~96 MB). Makes cell-type labels run-INDEPENDENT: they stop shifting with a run's sample mix.
- **Step 05r — reference transfer** (`05r_reference_transfer.R`) — `classifySingleR` (fine.tune) transfers the frozen labels onto a finished run, aligning genes to the model's gene set (zero-fills missing). Writes `annotation/reference_transfer_cells.csv.gz` (`cell_type_ref` per barcode) + `reference_transfer_composition.csv`. Additive: `cell_type_ref` sits alongside the de-novo `cell_type`. Gated on `REFERENCE_MODEL`; self-skips when unset, so non-bat runs are untouched.
- **Step 08c — cross-run benchmark** (`08c_benchmark_concordance.R`) — compares each anchor sample's frozen-reference composition in THIS run against the model's stored baseline; flags drift > `DRIFT_FLAG_PP` (5 pp). Emits a per-sample whole-blood signature (Neutrophil/RBC/Platelet %) as a sort readout. Writes `benchmark/concordance.csv`, `wholeblood_signature.csv`, `benchmark_report.md` (with collection-method + biological caveats). PASS: shared anchors reproduce within ~2.3–2.4 pp across the ES03-batch and ES17-batch integrations, resolving the original cross-run proportion discrepancy.
- **`apply_reference_labels()`** (`config.R`) — shared helper that swaps `cell_type` to `cell_type_ref` (keeping de-novo as `cell_type_denovo`) when the 05r CSV exists; records `options(scrna.label_source)`. Used by steps 06 and 09 so the PDFs show run-independent labels. Idempotent; no-op without a transfer.
- **HTML report (08b) frozen-reference view** — composition / cell-number bars and the UMAP default to `cell_type_ref` with a **Frozen reference vs De-novo toggle**; UMAP hover shows both calls; a new **Frozen-reference benchmark** nav section renders the concordance + whole-blood-signature tables. Falls back to de-novo when no transfer exists.
- **PDF label source** — `06_visualize.R` and `09_bootstrap_proportions.R` apply the frozen-reference labels; proportion plots print `Labels: frozen reference | de-novo` in the subtitle.
- **`SCRNA_RESULTS_DIR`** — env override pointing any step at a finished run dir, to re-render PDFs/reports without reconstructing `SCRNA_SAMPLE*`. Steps 06/07 recover the real sample list from the object / run-dir name.
- **Config knobs** — `REFERENCE_MODEL` (`SCRNA_REFERENCE_MODEL`), `ANCHOR_SAMPLES` (`SCRNA_ANCHORS`, default `Aksh1,ES332`), `DRIFT_FLAG_PP` (`SCRNA_DRIFT_PP`, 5), `REF_FINE_TUNE` (TRUE).
- **`run_pipeline.sh`** — 05r now runs immediately after step 05 (so 06/09 PDFs and the 08b HTML all get run-independent labels); 08c runs at the end. Both self-skip without a model and never fail the run.

#### Findings (bat whole blood)

- **ES17 cohort neutrophils recovered** — the de-novo annotation mislabeled the 24,737-cell neutrophil cluster as CD14+ Mono (cross-species SingleR vote flip at the human Neutrophil/Monocyte boundary). The frozen bat reference recovers them: ES18 63%, ES171 59%, ALL 35%. The original 0% was the artifact.
- **Both cohorts are unsorted whole blood (presort)** — high Neutrophil + RBC + Platelet together; confirmed not FACS- or Ficoll/PBMC-sorted. The earlier "assume postsort" guess for ES17 is superseded.
- **Aksh1 is a collection outlier** — terminal cardiac puncture (vs conscious vein draw) gives the lowest neutrophils and inflated lymphoid (B ~15%). Valid as a reproducibility anchor, not a biological baseline; use ES332 as the representative anchor.
- **Literature validation** (`docs/bat_neutrophil_literature.md`) — 25–63% neutrophils in captive bat whole blood is credible: neutrophil-dominant adult pteropodids, droplet-scRNA granulocyte dropout (so 0% is the artifact), captivity/handling stress amplification.

#### Files changed

- `pipeline/build_reference.R`, `pipeline/05r_reference_transfer.R`, `pipeline/08c_benchmark_concordance.R` — **new**
- `pipeline/config.R` — `apply_reference_labels()`; `REFERENCE_MODEL` / `ANCHOR_SAMPLES` / `DRIFT_FLAG_PP` / `REF_FINE_TUNE`; `SCRNA_RESULTS_DIR`
- `pipeline/08b_html_report.R` + `pipeline/report_template.Rmd` — ref-label default + toggle + benchmark section
- `pipeline/06_visualize.R`, `pipeline/09_bootstrap_proportions.R` — `apply_reference_labels()` + label-source subtitle; 06 takes its sample list from the object
- `pipeline/07_finalize_reports.R` — sample-list recovery under `SCRNA_RESULTS_DIR`
- `pipeline/run_pipeline.sh` — 05r after step 05; 08c at end
- `docs/frozen_reference_scope.md`, `docs/bat_neutrophil_literature.md` — **new**

---

## v0.8.0 — 2026-06-10

### Config Validation, Performance Hardening, Data Integrity Guard Rails, DX Improvements

#### New features

- **`pipeline/validate_config.R`** — pre-flight validator called by `run_pipeline.sh` before any step; catches broken `CLUSTER_CELLTYPE_MAP` color coverage, non-quoted integer keys, and missing `SAMPLE_PATHS` with actionable error messages; can also be run standalone via `Rscript pipeline/validate_config.R`; uses `commandArgs()` path resolution so it doesn't require a hard-coded repo path
- **`MARKERS$compute_integrated` flag** — gates `FindAllMarkers` in step 04 behind an explicit opt-in (default `FALSE`); previously ran unconditionally and added 20–30 minutes to every integration run; set to `TRUE` to write `integrated/integrated_cluster_markers.csv`
- **`SCRNA_BASE_DIR` env var** — overrides `BASE_DIR` without editing `config.R`; useful for running against data on NAS mounts or CI runners
- **`--help` / `-h` flag** on `run_pipeline.sh` — prints usage and exits

#### Performance

- **ScaleData HVG-only in step 05** — `ScaleData` now scales only the Highly Variable Genes (HVGs) instead of the full gene matrix; peak RAM drops from ~14 GB to ~1 GB on typical PBMC datasets; skips scaling entirely if `scale.data` is already present
- **Thread pinning** — `run_pipeline.sh` sets `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `BLAS_NUM_THREADS` to 1 to prevent CPU thrashing when `future` workers each spawn BLAS thread pools
- **Memory reclamation** — step 04 calls `rm(seu_list); gc()` after the merge step, freeing 2–4 GB per sample that was held redundantly alongside the merged object; step 05 adds `gc()` calls after the `saveRDS` loop to reclaim buffers before plot rendering
- **`FindAllMarkers` skip** — step 04 now defaults to not running the 20–30 min marker sweep; output is identical when `compute_integrated = FALSE`

#### Data integrity (7 guard rails — `05_annotate.R`, `01_load_qc.R`)

- **Contamination override now runs for both annotation paths** — previously only ran in the `else` (auto-annotation) branch, silently absorbing Neutrophil/HSPC/Basophil singletons into cluster-majority labels when `CLUSTER_CELLTYPE_MAP` was set
- **SingleR barcode-indexed assignment** — `singler_result$labels[colnames(merged)]` instead of positional assignment; guards against any worker ordering drift
- **`singler_delta` OOB guard** — `s[1] - s[2]` on a length-1 score vector now returns `0` instead of `NA`
- **`avg_expr` column guard** — skips cluster in subtype refinement loop if cluster name is absent from `avg_expr` columns (prevents `subscript out of bounds` crash on non-numeric cluster names)
- **Tie-breaking warning** — `which.max(table(singler_label_clean))` now emits `[TIE]` message when two labels tie at the same count so silent alphabetical wins are auditable
- **Empty factor level guard** — scType loop skips clusters with 0 cells instead of producing a 0-row `colSums` result that led to wrong/crash behaviour
- **`Read10X` modality guard** — fails fast with the names of available modalities if `"Gene Expression"` key is absent from a multi-modal matrix

#### DX

- **Timing messages** around long silent ops in steps 04 and 05 (`ScaleData`, `RunUMAP`, `FindNeighbors+FindClusters`, `SingleR`) with elapsed-seconds readout
- **PDF render warning** — `pdf_helpers.R` now emits a warning when `pdftools::pdf_render_page` fails for a specific page instead of silently returning `NULL` and producing a blank page
- Redundant `pdf_helpers.R` source call in step 07 removed

---

## v0.7.0 — 2026-05-21

### RBC/Platelet Detection, Monaco-Blind Cell Type Propagation, Composition Pagination

#### New features

- **scType Monaco-blind propagation** — `05_annotate.R`: after SingleR/manual assignment, any cluster scType labels as RBC, Platelet, Eosinophil, or Mast cell (types MonacoImmune lacks in its reference) is overridden in `cell_type`; emits `[scType override]` log message per cluster; uses `.sctype_monaco_blind` saved before the scType rm() cleanup block
- **RBC detection restored** — `MARKERS$RBC` (HBB, HBA1, HBA2, GYPA) added to scType `.gs_pos`; HBB/HBA1/HBA2/GYPA confirmed present in bat 10x feature matrix; previously cluster 10 (1,180 cells) was silently misidentified as B cell by MonacoImmune
- **Platelet detection restored** — `MARKERS$Platelet` now propagated to `cell_type` via Monaco-blind mechanism; previously cluster 14 (587 cells) was misidentified as CD14+ Mono
- **Eosinophil & Mast cell markers** — added to base `config.R` `MARKERS`: Eosinophil = (SIGLEC8, CCR3, EPX), Mast_cell = (TPSAB1, CPA3, MS4A2, KIT); added to scType `.gs_pos`; Monaco-blind propagation applies when clusters are identified by scType
- **HSPC added to scType** — `MARKERS$HSPC` (CD34, GATA2, AVP) now included in scType scoring; dual-detected by both scType and MonacoImmune "Progenitor cells" → SINGLER_NORM → HSPC
- **Composition plots paginated** — `06_visualize.R`: at most 3 samples per page; page titles include "(N/total)" suffix when >3 samples; `celltype_composition_combined.pdf` is a multi-page PDF merged via `.combine_pdfs()`
- **Overall_report multi-page composition** — `07_finalize_reports.R`: composition section loops over all pages of `celltype_composition_combined.pdf` (same pattern as split UMAP); each page titled "Cell Type Composition — Proportion & Count per Sample (N/total)"
- **`limitsize = FALSE` on individual feature plots** — `03_individual.R`: added to dotplot and feature plot ggsave calls; prevents crash when `ALL_MARKERS` grows beyond ~50 unique genes (feat_h formula can exceed ggplot2's 50-inch safety cap)

#### Bug fixes

- **scType propagation never ran** — `.cl_sctype` was cleaned up by `rm()` at line 267 before the propagation block at line 424; fixed by saving `.sctype_monaco_blind` before the rm() call
- **feat_h dimension crash** — adding Eosinophil/Mast cell markers pushed `ALL_MARKERS` from ~47 to ~54 genes; `feat_h = ceiling(54/4)*4 = 56` exceeded ggplot2's 50-inch cap; fixed with `limitsize = FALSE`

## v0.6.0 — 2026-05-08

### Bat 4-Sample Analysis, scType Annotation, Bootstrap & Rarefaction

#### New features

- **scType marker-based annotation** — `05_annotate.R` Part 1b runs scType alongside SingleR; genes weighted by inverse marker frequency (genes shared across fewer cell types score higher); scale.data HVG-safe filtering applied; outputs `sctype_labels_umap.pdf` and `singler_vs_sctype_comparison.csv`
- **Step 09 — Bootstrap proportion comparison** (`09_bootstrap_proportions.R`) — normalises cell type proportions across samples with different capture depths: 1000× bootstrap downsampling to smallest-n, multinomial 95% CI, pairwise chi-squared tile heatmap; outputs `bootstrap_proportions_report.pdf` + `bootstrap_summary.csv`
- **Step 10 — Rarefaction analysis** (`10_rarefaction.R`) — empirically determines minimum capture depth for stable proportion estimates; subsamples the largest sample at 11 depths (100–6000), 1000 draws each; plots CI width and RMSE vs n; fits theoretical a/√n curve; outputs `rarefaction_report.pdf` + `rarefaction_summary.csv`
- **Bat marker additions** — IDO1 added to Neutrophil panel (bat-specific granulocyte marker); CSF1R added to CD14_mono panel (monocyte colony-stimulating factor receptor)
- **4-sample bat run** — Sample6, Sample7, 10, ES03_newkit (20,709 total cells); cluster overrides verified by marker check:
  - Clusters 3, 11, 15 → CD8 T: CD3E 90.7%, NCAM1 0.5% — SingleR incorrectly called NK (CTL confirmed)
  - Clusters 13, 18 → DC: FLT3 72.8%/48%, HLA-DRA ~98%/80% — SingleR incorrectly called Monocyte
  - Cluster 17 → Neutrophil: S100A8/S100A9 100%, CSF3R 39.6% — SingleR called HSPC

#### Findings

- **NK → CD8 T correction**: what SingleR labels NK in bat whole blood are overwhelmingly CD3E+ (90.7%) cytotoxic T lymphocytes; true NK (NCAM1/CD56+) is only 0.0–0.7%; matches reference image label "Cytotoxic T lymphocytes / NK cells"
- **DC restoration**: clusters with FLT3 >48% and HLA-DR >80% are genuine DCs; SingleR absorbs them into Monocyte when CLUSTER_CELLTYPE_MAP is NULL
- **Neutrophil discrepancy vs human reference**: CSF3R 82% + CXCR2 50% confirm bat neutrophil identity; discrepancy with human-centric reference tools is expected — bat granulocyte signatures are absent from human SingleR references
- **Minimum capture depth**: rarefaction using ES03_newkit (n=6,520) as ground truth recommends ≥5,000 cells per sample for ±1% CI on all cell types including dominant populations (Monocyte)
- **Inter-bat composition differences**: all 6 pairwise chi-squared tests p << 0.001; Sample 10 is an outlier (CD4 T 43.6%, CD8 T 20.1%) consistent with published inter-individual immune variation in bats (Periasamy et al., 2019 PNAS)

#### Files changed

- `pipeline/config.R` — IDO1 to Neutrophil markers; CSF1R to CD14_mono markers; 4-sample CLUSTER_CELLTYPE_MAP (clusters 3/11/15/13/17/18)
- `pipeline/05_annotate.R` — Part 1b scType scoring block; marker sensitivity weighting; HVG-safe scale.data filtering
- `pipeline/09_bootstrap_proportions.R` — **new file**
- `pipeline/10_rarefaction.R` — **new file**

---

## v0.5.0 — 2026-05-07

### Bat Species Support + Configurable SingleR Reference

#### New features

- **`bat` / `human` species parameter** — `run_pipeline.sh` accepts a species keyword as first argument; exports `SCRNA_SPECIES`; `config.R` species override block activates bat-specific settings automatically.
- **Bat whole-blood overrides** — substitutes 3 marker genes absent from bat GTF (CST3 dropped, FCGR3A→FCGR2A, CEACAM8→CEACAM6); removes RBC and Neutrophil from `CONTAMINATION_TYPES` (both are expected in bat whole blood); bat-specific `SUBTYPE_MARKERS` for B cell and Monocyte lineages.
- **γδ T cell support** — `MARKERS$gamma_delta_T = c("TRDC", "TRGC1", "TRGC2")` added to bat override; `05_annotate.R` generates `feature_gamma_delta_T.pdf` when marker group is present.
- **Configurable SingleR reference** — `SINGLER_REF` parameter in `config.R` (default `"HumanPrimaryCellAtlas"`); bat override sets `"MonacoImmune"` which resolves 29 blood cell types including CD4 T, CD8 T, Treg, γδ T, classical/non-classical monocytes, pDC/mDC in a single pass.
- **Extended Monaco SINGLER_NORM** — 17 Monaco label mappings added (CD4+ T cells, CD8+ T cells, T regulatory cells, Vd2/non-Vd2 gd T cells, MAIT cells, NK cells, Plasmablasts, Classical/Intermediate/Non-classical monocytes, Plasmacytoid/Myeloid DC, Progenitor cells, Low-density neutrophils/basophils).
- **Higher bat clustering resolution** — bat block overrides `CLUSTER$resolutions`, `CLUSTER$default_res = 1.0`, `CLUSTER$compare_res` to improve resolution for whole-blood diversity.
- **Bat GTF merge** — `merge_gtf.py` merges `finalsort_gtf4MT_recalc.gtf` (30,011 genes) with TCR GTF; output `ESpe_merged_fullref.gtf` = 30,180 genes (169 new TCR loci added).

#### Fixes

- **`CLUSTER$resolutions` must include `default_res`** — bug where setting `CLUSTER$default_res` to a value not in `CLUSTER$resolutions` caused `RNA_snn_res.<X>` column to be absent, crashing `04_integrate.R` with `seurat_clusters not found`. Bat override block now sets all three CLUSTER resolution fields together.
- **Dynamic SingleR UMAP title** — title now reflects the actual reference used (`SingleR — HumanPrimaryCellAtlas` vs `SingleR — MonacoImmune`) instead of being hardcoded.

#### Files changed

- `pipeline/config.R` — `SINGLER_REF` default; bat block: `SINGLER_REF`, `CLUSTER$resolutions/default_res/compare_res`, `MARKERS$gamma_delta_T`, `ALL_MARKERS` refresh
- `pipeline/05_annotate.R` — dynamic ref loading; Monaco SINGLER_NORM entries; dynamic UMAP title; γδ T feature plot group
- `pipeline/run_pipeline.sh` — `bat`/`human` species keyword parsing; `SCRNA_SPECIES` export; species shown in pipeline header
- `/data/alvin/tmp/merge_gtf.py` — **new file** — GTF merge script; outputs `Bat/ESpe_merged_fullref.gtf`

---

## v0.4.0 — 2026-05-04

### Per-sample Cache, Contamination Detection, Step 08 Comparison Report

#### New features

- **Per-sample RDS cache** (`sample_cache/<name>/`) — steps 01–03 write `_filtered.rds`, `_singlets.rds`, `_seurat.rds` per sample on first run; subsequent integration runs with the same sample in a different combination load from cache, skipping re-computation. Invalidate by deleting `sample_cache/<name>/`.
- **Contamination / rare cell type detection** — `CONTAMINATION_TYPES` vector in `config.R` (Neutrophil, RBC, HSPC, Platelet, Basophil, Eosinophil, Mast cell); `05_annotate.R` overrides cluster majority-vote labels with per-cell SingleR labels for contamination types so they always appear on UMAP regardless of cluster size.
- **Extended SINGLER_NORM** — 30-entry lookup now covers the full granulocyte/erythroid lineage (Megakaryocyte, Neutrophil_-G-CSF, Pro-Myelocyte, Myelocyte, Basophils, Eosinophils, Mast_cells, Erythrocyte, BFU-E, CFU-E).
- **Contamination summary plot** — `annotation/contamination_summary.pdf` stacked bar of contamination cell counts per sample; added to `04-Annotation_report.pdf`.
- **Step 08 — Comparison report** (`08_comparison_report.R`) — standalone `Comparison_report.pdf` with cover page, sample quality & cell fate, doublet rates, cell type composition, integrated UMAP (all pages), and DE section. Run via `bash run_pipeline.sh <samples> 08`.

#### Fixes

- **Split-by-sample UMAP** — paginated at ≤2 panels per page; `pdftools::pdf_length()` loop in `07_finalize_reports.R` ensures all pages appear in `05-Integrated_report.pdf` and `Overall_report.pdf`.
- **Duplicate legend on split UMAP** — root cause: per-sample Seurat subsets have different cell type factor levels, producing separate legends. Fixed with `cowplot::get_legend()` on full merged object + `NoLegend()` on each per-sample panel assembled via `cowplot::plot_grid()`.
- **future.globals.maxSize** — raised from 4 GB to 16 GB (`future_mem_gb = 16L` in `config.R`) to prevent `future` worker memory limit errors with large merged objects (>40 K cells).

#### Files changed

- `pipeline/config.R` — `future_mem_gb = 16L`; `SAMPLE_CACHE_DIR`; `CONTAMINATION_TYPES`; new colors for Basophil/Eosinophil/Mast cell
- `pipeline/01_load_qc.R` — per-sample cache check/save
- `pipeline/02_doublets.R` — per-sample cache check/save
- `pipeline/03_individual.R` — per-sample cache check/save
- `pipeline/04_integrate.R` — cache fallback when individual RDS missing
- `pipeline/05_annotate.R` — contamination override loop; SINGLER_NORM extended; contamination summary plot
- `pipeline/06_visualize.R` — split UMAP: cowplot shared legend, 2-per-page pagination
- `pipeline/07_finalize_reports.R` — `pdftools::pdf_length()` loop for multi-page split UMAP
- `pipeline/08_comparison_report.R` — **new file**
- `pipeline/run_pipeline.sh` — step 08 added

---

## v0.3.0 — 2026-04-09

### Pipeline Refactor: Flexible Samples + Combined PDF Reports + Multi-core

#### New features
- **Dynamic sample folders** — pass 1 or more sample directories directly to `run_pipeline.sh`; sample names are auto-derived from folder names; `filter_matrix` subfolders are auto-detected
  ```bash
  bash run_pipeline.sh /data/alvin/Sample1 /data/alvin/Sample2
  bash run_pipeline.sh /data/alvin/Sample1                        # single-sample mode
  bash run_pipeline.sh /data/alvin/S1 /data/alvin/S2 /data/alvin/S3  # 3+ samples
  ```
- **Named results directory** — outputs written to `results_<SampleNames>/` (e.g. `results_Sample1andSample2/`) so runs with different sample sets never overwrite each other
- **N-sample integration** — Harmony integration now supports any number of samples (previously implicitly 2); single-sample mode skips Harmony automatically
- **Multi-core parallelism** — all heavy steps now use available CPU cores:
  - `scDblFinder` (step 02): `BiocParallel::MulticoreParam`
  - `FindAllMarkers`, `ScaleData` (steps 03, 04, 06): `future` multicore plan
  - `SingleR` (step 05): `BiocParallel::MulticoreParam`
  - Workers: `min(8, nCores - 2)` — configurable via `PARALLEL` in `config.R`

#### Combined PDF reports (step 07)
Five final reports written to `results_<samples>/`:

| File | Contents |
|------|----------|
| `report_qc.pdf` | QC violins and scatter plots per sample |
| `report_doublets.pdf` | Doublet UMAPs and score histograms per sample |
| `report_individual.pdf` | QC + HVG, elbow, PC heatmaps, UMAP, marker dot plots per sample |
| `report_annotation.pdf` | SingleR heatmap/UMAPs, marker feature plots, cell type UMAP |
| `report_integrated.pdf` | Harmony before/after, integrated UMAPs, dot plot, heatmap, violin, bar charts |

- Each figure has a **bold title banner** above the plot
- Small plots (histograms, elbow, bar charts) are **paired 2-per-portrait page**
- Large/wide plots get their own page at native aspect ratio
- Base-R plots (PC heatmaps, SingleR score heatmap) included via `qpdf`/`pdfunite`/`gs`

#### New step
- **Step 07** (`07_finalize_reports.R`) — merges partial per-step PDFs into the 5 final reports

#### Files changed
- `pipeline/config.R` — dynamic sample detection (`SCRNA_SAMPLE{N}` env vars), named `RESULTS_DIR`, `PARALLEL` config block, `save_report_pdf()` / `mark_small()` / `.combine_pdfs()` helpers
- `pipeline/01_load_qc.R` — dynamic sample loop, `plot_qc()` returns plot objects, saves `qc_report.pdf`
- `pipeline/02_doublets.R` — dynamic sample loop, BiocParallel, saves `doublets_report.pdf`
- `pipeline/03_individual.R` — dynamic sample loop, future multicore, saves `individual_report.pdf`
- `pipeline/04_integrate.R` — dynamic sample reads, single-sample path, future multicore, saves `integration_report.pdf`
- `pipeline/05_annotate.R` — BiocParallel for SingleR, saves `annotation_report.pdf`
- `pipeline/06_visualize.R` — future multicore, saves `visualization_report.pdf`
- `pipeline/07_finalize_reports.R` — **new file**
- `pipeline/run_pipeline.sh` — N-sample arg parsing, dynamic `LOG_DIR` / `RESULTS_DIR`, step 07 added

#### Dependencies added
```r
install.packages(c("future", "future.apply", "qpdf"))
BiocManager::install("BiocParallel")
```

---

## v0.2.0 — 2026-04-07

### Initial pipeline implementation

- Steps 01–06: Load & QC → Doublets → Individual → Integration → Annotation → Visualisation
- Seurat v5 + Harmony batch correction
- SingleR automated annotation (HumanPrimaryCellAtlasData)
- scDblFinder doublet detection
- Hardcoded for two samples: H1 and H2
- Individual PDF output per figure
