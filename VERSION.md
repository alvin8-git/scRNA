# Version History

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
