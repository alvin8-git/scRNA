# Version History

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
