# Pipeline Steps Reference

The core pipeline is a sequence of R scripts in `pipeline/` — steps 01–10 plus 06b (differential expression). Each step reads `.rds` objects written by earlier steps and writes its own outputs to the results directory. You can re-run from any step without reprocessing earlier ones.

Bat-wing-specific downstream analysis (steps 11–14: wing DEGs, pathway enrichment, CellChat, trajectory) lives separately under `pipeline/projects/bat_wing/` and is only run in `bat_wing` species mode; it sources the core `config.R` and is not part of the standard PBMC/whole-blood run.

---

## Pre-flight — `validate_config.R`

**Runs automatically** before step 01 (called by `run_pipeline.sh`). Can also be run standalone:

```bash
Rscript pipeline/validate_config.R
```

**What it checks:**

| Check | Failure message |
|-------|----------------|
| All types in `CLUSTER_CELLTYPE_MAP` have a colour in `CELLTYPE_COLORS` | `CLUSTER_CELLTYPE_MAP types missing from CELLTYPE_COLORS: <list>` |
| `CLUSTER_CELLTYPE_MAP` keys are quoted integers (`"0"`, not `0`) | `CLUSTER_CELLTYPE_MAP keys must be quoted integers` |
| Every path in `SAMPLE_PATHS` exists on disk | `Sample path does not exist: <path>` |

Exits with code 1 on failure. All checks run before reporting, so you see every problem at once rather than one at a time. Uses a `commandArgs()` path fallback so it resolves `config.R` correctly whether sourced by `run_pipeline.sh` or run directly as a top-level `Rscript`.

---

## Step 01 — `01_load_qc.R`

**Inputs:** 10x Genomics matrix folders (`filter_matrix/` or `filtered_feature_bc_matrix/`)

**Outputs** — the checkpoint goes to `individual/<sample>/` (plus `sample_cache/<sample>/`), the plots
and tables to `qc/`:

| File | Description |
|------|-------------|
| `individual/<sample>/<sample>_filtered.rds` | Seurat object after QC filtering |
| `qc/<sample>_violin_qc.pdf` | Gene count, UMI count, and %MT distributions |
| `qc/<sample>_scatter_qc.pdf` | UMI vs genes and UMI vs %MT scatter plots |
| `qc/cell_fate.csv` | Per-**sample** attrition funnel (loaded → post-QC → post-doublet, % retained) |
| `qc/qc_summary_table.csv` | Per-sample cell/gene counts and medians |
| `qc/qc_report.pdf` | The QC plots in one paginated report |

**What it does:** Reads 10x matrices with `Read10X()`, creates Seurat objects, computes QC metrics (`nFeature_RNA`, `nCount_RNA`, `percent.mt`), filters cells outside thresholds in `QC`, and writes per-sample `.rds` objects. There is no pre-filter `.rds` — only the filtered object is saved.

**Config keys:** `QC`, `QC$min_features`, `QC$max_features`, `QC$min_counts`, `QC$max_counts`, `QC$max_percent_mt`

---

## Step 02 — `02_doublets.R`

**Inputs:** `individual/<sample>/<sample>_filtered.rds` (from step 01)

**Outputs** — the checkpoint goes to `individual/<sample>/`, the plots to `doublets/`:

| File | Description |
|------|-------------|
| `individual/<sample>/<sample>_singlets.rds` | Seurat object with `scDblFinder.score` and `scDblFinder.class` in metadata, doublets removed |
| `doublets/<sample>_doublet_umap.pdf` | UMAP coloured by doublet/singlet classification |
| `doublets/<sample>_doublet_score_hist.pdf` | Distribution of doublet probability scores |
| `doublets/doublets_report.pdf` | The doublet plots in one paginated report |

**What it does:** Runs `scDblFinder()` on each sample independently, adds doublet scores to cell metadata, and removes cells classified as doublets.

**Config keys:** `DOUBLET$doublet_rate` (NULL = auto)

---

## Step 03 — `03_individual.R`

**Inputs:** `individual/<sample>/<sample>_singlets.rds` (from step 02)

**Outputs** (per sample, in `results_*/individual/<sample>/` AND `sample_cache/<sample>/`):

| File | Description |
|------|-------------|
| `<sample>_seurat.rds` | Normalised, clustered Seurat object with UMAP |
| `<sample>_umap_cluster.pdf` | UMAP coloured by cluster |
| `<sample>_umap_sample.pdf` | UMAP coloured by sample |
| `<sample>_umap_markers.pdf` | Canonical marker feature plots on UMAP |
| `<sample>_dotplot_markers.pdf` | Canonical markers × clusters (dot plot) |
| `<sample>_elbow.pdf` | PCA elbow plot |
| `<sample>_hvg.pdf` | Highly variable genes plot |
| `<sample>_pc_heatmaps.pdf` | Per-PC loading heatmaps |
| `<sample>_cluster_markers.csv` | Full marker gene table (Wilcoxon, all clusters) |

`individual/individual_report.pdf` collects the per-sample plots into one paginated report.

**What it does:** Normalises (`LogNormalize`), finds HVGs, scales, runs PCA and UMAP, clusters at all resolutions in `CLUSTER$resolutions`, and finds cluster markers. Results are also written to `sample_cache/` so multi-sample runs do not re-process the same sample twice.

**Config keys:** `NORM`, `DIM`, `CLUSTER`, `PLOT`

---

## Step 04 — `04_integrate.R`

**Inputs:** All `sample_cache/<sample>/<sample>_seurat.rds` objects

**Outputs** (in `results_*/integrated/`):

| File | Description |
|------|-------------|
| `integrated_seurat.rds` | Merged + Harmony-corrected Seurat object |
| `harmony_before_after.pdf` | UMAP before and after Harmony, coloured by sample |
| `integrated_umap_cluster.pdf`, `integrated_umap_sample.pdf`, `integrated_umap_split_sample.pdf` | Integrated UMAP variants |
| `integration_report.pdf` | The integration plots in one paginated report |
| `cluster_resolution_comparison.pdf` | Side-by-side UMAPs at `compare_res` resolutions |

**What it does:** Merges all per-sample objects, runs Harmony batch correction (or direct UMAP/clustering for a single sample), and saves an integrated Seurat object ready for annotation.

**Config keys:** `HARMONY`, `CLUSTER$compare_res`, `DIM`, `MARKERS$compute_integrated`

**`MARKERS$compute_integrated`:** When `TRUE`, runs `FindAllMarkers` after integration and writes `integrated/integrated_cluster_markers.csv`. Defaults to `FALSE` — skips the sweep, which saves 20–30 minutes on typical runs. Enable only when you need the full per-cluster marker table.

**Single-sample behaviour:** Skips Harmony; runs UMAP and clustering directly on the individual PCA.

---

## Step 05 — `05_annotate.R`

**Inputs:** `integrated/integrated_seurat.rds`

**Outputs** — the annotated object goes to `integrated/`, everything else to `annotation/`:

| File | Description |
|------|-------------|
| `integrated/integrated_annotated.rds` | Seurat object with `cell_type` in metadata |
| `cluster_annotation_table.csv` | Per-cluster majority label (column `final_cell_type`) and SingleR majority |
| `consensus_annotation.csv`, `singler_vs_sctype_comparison.csv` | SingleR vs ScType agreement tables |
| `contamination_summary.pdf` | Contamination-type prevalence per sample |
| `singler_scores_heatmap.pdf` | Per-cell SingleR score heatmap |
| `singler_delta_umap.pdf` | Annotation confidence (delta score) on UMAP |
| `canonical_markers_dotplot.pdf` | Canonical markers × clusters (use to fill `CLUSTER_CELLTYPE_MAP`) |
| `celltype_umap.pdf` | UMAP coloured by `cell_type` |
| `tcell_subclusters_umap.pdf` | T cell sub-cluster UMAP (if `SUBCLUSTER$enabled`) |
| `tcell_subclusters_dotplot.pdf` | T cell sub-cluster marker dot plot |
| `tcell_subcluster_summary.csv` | Sub-cluster cell counts and parent mapping |

**What it does:** Runs SingleR against the configured reference, normalises raw labels via `SINGLER_NORM` (56-entry mapping), applies per-cell contamination-type overrides, optionally applies `CLUSTER_CELLTYPE_MAP`, refines coarse T/B/mono labels using `SUBTYPE_MARKERS`, and writes `cell_type` to cell metadata. Prints a copy-pasteable `CLUSTER_CELLTYPE_MAP` to `logs/05_annotate.log`.

**Memory note:** `ScaleData` in this step scales only the Highly Variable Genes (HVGs) identified by `FindVariableFeatures`, not the full gene matrix. This reduces peak RAM from ~14 GB to ~1 GB on typical PBMC datasets. If `scale.data` is already present in the loaded object, scaling is skipped entirely.

**Config keys:** `SINGLER_REF`, `CLUSTER_CELLTYPE_MAP`, `CONTAMINATION_TYPES`, `SUBTYPE_MARKERS`, `SUBCLUSTER`, `MARKERS`

---

## Step 06 — `06_visualize.R`

**Inputs:** `integrated/integrated_annotated.rds`

**Outputs** (in `results_*/integrated/`):

| File | Description |
|------|-------------|
| `umap_triptych.pdf` | Cluster / sample / cell-type UMAP side-by-side |
| `umap_split_by_sample.pdf` | Per-sample UMAP panels (2 per page) |
| `integrated_umap_*.pdf` | Individual UMAP variants |
| `feature_<marker-group>.pdf` | Per-marker-group feature plots |
| `celltype_counts_bar.pdf` | Absolute cell counts per cell type per sample |
| `integrated_dotplot.pdf` | Canonical markers × cell type |
| `integrated_heatmap.pdf` | Top 3 markers per cluster (heatmap) |
| `celltype_proportions_bar.pdf` | Stacked bar: cell type proportions per sample |
| `celltype_composition_combined.pdf` | All composition plots combined |
| `violin_key_markers.pdf` | Key lineage markers violin per cell type |
| `visualization_report.pdf` | All above in one paginated report |

**What it does:** Generates publication-quality figures in 8 plot sets. Plot Set 8 is a text-summary set covering sorting effects, minUMI threshold effects, and B-cell population analysis (bat whole-blood runs only).

**Config keys:** `PLOT`, `SAMPLE_COLORS`, `CELLTYPE_COLORS`, `MARKERS`

---

## Step 06b — `06b_differential.R`

**Inputs:** `integrated/integrated_annotated.rds`

**Outputs** (in `results_*/differential/`):

| File | Description |
|------|-------------|
| `DE_<celltype>.csv` | Differentially expressed genes per cell type between samples |
| `volcano_<celltype>.pdf` | Volcano plot per cell type |
| `DE_all_celltypes.csv`, `DE_summary.csv` | Combined DE table and per-cell-type hit counts |
| `module_score_<module>.pdf` | Gene module scores across samples |
| `differential_report.pdf` | The DE plots in one paginated report |

**What it does:** Runs `FindMarkers()` for each cell type between samples defined by `SCRNA_CONDITION`. Skips automatically for single-sample runs.

**Config keys:** `SCRNA_CONDITION` (env var)

---

## Step 07 — `07_finalize_reports.R`

**Inputs:** All per-step PDF reports in `qc/`, `doublets/`, `individual/`, `annotation/`, `integrated/`

**Outputs** (at the run directory root — `DIRS$reports` is `RESULTS_DIR` itself, not a `reports/` subfolder):

| File | Contents |
|------|---------|
| `01-QC_report.pdf` | QC violin, scatter, cell fate plots |
| `02-Doublet_report.pdf` | Doublet score distributions and UMAP |
| `03-Individual_report.pdf` | Per-sample UMAPs, elbow, markers |
| `04-Annotation_report.pdf` | SingleR scores, delta, canonical dot plot |
| `05-Integrated_report.pdf` | Harmony UMAPs, triptych, heatmap, dot plot |
| `Overall_report.pdf` | A4-normalised curated cross-stage summary |

**What it does:** Merges per-step PDFs into five numbered category reports. Builds `Overall_report.pdf` by rasterising each source page with `pdftools`, normalising to A4, adding bold title banners and interpretation captions. Falls back to a simple PDF merge if `pdftools` is unavailable.

**Config keys:** `PLOT_CAPTIONS` (maps figure name patterns to `Good: … | Bad: …` captions)

---

## Step 08 — `08_comparison_report.R`

**Inputs:** `integrated/integrated_annotated.rds`, output files from steps 01–06b

**Outputs** (at the run directory root — `DIRS$reports` is `RESULTS_DIR` itself, not a `reports/` subfolder):

| File | Description |
|------|-------------|
| `Comparison_report.pdf` | Standalone cross-sample comparison PDF |

**What it does:** Generates a focused comparison document covering sample quality, doublet rates, cell-type composition, integration quality, and DE results across samples. Independent of step 07.

---

## Step 09 — `09_bootstrap_proportions.R`

**Inputs:** `integrated/integrated_annotated.rds`

**Outputs** (at the run directory root):

| File | Description |
|------|-------------|
| `bootstrap_proportions_report.pdf` | Bootstrap-normalised proportions with 95% CI error bars |
| `bootstrap_summary.csv` | Per-sample per-cell-type observed proportion + multinomial CI + bootstrap mean/CI |

**What it does:** Bootstraps cell-type proportions (1,000 resamples, each sample downsampled to the smallest) to produce multinomial 95% confidence intervals. Also runs pairwise chi-squared tests to identify statistically significant composition differences; those results appear in the PDF, not in a separate CSV.

**Not run by `run_pipeline.sh`** — invoke it directly against a finished run:

```bash
SCRNA_RESULTS_DIR=Results/results_<samples>_filtered Rscript pipeline/09_bootstrap_proportions.R
```

---

## Step 10 — `10_rarefaction.R`

**Inputs:** `integrated/integrated_annotated.rds`

**Outputs** (at the run directory root):

| File | Description |
|------|-------------|
| `rarefaction_report.pdf` | Proportion stability vs cell count, with the fitted CI ~ a/√n curve per cell type |
| `rarefaction_summary.csv` | Per-cell-type CI width, RMSE vs ground truth, and minimum stable cell count |

**What it does:** Treats the largest sample (by cell count) as ground truth, subsamples at increasing depths with 1,000 draws each, and reports the minimum *n* at which each cell type's empirical CI comes within 5% of its asymptote.

**Not run by `run_pipeline.sh`** — invoke it directly against a finished run:

```bash
SCRNA_RESULTS_DIR=Results/results_<samples>_filtered Rscript pipeline/10_rarefaction.R
```

---

## Intermediate `.rds` objects

Steps 01–03 write each checkpoint twice: into `individual/<sample>/` for the run, and into
`sample_cache/<sample>/` so later runs can skip the work.

| Object | Location | Written by | Read by |
|--------|----------|-----------|---------|
| `<sample>_filtered.rds` | `individual/<sample>/` + `sample_cache/<sample>/` | Step 01 | Step 02 |
| `<sample>_singlets.rds` | `individual/<sample>/` + `sample_cache/<sample>/` | Step 02 | Step 03 |
| `<sample>_seurat.rds` | `individual/<sample>/` + `sample_cache/<sample>/` | Step 03 | Step 04 |
| `integrated_seurat.rds` | `integrated/` | Step 04 | Step 05 |
| `integrated_annotated.rds` | `integrated/` | Step 05 | Steps 05r, 06, 06b, 07, 08, 08b, 09, 10 |

---

## Related

- [Configuration Reference](reference-config.md) — all config.R parameters
- [Output Files Reference](reference-outputs.md) — full output file inventory
- [How to Re-run from a Specific Step](howto-rerun-steps.md)
- [How to Override Cell Type Annotations](howto-override-annotations.md)
