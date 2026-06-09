# Pipeline Steps Reference

The pipeline consists of 10 R scripts in `pipeline/`. Each step reads `.rds` objects written by earlier steps and writes its own outputs to the results directory. You can re-run from any step without reprocessing earlier ones.

---

## Step 01 — `01_load_qc.R`

**Inputs:** 10x Genomics matrix folders (`filter_matrix/` or `filtered_feature_bc_matrix/`)

**Outputs** (per sample, in `results_*/qc/<sample>/`):

| File | Description |
|------|-------------|
| `<sample>_seurat_raw.rds` | Seurat object before QC filtering |
| `<sample>_seurat.rds` | Seurat object after QC filtering |
| `qc_violin.pdf` | Gene count, UMI count, and %MT distributions |
| `qc_scatter.pdf` | UMI vs genes and UMI vs %MT scatter plots |
| `cell_fate.csv` | Per-cell QC pass/fail with reason |

**What it does:** Reads 10x matrices with `Read10X()`, creates Seurat objects, computes QC metrics (`nFeature_RNA`, `nCount_RNA`, `percent.mt`), filters cells outside thresholds in `QC`, and writes per-sample `.rds` objects.

**Config keys:** `QC`, `QC$min_features`, `QC$max_features`, `QC$min_counts`, `QC$max_counts`, `QC$max_percent_mt`

---

## Step 02 — `02_doublets.R`

**Inputs:** `qc/<sample>/<sample>_seurat.rds` (from step 01)

**Outputs** (per sample, in `results_*/doublets/<sample>/`):

| File | Description |
|------|-------------|
| `<sample>_seurat.rds` | Seurat object with `scDblFinder.score` and `scDblFinder.class` in metadata |
| `doublet_umap.pdf` | UMAP coloured by doublet/singlet classification |
| `doublet_score_hist.pdf` | Distribution of doublet probability scores |

**What it does:** Runs `scDblFinder()` on each sample independently, adds doublet scores to cell metadata, and removes cells classified as doublets.

**Config keys:** `DOUBLET$db_rate` (NULL = auto)

---

## Step 03 — `03_individual.R`

**Inputs:** `doublets/<sample>/<sample>_seurat.rds` (from step 02)

**Outputs** (per sample, in `results_*/individual/<sample>/` AND `sample_cache/<sample>/`):

| File | Description |
|------|-------------|
| `<sample>_seurat.rds` | Normalised, clustered Seurat object with UMAP |
| `<sample>_umap_clusters.pdf` | UMAP coloured by cluster |
| `<sample>_umap_celltype.pdf` | UMAP coloured by cell count per cluster |
| `<sample>_top_markers.pdf` | Top 5 DE markers per cluster (dot plot) |
| `<sample>_elbow.pdf` | PCA elbow plot |
| `<sample>_variable_genes.pdf` | Highly variable genes plot |
| `<sample>_markers.csv` | Full marker gene table (Wilcoxon, all clusters) |

**What it does:** Normalises (`LogNormalize`), finds HVGs, scales, runs PCA and UMAP, clusters at all resolutions in `CLUSTER$resolutions`, and finds cluster markers. Results are also written to `sample_cache/` so multi-sample runs do not re-process the same sample twice.

**Config keys:** `NORM`, `DIMS`, `CLUSTER`, `PLOT`

---

## Step 04 — `04_integrate.R`

**Inputs:** All `sample_cache/<sample>/<sample>_seurat.rds` objects

**Outputs** (in `results_*/integrated/`):

| File | Description |
|------|-------------|
| `integrated_seurat.rds` | Merged + Harmony-corrected Seurat object |
| `integration_umap_before.pdf` | UMAP before Harmony (coloured by sample) |
| `integration_umap_after.pdf` | UMAP after Harmony (coloured by sample) |
| `cluster_resolution_comparison.pdf` | Side-by-side UMAPs at `compare_res` resolutions |

**What it does:** Merges all per-sample objects, runs Harmony batch correction (or direct UMAP/clustering for a single sample), and saves an integrated Seurat object ready for annotation.

**Config keys:** `HARMONY`, `CLUSTER$compare_res`, `DIMS`

**Single-sample behaviour:** Skips Harmony; runs UMAP and clustering directly on the individual PCA.

---

## Step 05 — `05_annotate.R`

**Inputs:** `integrated/integrated_seurat.rds`

**Outputs** (in `results_*/annotation/`):

| File | Description |
|------|-------------|
| `integrated_annotated.rds` | Seurat object with `final_cell_type` in metadata |
| `singler_scores_heatmap.pdf` | Per-cell SingleR score heatmap |
| `singler_delta_umap.pdf` | Annotation confidence (delta score) on UMAP |
| `canonical_markers_dotplot.pdf` | Canonical markers × clusters (use to fill `CLUSTER_CELLTYPE_MAP`) |
| `annotation_umap.pdf` | UMAP coloured by `final_cell_type` |
| `tcell_subclusters_umap.pdf` | T cell sub-cluster UMAP (if `SUBCLUSTER$enabled`) |
| `tcell_subclusters_dotplot.pdf` | T cell sub-cluster marker dot plot |
| `tcell_subcluster_summary.csv` | Sub-cluster cell counts and parent mapping |

**What it does:** Runs SingleR against the configured reference, normalises raw labels via `SINGLER_NORM` (30-entry mapping), applies per-cell contamination-type overrides, optionally applies `CLUSTER_CELLTYPE_MAP`, refines coarse T/B/mono labels using `REFINEMENT_MARKERS`, and writes `final_cell_type` to cell metadata. Prints a copy-pasteable `CLUSTER_CELLTYPE_MAP` to `logs/05_annotate.log`.

**Config keys:** `SINGLER_REF`, `CLUSTER_CELLTYPE_MAP`, `CONTAMINATION_TYPES`, `REFINEMENT_MARKERS`, `SUBCLUSTER`, `MARKERS`

---

## Step 06 — `06_visualize.R`

**Inputs:** `annotation/integrated_annotated.rds`

**Outputs** (in `results_*/integrated/`):

| File | Description |
|------|-------------|
| `umap_triptych.pdf` | Cluster / sample / cell-type UMAP side-by-side |
| `umap_split_by_sample.pdf` | Per-sample UMAP panels (2 per page) |
| `integrated_umap_*.pdf` | Individual UMAP variants |
| `canonical_markers_feature_*.pdf` | Per-cell-type feature plots |
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

**Inputs:** `annotation/integrated_annotated.rds`

**Outputs** (in `results_*/differential/`):

| File | Description |
|------|-------------|
| `de_<celltype>.csv` | Differentially expressed genes per cell type between samples |
| `volcano_<celltype>.pdf` | Volcano plot per cell type |
| `module_scores_inflammation.pdf` | Inflammatory gene module scores across samples |

**What it does:** Runs `FindMarkers()` for each cell type between samples defined by `SCRNA_CONDITION`. Skips automatically for single-sample runs.

**Config keys:** `SCRNA_CONDITION` (env var)

---

## Step 07 — `07_finalize_reports.R`

**Inputs:** All per-step PDF reports in `qc/`, `doublets/`, `individual/`, `annotation/`, `integrated/`

**Outputs** (in `results_*/reports/`):

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

**Inputs:** `annotation/integrated_annotated.rds`, output files from steps 01–06b

**Outputs** (in `results_*/reports/`):

| File | Description |
|------|-------------|
| `Comparison_report.pdf` | Standalone cross-sample comparison PDF |

**What it does:** Generates a focused comparison document covering sample quality, doublet rates, cell-type composition, integration quality, and DE results across samples. Independent of step 07.

---

## Step 09 — `09_bootstrap_proportions.R`

**Inputs:** `annotation/integrated_annotated.rds`

**Outputs** (in `results_*/proportions/`):

| File | Description |
|------|-------------|
| `bootstrap_proportions.pdf` | Bootstrap-normalised proportions with 95% CI error bars |
| `bootstrap_proportions.csv` | Per-sample per-cell-type mean, lower CI, upper CI |
| `pairwise_chisq.csv` | Pairwise chi-squared tests between samples |

**What it does:** Bootstraps cell-type proportions (1,000 resamples) to produce multinomial 95% confidence intervals. Runs pairwise chi-squared tests to identify statistically significant composition differences.

---

## Step 10 — `10_rarefaction.R`

**Inputs:** `annotation/integrated_annotated.rds`

**Outputs** (in `results_*/rarefaction/`):

| File | Description |
|------|-------------|
| `rarefaction_curves.pdf` | Proportion stability vs cell count per sample |
| `min_capture_depth.csv` | Minimum cells needed for stable proportions per cell type |

**What it does:** Downsamples each sample to increasing cell counts and measures proportion estimate stability. Identifies the minimum cell count required for each cell type to converge to a stable proportion estimate.

---

## Intermediate `.rds` objects

| Object | Location | Written by | Read by |
|--------|----------|-----------|---------|
| `<sample>_seurat_raw.rds` | `qc/<sample>/` | Step 01 | — |
| `<sample>_seurat.rds` | `qc/<sample>/` → `doublets/<sample>/` | Steps 01, 02 | Steps 02, 03 |
| `<sample>_seurat.rds` (cached) | `sample_cache/<sample>/` | Step 03 | Step 04 |
| `integrated_seurat.rds` | `integrated/` | Step 04 | Step 05 |
| `integrated_annotated.rds` | `annotation/` | Step 05 | Steps 06, 06b, 07, 08, 09, 10 |

---

## Related

- [Configuration Reference](reference-config.md) — all config.R parameters
- [Output Files Reference](reference-outputs.md) — full output file inventory
- [How to Re-run from a Specific Step](howto-rerun-steps.md)
- [How to Override Cell Type Annotations](howto-override-annotations.md)
