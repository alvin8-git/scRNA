# Output Files Reference

All outputs land under `Results/results_<sample-names>_<matrix-tag>/`. For example, a run on samples `ES03` and `ES12` with filtered matrices produces `Results/results_ES03-ES12_filtered/`.

---

## Directory layout

```
results_<samples>_filtered/
├── logs/                        # per-step log files (stdout + stderr)
│   ├── 01_load_qc.log
│   ├── 02_doublets.log
│   ├── 03_individual.log
│   ├── 04_integrate.log
│   ├── 05_annotate.log          # includes suggested CLUSTER_CELLTYPE_MAP
│   ├── 06_visualize.log
│   └── 07_finalize_reports.log
├── qc/
│   └── <sample>/
│       ├── <sample>_seurat_raw.rds
│       ├── <sample>_seurat.rds
│       ├── qc_violin.pdf
│       ├── qc_scatter.pdf
│       └── cell_fate.csv
├── doublets/
│   └── <sample>/
│       ├── <sample>_seurat.rds
│       ├── doublet_umap.pdf
│       └── doublet_score_hist.pdf
├── individual/
│   └── <sample>/
│       ├── <sample>_seurat.rds
│       ├── <sample>_umap_clusters.pdf
│       ├── <sample>_umap_celltype.pdf
│       ├── <sample>_top_markers.pdf
│       ├── <sample>_elbow.pdf
│       ├── <sample>_variable_genes.pdf
│       └── <sample>_markers.csv
├── integrated/
│   ├── integrated_umap_before.pdf
│   ├── integrated_umap_after.pdf
│   ├── cluster_resolution_comparison.pdf
│   ├── umap_triptych.pdf
│   ├── umap_split_by_sample.pdf
│   ├── integrated_umap_cluster.pdf
│   ├── integrated_umap_sample.pdf
│   ├── integrated_umap_celltype.pdf
│   ├── canonical_markers_feature_<celltype>.pdf
│   ├── integrated_dotplot.pdf
│   ├── integrated_heatmap.pdf
│   ├── integrated_cluster_markers.csv
│   ├── celltype_proportions_bar.pdf
│   ├── celltype_composition_combined.pdf
│   ├── violin_key_markers.pdf
│   └── visualization_report.pdf
├── annotation/
│   ├── integrated_annotated.rds  ← primary analysis object
│   ├── singler_scores_heatmap.pdf
│   ├── singler_delta_umap.pdf
│   ├── annotation_umap.pdf
│   ├── canonical_markers_dotplot.pdf
│   ├── tcell_subclusters_umap.pdf
│   ├── tcell_subclusters_dotplot.pdf
│   └── tcell_subcluster_summary.csv
├── differential/                 # step 06b, multi-sample only
│   ├── de_<celltype>.csv
│   └── volcano_<celltype>.pdf
├── proportions/                  # step 09
│   ├── bootstrap_proportions.pdf
│   ├── bootstrap_proportions.csv
│   └── pairwise_chisq.csv
├── rarefaction/                  # step 10
│   ├── rarefaction_curves.pdf
│   └── min_capture_depth.csv
└── reports/
    ├── 01-QC_report.pdf
    ├── 02-Doublet_report.pdf
    ├── 03-Individual_report.pdf
    ├── 04-Annotation_report.pdf
    ├── 05-Integrated_report.pdf
    ├── Overall_report.pdf        ← A4-normalised summary with captions
    └── Comparison_report.pdf     ← step 08 standalone comparison
```

---

## Key files in detail

### `annotation/integrated_annotated.rds`

The primary analysis object. A Seurat object with all cells from all samples. Key metadata columns:

| Column | Type | Description |
|--------|------|-------------|
| `sample` | character | Sample name (matches `SAMPLE_NAMES`) |
| `seurat_clusters` | factor | Cluster assignment at `default_res` |
| `singler_label` | character | Raw SingleR label after `SINGLER_NORM` mapping |
| `final_cell_type` | character | Final annotation (SingleR + manual overrides) |
| `scDblFinder.score` | numeric | Doublet probability [0, 1] |
| `scDblFinder.class` | character | `"singlet"` or `"doublet"` |

Load in R:
```r
library(Seurat)
seu <- readRDS("results_ES03-ES12_filtered/annotation/integrated_annotated.rds")
table(seu$final_cell_type)
table(seu$sample, seu$final_cell_type)
```

### `logs/05_annotate.log`

After annotation, scan this file for the suggested `CLUSTER_CELLTYPE_MAP`:

```
CLUSTER_CELLTYPE_MAP <- c(
  "0"  = "CD4 T (naive)",
  "1"  = "NK",
  ...
)
```

Copy-paste into `pipeline/config.R`, correct any mislabelled clusters, then re-run steps 05–07.

### `reports/Overall_report.pdf`

A4-normalised summary report. Each page shows:
- **Bold title banner** (top) — figure name
- **Figure** (centre)
- **Interpretation caption** (bottom) — `Good: … | Bad: …` guide

Requires `pdftools` R package. If unavailable, a simple PDF merge is produced instead (no normalisation or captions).

### `qc/<sample>/cell_fate.csv`

One row per barcode. Columns: `barcode`, `nFeature_RNA`, `nCount_RNA`, `percent.mt`, `fate` (`"kept"` or `"removed"`), `reason` (`"low_genes"`, `"high_genes"`, `"low_counts"`, `"high_counts"`, `"high_mt"`, or `"doublet"`).

Useful for auditing how many cells were removed and why.

### `individual/<sample>/<sample>_markers.csv`

All cluster markers from Wilcoxon rank-sum test. Columns: `gene`, `p_val`, `avg_log2FC`, `pct.1`, `pct.2`, `p_val_adj`, `cluster`. Filtered to `p_val_adj < 0.05` and `avg_log2FC > 0.25`.

### `differential/de_<celltype>.csv`

DE results between conditions (requires `SCRNA_CONDITION`). Columns: `gene`, `avg_log2FC`, `p_val_adj`, `pct.1`, `pct.2`, `cell_type`, `comparison`.

### `proportions/bootstrap_proportions.csv`

Columns: `sample`, `cell_type`, `mean_prop`, `lower_95ci`, `upper_95ci`. Proportions are from 1,000 bootstrap resamples.

---

## Sample cache

`sample_cache/<sample>/<sample>_seurat.rds` — written by step 03, shared across all integration runs involving that sample. If you run `ES03 + ES12` and later run `ES03 + ES14`, ES03 is loaded from cache rather than re-processed.

Delete a subdirectory to force reprocessing:
```bash
rm -rf sample_cache/ES03/
```

---

## Related

- [Pipeline Steps Reference](reference-pipeline-steps.md) — what each step produces
- [Configuration Reference](reference-config.md) — controlling output paths and parameters
