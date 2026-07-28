# Output Files Reference

All outputs land under `Results/results_<sample-names>_<matrix-tag>/`. For example, a run on samples `ES03` and `ES12` with filtered matrices produces `Results/results_ES03-ES12_filtered/`.

---

## Directory layout

The final PDFs sit at the top level of the run directory, not in a `reports/` subfolder — `DIRS$reports`
is the run directory itself. Only step 08b's interactive HTML gets its own `reports/` folder.

```
results_<samples>_filtered/
├── 01-QC_report.pdf              ← step 07 merges the per-step PDFs into these five
├── 02-Doublet_report.pdf
├── 03-Individual_report.pdf
├── 04-Annotation_report.pdf
├── 05-Integrated_report.pdf
├── Overall_report.pdf            ← A4-normalised summary with captions
├── Comparison_report.pdf         ← step 08, run on request
├── bootstrap_proportions_report.pdf  ← step 09, run on request
├── bootstrap_summary.csv
├── rarefaction_report.pdf        ← step 10, run on request
├── rarefaction_summary.csv
├── logs/                         # per-step logs, named <step>_<script>.log
│   ├── 01_01_load_qc.log
│   ├── 02_02_doublets.log
│   ├── 03_03_individual.log
│   ├── 04_04_integrate.log
│   ├── 05_05_annotate.log        # includes suggested CLUSTER_CELLTYPE_MAP
│   ├── 06_06_visualize.log
│   ├── 06b_06b_differential.log
│   └── 07_07_finalize_reports.log
├── qc/
│   ├── <sample>_violin_qc.pdf
│   ├── <sample>_scatter_qc.pdf
│   ├── cell_fate.csv
│   ├── qc_summary_table.csv
│   └── qc_report.pdf
├── doublets/
│   ├── <sample>_doublet_umap.pdf
│   ├── <sample>_doublet_score_hist.pdf
│   └── doublets_report.pdf
├── individual/
│   ├── individual_report.pdf
│   └── <sample>/
│       ├── <sample>_filtered.rds  ← step 01 checkpoint
│       ├── <sample>_singlets.rds  ← step 02 checkpoint
│       ├── <sample>_seurat.rds    ← step 03 checkpoint
│       ├── <sample>_elbow.pdf
│       ├── <sample>_hvg.pdf
│       ├── <sample>_pc_heatmaps.pdf
│       ├── <sample>_umap_cluster.pdf
│       ├── <sample>_umap_sample.pdf
│       ├── <sample>_umap_markers.pdf
│       ├── <sample>_dotplot_markers.pdf
│       └── <sample>_cluster_markers.csv
├── integrated/
│   ├── integrated_seurat.rds      ← step 04 output
│   ├── integrated_annotated.rds   ← step 05 output; primary analysis object
│   ├── harmony_before_after.pdf
│   ├── cluster_resolution_comparison.pdf
│   ├── integrated_umap_cluster.pdf
│   ├── integrated_umap_sample.pdf
│   ├── integrated_umap_celltype.pdf
│   ├── integrated_umap_split_sample.pdf
│   ├── umap_triptych.pdf
│   ├── umap_split_by_sample.pdf
│   ├── feature_<marker-group>.pdf
│   ├── integrated_dotplot.pdf
│   ├── integrated_heatmap.pdf
│   ├── integrated_cluster_markers.csv   # only when MARKERS$compute_integrated = TRUE
│   ├── celltype_proportions_bar.pdf
│   ├── celltype_counts_bar.pdf
│   ├── celltype_composition_combined.pdf
│   ├── integration_report.pdf
│   └── visualization_report.pdf
├── annotation/
│   ├── celltype_umap.pdf
│   ├── canonical_markers_dotplot.pdf   # read this to fill CLUSTER_CELLTYPE_MAP
│   ├── singler_labels_umap.pdf
│   ├── singler_scores_heatmap.pdf
│   ├── singler_delta_umap.pdf
│   ├── sctype_labels_umap.pdf
│   ├── singler_vs_sctype_comparison.csv
│   ├── consensus_annotation.csv
│   ├── contamination_summary.pdf
│   ├── cluster_annotation_table.csv    # per-cluster majority label
│   ├── tcell_subclusters_umap.pdf
│   ├── tcell_subclusters_dotplot.pdf
│   ├── tcell_subcluster_summary.csv
│   ├── reference_transfer_cells.csv.gz      # step 05r, only with REFERENCE_MODEL set
│   ├── reference_transfer_composition.csv
│   └── annotation_report.pdf
├── differential/                 # step 06b, multi-sample only
│   ├── DE_<celltype>.csv
│   ├── DE_all_celltypes.csv
│   ├── DE_summary.csv
│   ├── volcano_<celltype>.pdf
│   ├── module_score_<module>.pdf
│   └── differential_report.pdf
├── benchmark/                    # step 08c, only with REFERENCE_MODEL set
│   ├── concordance.csv
│   ├── wholeblood_signature.csv
│   └── benchmark_report.md
└── reports/                      # step 08b
    ├── <run>_report.html         ← self-contained interactive report
    └── build_report.log
```

---

## Key files in detail

### `integrated/integrated_annotated.rds`

The primary analysis object. A Seurat object with all cells from all samples. Key metadata columns:

| Column | Type | Description |
|--------|------|-------------|
| `sample` | character | Sample name (matches `SAMPLE_NAMES`) |
| `seurat_clusters` | factor | Cluster assignment at `default_res` |
| `singler_label` | character | Raw SingleR label, before normalisation |
| `singler_label_clean` | character | SingleR label after pruning + `SINGLER_NORM` mapping |
| `singler_pruned` | character | `NA` where SingleR pruned the call as low-confidence |
| `singler_delta` | numeric | Score gap between the best and next-best label |
| `sctype_label` | character | Independent ScType marker-based call, for cross-checking |
| `cell_type` | character | **The label every downstream step uses.** SingleR + contamination overrides + `CLUSTER_CELLTYPE_MAP` + sub-type refinement |
| `scDblFinder.score` | numeric | Doublet probability [0, 1] |
| `scDblFinder.class` | character | `"singlet"` or `"doublet"` |

Two more columns appear only after step 05r has run (`REFERENCE_MODEL` set): `cell_type_ref`, the
run-independent frozen-reference label, and `cell_type_denovo`, preserving the original de-novo call.
When that output exists, `config.R`'s `apply_reference_labels()` promotes `cell_type_ref` into
`cell_type`, so plots and tables switch to the frozen labels without any step needing to change.

There is no `final_cell_type` column on the Seurat object — that name belongs to the per-cluster
majority-label column of `annotation/cluster_annotation_table.csv`.

Load in R:
```r
library(Seurat)
seu <- readRDS("results_ES03-ES12_filtered/integrated/integrated_annotated.rds")
table(seu$cell_type)
table(seu$sample, seu$cell_type)
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

### `Overall_report.pdf`

A4-normalised summary report. Each page shows:
- **Bold title banner** (top) — figure name
- **Figure** (centre)
- **Interpretation caption** (bottom) — `Good: … | Bad: …` guide

Requires `pdftools` R package. If unavailable, a simple PDF merge is produced instead (no normalisation or captions).

### `qc/cell_fate.csv`

The cell-attrition funnel, **one row per sample** (not per barcode). Columns: `Sample`,
`GEM_barcodes`, `After_load`, `Removed_load_low_genes`, `After_QC`, `Removed_QC`,
`After_doublet_removal`, `Removed_doublets`, `Total_removed`, `Pct_retained`.

Useful for auditing how many cells were removed at each stage. `qc/qc_summary_table.csv` sits beside it
with the per-sample medians: `Sample`, `Cells`, `Genes`, `Median_nFeature`, `Median_nCount`,
`Median_pct_mt`.

### `individual/<sample>/<sample>_cluster_markers.csv`

All cluster markers from the Wilcoxon rank-sum test. Columns: `p_val`, `avg_log2FC`, `pct.1`, `pct.2`, `p_val_adj`, `cluster`, `gene`.

### `differential/DE_<celltype>.csv`

DE results between conditions (requires `SCRNA_CONDITION`). Columns: `gene`, `avg_log2FC`, `p_val_adj`, `pct.1`, `pct.2`, `cell_type`, `comparison`.

### `bootstrap_summary.csv`

Written to the run directory root by step 09. Columns: `sample`, `cell_type`, `total`, `n`, `prop`,
`ci_lo`, `ci_hi` (observed proportion with its analytical multinomial 95% CI), then `boot_mean`,
`boot_lo`, `boot_hi` from 1,000 bootstrap resamples down to the smallest sample.

### `annotation/cluster_annotation_table.csv`

One row per cluster. Columns: `seurat_clusters`, `n_cells`, `final_cell_type` (the majority
`cell_type` in that cluster), `singler_majority`. This is the only place the name
`final_cell_type` is used.

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
