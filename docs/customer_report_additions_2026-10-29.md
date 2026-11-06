# Customer report: candidate additions (2026-10-29)

Scope: what else the `--customer` HTML report (`pipeline/report_template.Rmd`, built by
`pipeline/08b_html_report.R`) could show, ranked by value to a customer with no
bioinformatics staff. Every item is rendered from data (plotly / DT / CSS), never an embedded
screenshot of a PDF page. Sizes are estimates for an 8-sample cm run.

## Size budget

Measured on the 2026-10-28 builds:

| Report | Total | Widget JSON | Embedded images | Rest (plotly.js, DT, CSS) |
|---|---|---|---|---|
| H1 8-sample Leiden | 11.4 MB | 5.9 MB | 0.9 MB | ~4.6 MB |
| ASF5 8-sample | 11.1 MB | 5.7 MB | 0.8 MB | ~4.6 MB |

That leaves roughly **8.5 MB** of headroom under 20 MB. The two largest existing widgets are
the UMAP grid (2.1 MB) and the QC violins (0.8 MB). The JS libraries are paid once, so new
plotly/DT sections cost only their data.

## What each run already has on disk

| Data | Path | Used by the customer report today? |
|---|---|---|
| Pseudobulk edgeR QL DE, every cell type x every timepoint contrast | `differential/DE_all_celltypes.csv` (~870k rows) | No. 08b loads it into `de` but no template chunk renders it |
| Proportion tests per cell type x contrast | `composition_test.csv` | No |
| Bootstrap CIs on proportions | `bootstrap_summary.csv` | No |
| Cluster markers | `export/cluster_markers.csv` | Dot plot only |
| Pseudobulk CPM | `export/pseudobulk_cpm.csv.gz` | Correlation heatmap only |
| Excel bundle of the above | `export/<run>_analysis.xlsx` (45 KB) | No |
| Per-cell table (barcode, sample, cluster, labels, UMAP) | embedded in the HTML as `<run>_per_cell.csv.gz`; same columns as `export/per_cell_assignments.csv.gz` minus QC extras, UMAP rounded to 3 decimals | Yes (Downloads section) |
| Trajectory, pathways, CellChat | `trajectory/`, `pathways/`, `cellchat/` | Empty in all four runs |

## Replicate check (decides what may carry a p-value)

| Run | Inferential contrasts | Descriptive only |
|---|---|---|
| H1 8-sample, ASF5, CLMC | all 6 (2 vs 2 libraries) | none |
| H1 7-sample SoupX (H1D20_1 dropped) | D11vD0, D30vD0, D30vD11 | every contrast involving D20 (n = 1) |

The DE file already carries an `inferential` column; any new section must grey out or hide
p-values where it is FALSE, the same rule the composition panels follow.

## Ranked candidates

| # | Addition | Customer value | Data present? | Validity | Est. size | How |
|---|---|---|---|---|---|---|
| 1 | **Timepoint DE volcano + gene table** per cell type, dropdown for contrast | High: "which genes change between days, in which cell type" | Yes (`DE_all_celltypes.csv`) | Valid: sample-level pseudobulk edgeR QL, 2 vs 2. Descriptive contrasts shown without p-values | 1.0–1.5 MB | One plotly scattergl with a dropdown; per panel all genes at FDR < 0.05 & abs(log2FC) ≥ 1 capped at the top 300, plus 300 sampled background genes; DT table of the top 50 per panel |
| 2 | **Differentiation efficiency panel**: % TNNT2+ cells per sample, CM subtype split, MYH7/MYH6 and MYL2/MYL7 ratios | High: the number a cardio lab reports (flow-cytometry cTnT equivalent), label-independent | Computable from the object at build time | Descriptive per sample; replicate spread shown | < 50 KB | DT table + small plotly bar |
| 3 | **Composition statistics**: proportion change with FDR and bootstrap 95% CI | High: turns the stacked bars into "significantly up/down" | Yes (`composition_test.csv`, `bootstrap_summary.csv`) | t-test on 2 vs 2 proportions: low power, say so; descriptive rows marked | < 30 KB | DT with coloured FDR; forest-style plotly of CIs |
| 4 | **Pathway enrichment (GSEA)** per cell type x contrast on the DE statistic | Medium-high: biology in words, not gene lists | No: needs a run of fgsea/clusterProfiler `gseGO` (both installed, plus org.Hs.eg.db). Hallmark/MSigDB not installed (`msigdbr` absent) | Valid if ranked by the edgeR signed statistic over all tested genes; GO terms collapsed for redundancy | 150–300 KB | DT of top terms + plotly dot plot |
| 5 | **Pseudotime (cardiac lineage)** | Medium-high for a differentiation series | No: `pipeline/projects/cm/trajectory.R` (monocle3 1.4.27 installed) must run first, ~10–20 min per run | A hypothesis, not lineage tracing. Validate against the real collection day (pseudotime should rise D0 → D30); report that correlation | 300–500 KB | Rasterised UMAP coloured by pseudotime (same scattermore route as the QC panels), plotly gene-trend lines for ~18 markers, pseudotime density by day |
| 6 | **Gene explorer**: type a gene, see mean expression per cell type x sample | Medium: answers the customer's own gene questions | Yes (pseudobulk CPM) | Descriptive | 1–2 MB for ~3,000 genes (all 20k genes ≈ 8 MB, too big) | Compact JSON + plotly heatmap redrawn on search |
| 7 | **Data downloads** embedded (Excel bundle) | Medium: no server needed to hand over tables | Yes (45 KB xlsx) | n/a | ~60 KB | `data:` link; per-cell assignments (2.8 MB) only if asked |
| 8 | **Cell-cycle phase** per cell type and day | Low-medium: proliferation should fall with maturation | No: `CellCycleScoring` at build time (seconds) | Descriptive | < 30 KB | Stacked plotly bar |

Items 1–3 + 7 together ≈ 1.6 MB; all eight ≈ 4–5 MB, still under 20 MB per report.

## Not recommended

- **Volcano plot of highly variable genes.** HVGs have no fold change and no p-value; they
  are chosen for variance across all cells, with no comparison behind them. A volcano needs a
  contrast. The valid plot is item 1 (timepoint DE). An HVG mean-variance plot is possible but
  of no use to a customer.
- **Cell-level DE p-values** (Wilcoxon per cell): pseudoreplication; thousands of cells from
  one library are not independent replicates. Item 1 already uses sample-level pseudobulk.
- **CellChat**: not installed (needs a new dependency) and `cellchat/` is empty.
- **RNA velocity**: needs spliced/unspliced counts; the C4 `filter_matrix` has only total counts.

## Recommendation

Build 1, 2, 3 and 7 first: all four use data already on disk or seconds of computation, and
together add about 1.6 MB. Then 4 and 5, which need a pipeline step per run before the
report can show them.

## Analysis object supplement (08f, 2026-11-06)

The HTML report stays the main deliverable; the object is too large to embed or email, so `08f` writes
it next to the tables as a supplement, in `export/`:

| File | Contents |
|---|---|
| `<run>_seurat_slim.rds` | Seurat object, `xz`: RNA raw counts, UMAP, metadata `orig.ident sample nCount_RNA nFeature_RNA percent.mt seurat_clusters cell_type` (+ `cell_type_ref` / `cell_type_denovo` when present). ~31 MB at 11.8k cells, ~131 MB at 47.9k (xz took 93-144 s at 11.8k cells) |
| `<run>.h5ad` | AnnData: `X` = raw counts (CSR), `obs` = same columns, `var` = gene names, `obsm/X_umap`. Written by hand with `hdf5r` (`pipeline/object_export.R`); no anndata/SeuratDisk needed. Uncompressed (gzip chunks made it larger): 62 MB at 11.8k cells, written in ~5 s |

`xz` is ~8x slower than gzip, so large runs spend minutes here. Skip with `SCRNA_NO_OBJECT_EXPORT=1` or
`Rscript pipeline/08f_customer_export.R <run> --no-object`. `SCRNA_EXPORT_DIR` redirects the whole `export/`
folder (used for read-only test runs). A failure in the object step is logged and does not stop the CSV/xlsx export.
