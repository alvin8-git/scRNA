# Pipeline Documentation

## Overview

This pipeline analyses human PBMC single-cell RNA-seq data produced by the DNB C4 platform (dnbc4tools v3.0). It performs quality control, doublet detection, individual sample analysis, Harmony batch integration, automated and manual cell type annotation, and publication-quality visualisation using Seurat v5 in R.

---

## Data Format Notes (DNB C4 vs 10X)

The DNB C4 output is compatible with Seurat's `Read10X()` but has one key difference:

| | 10X Genomics | DNB C4 |
|--|--|--|
| `features.tsv.gz` columns | 2 (ENSEMBL ID, gene symbol) | 3 (ENSEMBL ID, gene symbol, feature type) |
| Barcode format | `ACGTACGT-1` | `CELL1_N3` |
| Matrix format | MTX (gzip) | MTX (gzip) — identical |

**Critical:** always pass `gene.column = 2` to `Read10X()` to use gene symbols instead of ENSEMBL IDs with version numbers (e.g. `ENSG00000243485.5`).

---

## Pipeline Steps

### Step 01 — Load & QC (`01_load_qc.R`)

**Input:** `H1/filter_matrix/`, `H2/filter_matrix/`  
**Output:** `results/qc/*.pdf`, `results/individual/H{1,2}/H{1,2}_filtered.rds`

Loads each sample with `Read10X(gene.column = 2)`, adds sample metadata, computes:
- `nFeature_RNA` — number of detected genes per cell
- `nCount_RNA` — total UMI count per cell
- `percent.mt` — percentage of reads mapping to mitochondrial genes (`^MT-`)
- `percent.ribo` — percentage mapping to ribosomal genes (`^RP[SL]`) (diagnostic only)

Default QC thresholds (adjust in `config.R`):

| Metric | Min | Max | Rationale |
|--------|-----|-----|-----------|
| nFeature_RNA | 200 | 5000 | Removes empty droplets and potential doublets |
| nCount_RNA | 500 | 25000 | Removes debris and outlier high-UMI cells |
| percent.mt | — | 20% | PBMC can legitimately reach 15-20% |

**Inspect** `H1_violin_qc.pdf` and `H2_violin_qc.pdf` before accepting defaults. For H1 (only 386 input cells), relax `max_features` to 6000 if a genuine population falls above 5000.

Barcodes are prefixed with the sample name (`RenameCells(add.cell.id = sample_name)`) so `H1_CELL1_N3` and `H2_CELL1_N3` remain distinct after merging.

---

### Step 02 — Doublet Detection (`02_doublets.R`)

**Input:** `H{1,2}_filtered.rds`  
**Output:** `results/doublets/*.pdf`, `H{1,2}_singlets.rds`

Uses **scDblFinder** (Bioconductor), which is natively compatible with Seurat v5. The tool:
1. Converts the Seurat object to `SingleCellExperiment`
2. Generates artificial doublets by random cell fusion
3. Trains a gradient boosting classifier to score each cell
4. Labels cells as `singlet` or `doublet`

Expected doublet rates used:
- H1 (386 cells): **3.1%** → ~12 expected doublets
- H2 (967 cells): **7.7%** → ~75 expected doublets

Adjust rates in `config.R → DOUBLET$doublet_rate`.

**Outputs per sample:**
- `*_doublet_umap.pdf` — UMAP with doublets highlighted in red
- `*_doublet_score_hist.pdf` — score distribution (bimodal = clean separation)

---

### Step 03 — Individual Analysis (`03_individual.R`)

**Input:** `H{1,2}_singlets.rds`  
**Output:** `results/individual/H{1,2}/*.pdf`, `H{1,2}_seurat.rds`

Full per-sample Seurat workflow:

1. **Normalisation** — `LogNormalize` with scale factor 10,000
2. **HVG selection** — `FindVariableFeatures(nfeatures = 2000, method = "vst")`
3. **Scaling** — `ScaleData(vars.to.regress = "percent.mt")` removes mitochondrial confounding
4. **PCA** — 30 PCs computed; elbow plot + PC heatmaps saved
5. **UMAP** — `dims = 1:20`, `seed = 42`
6. **Clustering** — `FindNeighbors` + `FindClusters` at resolutions 0.3, 0.4, 0.5, 0.6, 0.8; default `res = 0.5`
7. **Marker genes** — `FindAllMarkers(only.pos = TRUE, test = "wilcox", min.pct = 0.1, logfc.threshold = 0.25)`

For H1 (~350 cells), consider using `resolution = 0.4` to avoid over-clustering.

---

### Step 04 — Harmony Integration (`04_integrate.R`)

**Input:** `H{1,2}_seurat.rds`  
**Output:** `results/integrated/*.pdf`, `integrated_seurat.rds`

1. **Merge** — `merge(seu_H1, seu_H2)` then `JoinLayers()` (required in Seurat v5 after merge)
2. **Re-normalise** on merged object (normalisation must be done after merge)
3. **PCA** on merged HVGs
4. **Harmony** — runs on PCA embedding, returns a `harmony` reduction
   - `theta = 2`: diversity penalty; increase to 3-4 if samples don't intermix
   - `lambda = 1`: ridge regression penalty
   - `nclust = 50`: k-means clusters (~nCells/25 rule of thumb)
5. **UMAP** on Harmony embedding
6. **Clustering** on Harmony-corrected neighbours

**Key diagnostic:** `harmony_before_after.pdf` — H1/H2 should intermix after correction while cell type clusters remain separated. If cell types blur together, reduce `theta`.

---

### Step 05 — Annotation (`05_annotate.R`)

**Input:** `integrated_seurat.rds`  
**Output:** `results/annotation/*.pdf`, `cluster_annotation_table.csv`, `integrated_annotated.rds`

**Part 1: SingleR automated annotation**
- Reference: `HumanPrimaryCellAtlasData` (broad PBMC cell types)
- Run at single-cell level with `fine.tune = TRUE`
- Low-confidence calls pruned (`prune = TRUE`)
- Delta score (confidence margin) saved as `singler_delta`

**Part 2: Canonical marker dot plot**
- Dot plot of all canonical PBMC markers across clusters
- Use this to manually assign `CLUSTER_CELLTYPE_MAP` in `config.R`

**Part 3: Label assignment**
- If `CLUSTER_CELLTYPE_MAP` is set: uses manual labels
- Otherwise: uses SingleR majority-vote per cluster as fallback

**Canonical marker genes used:**

| Cell Type | Markers |
|-----------|---------|
| T cell (pan) | CD3D, CD3E |
| CD4+ T | CD4, IL7R, CCR7 |
| CD8+ T | CD8A, CD8B, GZMK |
| NK | NKG7, GNLY, KLRD1 |
| B cell | MS4A1, CD79A, CD19 |
| CD14+ Mono | CD14, LYZ, CST3, S100A8 |
| FCGR3A+ Mono | FCGR3A, MS4A7 |
| DC | FCER1A, CLEC9A |
| Platelet | PPBP, PF4 |

---

### Step 06 — Visualisation (`06_visualize.R`)

**Input:** `integrated_annotated.rds`  
**Output:** `results/integrated/*.pdf`

All figures saved at 300 DPI as PDF.

| File | Description |
|------|-------------|
| `umap_triptych.pdf` | 3-panel: cluster / sample / cell type |
| `integrated_umap_cluster.pdf` | UMAP coloured by Seurat cluster |
| `integrated_umap_sample.pdf` | UMAP coloured by H1 vs H2 |
| `integrated_umap_celltype.pdf` | UMAP coloured by cell type |
| `umap_split_by_sample.pdf` | Side-by-side H1 vs H2 UMAP |
| `harmony_before_after.pdf` | Batch effect before/after Harmony |
| `feature_T_cells.pdf` | Feature plots: T cell markers |
| `feature_NK.pdf` | Feature plots: NK markers |
| `feature_B_cells.pdf` | Feature plots: B cell markers |
| `feature_monocytes.pdf` | Feature plots: monocyte markers |
| `feature_DC_platelet.pdf` | Feature plots: DC + platelet markers |
| `feature_all_canonical_markers.pdf` | All markers in one composite |
| `integrated_dotplot.pdf` | Dot plot: cell type × canonical markers |
| `integrated_heatmap.pdf` | Heatmap: top 3 markers per cluster |
| `celltype_proportions_bar.pdf` | Cell type proportions per sample |
| `celltype_counts_bar.pdf` | Cell type counts per sample |
| `violin_key_markers.pdf` | Violin plots for 8 lineage markers |

---

## Configuration Reference (`config.R`)

All parameters are centralised in `config.R`. Every script sources this file at startup.

### Key parameters to adjust

```r
# QC thresholds
QC$min_features   <- 200    # lower = keep more cells
QC$max_features   <- 5000   # raise to 6000 if legitimate cells exceed this
QC$max_percent_mt <- 20     # lower to 15 for stricter filtering

# Clustering resolution
CLUSTER$default_res <- 0.5  # raise for more clusters, lower for fewer

# Harmony diversity penalty
HARMONY$theta <- 2           # raise to 3-4 if samples don't mix; lower to 1 if over-integrated

# Manual annotation (fill after inspecting canonical_markers_dotplot.pdf)
CLUSTER_CELLTYPE_MAP <- c(
  "0" = "CD4 T",
  "1" = "CD14+ Mono",
  # ...
)
```

---

## Re-running Individual Steps

```bash
conda activate scrna_seurat

# Re-run a single step
bash /data/alvin/scRNA/pipeline/run_pipeline.sh 05

# Re-run annotation and visualisation only (most common after manual labelling)
bash /data/alvin/scRNA/pipeline/run_pipeline.sh 05 06

# Re-run from doublet detection onwards
bash /data/alvin/scRNA/pipeline/run_pipeline.sh 02 03 04 05 06
```

---

## Saved R Objects

| File | Content | Size |
|------|---------|------|
| `results/individual/H1/H1_filtered.rds` | H1 post-QC Seurat object | ~5 MB |
| `results/individual/H1/H1_singlets.rds` | H1 post-doublet-removal | ~5 MB |
| `results/individual/H1/H1_seurat.rds` | H1 fully processed (UMAP, clusters) | ~15 MB |
| `results/individual/H2/H2_*.rds` | Same for H2 | — |
| `results/integrated/integrated_seurat.rds` | Merged + Harmony integrated | ~30 MB |
| `results/integrated/integrated_annotated.rds` | Final annotated object | ~30 MB |

Load any object for interactive exploration:
```r
library(Seurat)
seu <- readRDS("/data/alvin/scRNA/results/integrated/integrated_annotated.rds")
DimPlot(seu, group.by = "cell_type")
```
