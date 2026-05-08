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
- Reference controlled by `SINGLER_REF` in `config.R`:
  - `"HumanPrimaryCellAtlas"` (default) — broad PBMC cell types (~30 types)
  - `"MonacoImmune"` — blood-optimised, resolves CD4/CD8/Treg/γδ T, monocyte subtypes, pDC/mDC (29 types); set automatically by `bat` species override
- Run at single-cell level with `fine.tune = TRUE`
- Low-confidence calls pruned (`prune = TRUE`)
- Delta score (confidence margin) saved as `singler_delta`
- Raw labels normalised via `SINGLER_NORM` lookup — covers HumanPrimaryCell granulocyte/erythroid lineages + full Monaco label.main set

**Part 1b: scType marker-based scoring**
- Runs alongside SingleR with no reference transcriptome — uses `MARKERS` gene lists directly
- Genes weighted by inverse frequency (markers found in fewer cell types score higher, reducing false positives from shared markers e.g. FCGR3B shared between Neutrophil and FCGR3A_mono)
- Filters genes against `rownames(scale.data)` — HVG-safe (scale.data only contains the 2000 HVGs)
- Outputs: `sctype_labels_umap.pdf`, `singler_vs_sctype_comparison.csv`
- `sctype_label` added to metadata for cross-reference with SingleR

**Part 2: Canonical marker dot plot**
- Dot plot of all canonical PBMC markers across clusters
- Use this to manually assign `CLUSTER_CELLTYPE_MAP` in `config.R`

**Part 3: Label assignment**
- If `CLUSTER_CELLTYPE_MAP` is set: uses manual labels
- Otherwise: uses SingleR majority-vote per cluster as fallback

**Part 4 (contamination override)**
- After majority-vote assignment, per-cell SingleR labels for types listed in `CONTAMINATION_TYPES` (Neutrophil, RBC, HSPC, Platelet, Basophil, Eosinophil, Mast cell) override the cluster label
- This ensures scattered contamination cells appear on UMAP even if they are a minority in their cluster
- Override count logged per type: `[CONTAM] Preserved N <type> cell(s)`
- `contamination_summary.pdf` — stacked bar of contamination counts per sample

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

> **Note on split-by-sample UMAP:** panels are produced at ≤2 per page. A shared cell-type legend (extracted from the full merged object via `cowplot::get_legend()`) is placed as a right column; each per-sample panel uses `NoLegend()` to prevent duplicate legends.

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

### Step 09 — Bootstrap Proportion Comparison (`09_bootstrap_proportions.R`)

**Input:** `integrated_annotated.rds`
**Output:** `reports/bootstrap_proportions_report.pdf`, `reports/bootstrap_summary.csv`

Normalises cell type proportion estimates across samples with different capture depths:

1. **Observed proportions** — per sample with multinomial 95% CI (Wilson-score per cell type)
2. **Bootstrap downsampling** — 1000 draws at n_min (smallest sample size) without replacement; removes capture-depth bias
3. **Per-cell-type CI plot** — dot + error bar per sample × cell type (faceted)
4. **Pairwise chi-squared** — all sample pairs tested; tile heatmap with p-values

> Note: bootstrap mean ≈ observed proportion (same point estimate, different uncertainty). The value of bootstrap is the CI at a normalised depth, not a shifted mean.

Run: `bash pipeline/run_pipeline.sh <samples> 09`

---

### Step 10 — Rarefaction Analysis (`10_rarefaction.R`)

**Input:** `integrated_annotated.rds`
**Output:** `reports/rarefaction_report.pdf`, `reports/rarefaction_summary.csv`

Empirically determines minimum capture depth for stable proportion estimates:

1. Uses the **largest sample** as ground truth
2. Subsamples at increasing depths (100, 250, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000)
3. 1000 bootstrap draws per depth; computes 95% CI width and RMSE vs ground truth
4. Fits theoretical CI ~ a/√n per cell type
5. Reports minimum n where CI ≤ 2× asymptotic value

**Findings for bat whole blood (ES03_newkit ground truth, n=6,520):**

| Precision target | Min cells needed |
|---|---|
| ±5% (crude) | ~500 |
| ±2% (moderate) | ~2,000 |
| ±1% (good) | ~4,000–5,000 |
| ±0.5% (high) | 6,000+ |

**Recommendation: ≥5,000 cells per bat whole blood sample** for full ±1% precision across all cell types.

Run: `bash pipeline/run_pipeline.sh <samples> 10`

---

### Step 08 — Comparison Report (`08_comparison_report.R`)

**Input:** existing output files in `RESULTS_DIR` (QC CSVs, doublet PDFs, composition CSVs, integrated PDFs, DE files)  
**Output:** `Comparison_report.pdf` in `RESULTS_DIR`

Standalone cross-sample summary report with sections:

| Section | Source files |
|---------|-------------|
| Cover page | pipeline metadata |
| Sample Quality & Cell Fate | `qc/cell_fate.csv`, QC violins/scatter PDFs |
| Doublet Rates | doublet UMAP + histogram PDFs per sample |
| Cell Type Composition | `integrated/celltype_composition_combined.pdf` |
| Integrated UMAP | all pages of `integrated/umap_split_by_sample.pdf` |
| Key Marker Expression | `integrated/integrated_dotplot.pdf`, violin PDFs |

Run it alone: `bash run_pipeline.sh <samples> 08`

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

## Bat (*Eonycteris spelaea*) Whole-Blood Support

Pass `bat` as the first argument to `run_pipeline.sh` to activate species-specific overrides:

```bash
bash pipeline/run_pipeline.sh bat /path/ES03 /path/ES12
```

The `bat` keyword exports `SCRNA_SPECIES=bat`; `config.R` reads it and automatically applies:

| Setting | Human default | Bat override |
|---------|--------------|-------------|
| `SINGLER_REF` | `HumanPrimaryCellAtlas` | `MonacoImmune` |
| `CLUSTER$default_res` | 0.5 | 1.0 |
| `CLUSTER$resolutions` | 0.3–0.8 | 0.3, 0.5, 0.8, 1.0 |
| `MARKERS$CD14_mono` | CD14, LYZ, CST3, S100A8 | CD14, LYZ, S100A8, S100A9, **CSF1R** (CST3 absent; CSF1R added) |
| `MARKERS$FCGR3A_mono` | FCGR3A, MS4A7 | FCGR2A, FCGR3B, MS4A7 (FCGR3A absent) |
| `MARKERS$Neutrophil` | FCGR3B, CSF3R, CXCR2, CEACAM8 | FCGR3B, CSF3R, CXCR2, CEACAM6, **IDO1** (CEACAM8 absent; IDO1 bat-specific) |
| `MARKERS$gamma_delta_T` | — | TRDC, TRGC1, TRGC2 |
| `CONTAMINATION_TYPES` | Neutrophil, RBC, HSPC, Platelet, … | Basophil, Eosinophil, Mast cell only |
| B cell subtype markers | Uses IGHD, IGHM, IGHG1 | Uses TCL1A, IL4R, CD24, FCER2 (isotypes absent) |

**Bat reference GTF:** `Bat/ESpe_merged_fullref.gtf` — 30,180 genes (30,011 original + 169 new TCR loci merged from `ESpe_Peaks2UTRed_genome_tcellgenes_ZF_Dec2024_v2.gtf`). The new GTF has only 207 `gene`-level feature lines (TCR genes only); the remaining ~30k genes are represented in transcript/exon/CDS records — this does not affect GEM quantification since STAR counts via exon features.

**Expected bat whole-blood cell types** (4-sample run Sample6/Sample7/10/ES03_newkit, res=0.5, MonacoImmune):

| Cell type | Typical range | Notes |
|---|---|---|
| Monocyte | 23–80% | Dominant; bat whole blood is constitutively myeloid-high |
| CD4 T | 7–44% | Wide inter-individual variation; Sample10 outlier (44%) |
| CD8 T | 0.4–20% | Confirmed by CD3E 90.7%, NCAM1 0.5% — SingleR calls NK incorrectly |
| Neutrophil | 1–10% | CSF3R 82%, CXCR2 50%, IDO1 44% — genuine; human references under-call |
| DC | 1–3% | FLT3 73%, HLA-DR 98% — SingleR absorbs into Monocyte without override |
| B cell | 0.7–4% | Low in bat whole blood; consistent across samples |
| HSPC | 1–4% | Haematopoietic progenitors; normal in bat whole blood |
| NK (true) | 0–1% | NCAM1/CD56+; minor population |

**Mandatory CLUSTER_CELLTYPE_MAP overrides for bat whole blood** (verify cluster numbers per run):
- NK-labelled clusters with CD3E >80% → **CD8 T** (CTL, not NK)
- Monocyte-labelled clusters with FLT3 >40% and HLA-DR >70% → **DC**
- HSPC-labelled clusters with S100A8/A9 100% and CSF3R >30% → **Neutrophil**

**Minimum recommended capture depth:** ≥5,000 cells per sample for ±1% CI on all populations (rarefaction analysis using ES03_newkit as ground truth).

---

## Re-running Individual Steps

```bash
conda activate scrna_seurat

# Re-run a single step
bash /data/alvin/scRNA/pipeline/run_pipeline.sh 05

# Re-run annotation and visualisation only (most common after manual labelling)
bash /data/alvin/scRNA/pipeline/run_pipeline.sh /path/S1 /path/S2 05 06 07

# Re-run from doublet detection onwards
bash /data/alvin/scRNA/pipeline/run_pipeline.sh /path/S1 /path/S2 02 03 04 05 06 07

# Re-generate comparison report only (all earlier steps already done)
bash /data/alvin/scRNA/pipeline/run_pipeline.sh /path/S1 /path/S2 08

# Re-run visualisation + both reports
bash /data/alvin/scRNA/pipeline/run_pipeline.sh /path/S1 /path/S2 06 07 08
```

> Steps 01–03 skip automatically via per-sample RDS cache (`sample_cache/<name>/`) when the same sample has been processed in a previous run. Delete `sample_cache/<name>/` to force reprocessing.

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
