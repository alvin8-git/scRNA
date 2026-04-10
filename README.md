# scRNA-seq Seurat Pipeline

[![R](https://img.shields.io/badge/R-4.3.3-276DC3?logo=r)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-v5.1.0-blue)](https://satijalab.org/seurat/)
[![Harmony](https://img.shields.io/badge/Harmony-batch_correction-green)](https://github.com/immunogenomics/harmony)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A modular, end-to-end single-cell RNA-seq analysis pipeline built on **Seurat v5**. Supports single-sample and multi-sample (Harmony-integrated) workflows with automated cell type annotation, publication-quality figures, and combined PDF reports.

Developed for human PBMC data from the **DNB C4** sequencing platform but adaptable to any species or tissue type.

---

## Features

- **Automatic sample detection** — pass one or more sample folder paths; multi-sample runs trigger Harmony batch correction automatically
- **7-step modular pipeline** — each step saves intermediate `.rds` objects so you can re-run from any point
- **Parallel execution** — `future` multicore for Seurat steps, `BiocParallel` for scDblFinder and SingleR
- **Automated cell type annotation** — SingleR (HumanPrimaryCellAtlasData) with optional manual override via `CLUSTER_CELLTYPE_MAP`
- **Combined PDF reports** — one report per category (QC, doublets, individual, annotation, integrated) with titled figure pages
- **Dynamic results directory** — named after input samples (e.g. `results_H1andH2/`)

---

## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 01 | `01_load_qc.R` | Load 10x matrices, compute QC metrics, filter low-quality cells |
| 02 | `02_doublets.R` | Doublet detection with scDblFinder |
| 03 | `03_individual.R` | Per-sample: normalize, HVG, PCA, UMAP, cluster, markers |
| 04 | `04_integrate.R` | Merge samples + Harmony batch correction (or single-sample UMAP) |
| 05 | `05_annotate.R` | SingleR automated annotation + canonical marker plots |
| 06 | `06_visualize.R` | Publication-quality UMAP, dot plots, heatmaps, violin plots |
| 07 | `07_finalize_reports.R` | Merge all PDFs into 5 final combined reports |

---

## Requirements

- **Linux** (tested on Ubuntu 20.04)
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://github.com/mamba-org/mamba)
- ~8 GB RAM minimum; 16–32 GB recommended for multi-sample runs
- ~5 GB disk for conda environment

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/alvin8-git/scRNA.git
cd scRNA
```

### 2. Create the conda environment

```bash
# Uses mamba for faster solving (falls back to conda if mamba not found)
bash pipeline/setup_env.sh
```

This installs R 4.3.3 with all required packages:
`Seurat`, `harmony`, `scDblFinder`, `SingleR`, `celldex`, `BiocParallel`, `ggplot2`, `patchwork`, `cowplot`, `pheatmap`, `viridis`, `dplyr`

Installation takes **10–20 minutes** on first run.

### 3. Activate the environment

```bash
conda activate scrna_seurat
```

---

## Usage

### Input data format

Each sample folder should contain a 10x Genomics-compatible matrix in one of these subfolder layouts:

```
SampleA/
└── filter_matrix/          # preferred: pre-filtered barcodes
    ├── barcodes.tsv.gz
    ├── features.tsv.gz
    └── matrix.mtx.gz
```

Or `filtered_feature_bc_matrix/` — the pipeline detects both automatically.

### Run the pipeline

```bash
# Single sample
bash pipeline/run_pipeline.sh /path/to/SampleA

# Two samples (Harmony integration enabled automatically)
bash pipeline/run_pipeline.sh /path/to/SampleA /path/to/SampleB

# Any number of samples
bash pipeline/run_pipeline.sh /path/to/S1 /path/to/S2 /path/to/S3

# Specific steps only (e.g. re-run annotation and visualization)
bash pipeline/run_pipeline.sh /path/to/SampleA /path/to/SampleB 05 06 07
```

### Test run (with example data)

If you have the H1 and H2 example data placed under `H1/` and `H2/`:

```bash
conda activate scrna_seurat
bash pipeline/run_pipeline.sh H1 H2
```

Expected output in `results_H1andH2/`:

```
results_H1andH2/
├── logs/                        # per-step log files
├── qc/                          # QC violin & scatter plots
├── doublets/                    # scDblFinder score plots
├── individual/                  # per-sample UMAP, markers, dot plots
├── integrated/                  # Harmony UMAP, heatmap, violin plots
├── annotation/                  # SingleR scores, marker feature plots
├── report_qc.pdf                # combined QC report
├── report_doublets.pdf
├── report_individual.pdf
├── report_annotation.pdf
└── report_integrated.pdf
```

Total runtime on the H1+H2 dataset (~1,200 cells): **~15–25 minutes** depending on CPU cores.

---

## Configuration

All parameters are centralised in `pipeline/config.R`. Open it in a text editor before running the pipeline.

---

### QC Thresholds

```r
QC <- list(
  min_features   = 200,    # min unique genes per cell
  max_features   = 5000,   # max unique genes per cell  (doublet proxy)
  min_counts     = 500,    # min total UMIs per cell
  max_counts     = 25000,  # max total UMIs per cell
  max_percent_mt = 20      # max mitochondrial gene %
)
```

**Recommended values by tissue:**

| Parameter | Human PBMC (recommended) | Non-PBMC tissue | Notes |
|-----------|--------------------------|-----------------|-------|
| `min_features` | 300–500 | 200–500 | 200 is permissive; raise to exclude low-quality cells |
| `max_features` | 4000–5000 | 6000–8000 | Cells above threshold are likely doublets |
| `min_counts` | 500–1000 | 500–1000 | Raise for high-depth libraries |
| `max_counts` | 15000–25000 | 20000–40000 | Generous upper cap |
| `max_percent_mt` | **10–15%** | 20–40% | **Current default (20%) is too permissive for PBMCs** — dying/stressed cells up to 20% will be retained. Lower to 10–15% for cleaner PBMC data. Neurons and cardiomyocytes naturally have higher MT% so keep higher for those. |

> **Tip:** Always inspect `report_qc.pdf` before finalising thresholds. The violin and scatter plots show where the natural breaks in your data are — adjust thresholds to those breaks rather than using fixed defaults.

---

### Cell Type Detection

The pipeline detects cell types via **two independent routes**:

**1. SingleR (automated, unrestricted)**
Scores every cell against the full `HumanPrimaryCellAtlasData` reference (~30+ types) using the entire transcriptome — not just canonical markers. Cell types detectable include:

- T cells (CD4, CD8, Tregs), NK cells, NKT cells
- B cells, Plasma cells
- Monocytes (CD14+, CD16+/FCGR3A+), Macrophages
- Dendritic cells (mDC, pDC), Mast cells
- Neutrophils, Eosinophils, Basophils *(if present — rare in healthy PBMC)*
- Megakaryocytes/Platelets, Erythrocytes
- Endothelial, Epithelial, Fibroblasts *(tissue samples)*

> Rare populations like Neutrophils will be detected and labelled by SingleR even if they are not in your canonical `MARKERS` list.

**2. Manual annotation via `CLUSTER_CELLTYPE_MAP`**
You assign any label you choose after inspecting the dot plot. The canonical markers are only used for visualisation — cluster discovery uses all 2,000 HVGs.

---

### Canonical Marker Genes

These are used only for **feature plots and dot plots** (visual QC). They do not affect clustering or SingleR annotation.

```r
MARKERS <- list(
  T_pan       = c("CD3D", "CD3E"),
  CD4_T       = c("CD4", "IL7R", "CCR7"),
  CD8_T       = c("CD8A", "CD8B", "GZMK"),
  NK          = c("NKG7", "GNLY", "KLRD1"),
  B_cell      = c("MS4A1", "CD79A", "CD19"),
  CD14_mono   = c("CD14", "LYZ", "CST3", "S100A8"),
  FCGR3A_mono = c("FCGR3A", "MS4A7"),
  DC          = c("FCER1A", "CLEC9A"),
  Platelet    = c("PPBP", "PF4")
)
```

The full marker list in `config.R` now includes all major populations plus contamination indicators:

| Group | Markers | Interpretation |
|-------|---------|----------------|
| CD4 T | CD4, IL7R, CCR7 | Expected ~40–50% of healthy PBMC |
| CD8 T | CD8A, CD8B, GZMK | Expected ~20–30% |
| **Treg** | FOXP3, IL2RA, CTLA4 | Small (~5–10% of T cells); expanded in cancer/autoimmunity |
| NK | NKG7, GNLY, KLRD1 | Expected ~5–15% |
| B cell | MS4A1, CD79A, CD19 | Expected ~5–10% |
| **Plasma** | MZB1, JCHAIN, SDC1 | Rare in healthy blood; elevated in infection or myeloma |
| CD14+ Mono | CD14, LYZ, S100A8 | Expected ~10–15% |
| FCGR3A+ Mono | FCGR3A, MS4A7 | Expected ~2–5% |
| **Neutrophil** | FCGR3B, CSF3R, CXCR2, CEACAM8 | Rare in healthy PBMC; **elevated = inflammation or poor/slow Ficoll separation** |
| DC | FCER1A, CLEC9A | Rare (~1%); absence after filtering suggests over-filtering |
| Platelet | PPBP, PF4 | Small cluster is normal; large cluster = **poor PBMC separation** |
| **RBC** | HBB, HBA1, HBA2, GYPA | **Contamination** — should be absent; present = failed Ficoll gradient |
| **HSPC** | CD34, GATA2, AVP | **Contamination** — haematopoietic progenitors; seen in mobilised or bone marrow samples |

**Quality indicators to check in `report_annotation.pdf`:**

- ✅ **Good:** T cells ~60–70%, NK ~10%, B cells ~5–10%, Monocytes ~15%
- ✅ **Good:** No RBC or HSPC clusters detected
- ⚠️ **Warning:** Neutrophil cluster present — donor may be inflamed, or sample processed slowly (neutrophils lyse quickly at room temperature)
- ⚠️ **Warning:** Large platelet cluster (>5%) — platelet aggregation during PBMC isolation
- ❌ **Contamination:** RBC cluster (HBB/HBA1+ cells) — Ficoll gradient failed; repeat isolation or filter HBB-high cells
- ❌ **Contamination:** HSPC cluster (CD34+) — sample is not pure PBMC; may be bone marrow aspirate or G-CSF mobilised blood

---

### Clustering Resolution

```r
CLUSTER <- list(
  resolutions = c(0.3, 0.4, 0.5, 0.6, 0.8),  # all tested; stored in metadata
  default_res = 0.5,                            # used for downstream analysis
  algorithm   = 1                               # 1 = Louvain, 4 = Leiden
)
```

| Resolution | Typical clusters (PBMC ~1000 cells) | Use when |
|-----------|--------------------------------------|----------|
| 0.3 | 8–10 | Starting point; broad populations |
| 0.5 | 12–15 | **Default — good for most PBMC runs** |
| 0.8 | 18–22 | Needed to split CD4/CD8 T or mono subtypes |

> Inspect the elbow plot in `report_individual.pdf` and the cluster UMAP at multiple resolutions before finalising.

---

### Manual Annotation Workflow

After the first full run:

1. Open `results_*/annotation/canonical_markers_dotplot.pdf`
2. For each cluster, identify which lineage markers are expressed (large bright dots)
3. Fill `CLUSTER_CELLTYPE_MAP` in `config.R`:
   ```r
   CLUSTER_CELLTYPE_MAP <- c(
     "0" = "CD4 T",
     "1" = "CD8 T",
     "2" = "NK",
     "3" = "B cell",
     "4" = "CD14+ Mono",
     "5" = "FCGR3A+ Mono",
     "6" = "DC",
     "7" = "Platelet"
   )
   ```
4. Re-run annotation and visualisation only:
   ```bash
   bash pipeline/run_pipeline.sh /path/to/S1 /path/to/S2 05 06 07
   ```

> If `CLUSTER_CELLTYPE_MAP` is `NULL`, SingleR majority labels are used automatically as a fallback.

---

## Adapting for Other Species or Tissues

| Parameter | Location | Example |
|-----------|----------|---------|
| MT gene pattern | `config.R` line ~34 | `"^MT-"` → `"^mt-"` (mouse) |
| `max_percent_mt` | `config.R` `QC` list | 20 → 35 (heart/brain tissue) |
| SingleR reference | `05_annotate.R` line ~36 | `HumanPrimaryCellAtlasData()` → `MouseRNAseqData()` or `ImmGenData()` |
| Canonical markers | `config.R` `MARKERS` list | Replace with tissue-specific genes |
| Cell type colours | `config.R` `CELLTYPE_COLORS` | Update names to match expected types |
| Harmony theta | `config.R` `HARMONY$theta` | 2 → 3–5 for strong batch effects |

---

## License

MIT License — see [LICENSE](LICENSE) for details.
