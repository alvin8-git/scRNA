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

All parameters are in `pipeline/config.R`. Key sections to edit:

```r
# QC thresholds
QC <- list(
  min_features    = 200,
  max_features    = 6000,
  max_percent_mt  = 20
)

# Clustering resolution
CLUSTER <- list(default_res = 0.5, resolutions = c(0.3, 0.5, 0.8))

# Manual cell type annotation (fill after inspecting annotation/canonical_markers_dotplot.pdf)
CLUSTER_CELLTYPE_MAP <- NULL   # e.g. c("0"="CD4 T", "1"="CD8 T", ...)
```

### Manual annotation workflow

After the first full run, inspect `results_*/annotation/canonical_markers_dotplot.pdf`, then:

1. Fill `CLUSTER_CELLTYPE_MAP` in `config.R`
2. Re-run annotation and visualization only:
   ```bash
   bash pipeline/run_pipeline.sh /path/to/S1 /path/to/S2 05 06 07
   ```

---

## Adapting for other species or tissues

| Parameter | Location | Example change |
|-----------|----------|----------------|
| MT gene pattern | `config.R` `QC$mt_pattern` | `"^MT-"` → `"^mt-"` (mouse) |
| SingleR reference | `05_annotate.R` line ~36 | `HumanPrimaryCellAtlasData()` → `MouseRNAseqData()` |
| Canonical markers | `config.R` `MARKERS` list | Replace with tissue-specific genes |
| Cell type colours | `config.R` `CELLTYPE_COLORS` | Update to match expected types |

---

## License

MIT License — see [LICENSE](LICENSE) for details.
