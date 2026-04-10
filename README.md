# PBMC scRNA-seq Analysis Pipeline

Single-cell RNA-seq analysis of human PBMC samples (H1, H2) from DNB C4 sequencing platform, using Seurat v5 and Harmony integration.

## Quick Start

```bash
# 1. Activate environment
conda activate scrna_seurat

# 2. Run full pipeline
bash /data/alvin/scRNA/pipeline/run_pipeline.sh

# 3. After inspecting annotation plots, fill in cell type map in config.R, then:
bash /data/alvin/scRNA/pipeline/run_pipeline.sh 05 06
```

## Project Structure

```
/data/alvin/scRNA/
├── H1/                          # DNB C4 output — Donor 1
│   ├── filter_matrix/           # Pre-filtered expression matrix (386 cells)
│   └── raw_matrix/              # All barcodes (unfiltered)
├── H2/                          # DNB C4 output — Donor 2
│   ├── filter_matrix/           # Pre-filtered expression matrix (967 cells)
│   └── raw_matrix/
├── pipeline/                    # Analysis scripts
│   ├── config.R                 # Central configuration (edit this)
│   ├── 01_load_qc.R
│   ├── 02_doublets.R
│   ├── 03_individual.R
│   ├── 04_integrate.R
│   ├── 05_annotate.R
│   ├── 06_visualize.R
│   ├── run_pipeline.sh
│   └── setup_env.sh
├── results/                     # All outputs (generated)
│   ├── qc/
│   ├── doublets/
│   ├── individual/
│   ├── integrated/
│   └── annotation/
└── 20260401_pbmc.pptx           # Reference slides
```

## Dataset

| Sample | Donor | Cells (raw) | Cells (post-QC) | Platform |
|--------|-------|-------------|-----------------|----------|
| H1 | Human PBMC | 386 | 359 | DNB C4 / dnbc4tools v3.0 |
| H2 | Human PBMC | 967 | 848 | DNB C4 / dnbc4tools v3.0 |

**Total after QC + doublet removal:** 1,207 cells across 7 clusters

## Cell Type Composition (SingleR, current run)

| Cell Type | Count | % |
|-----------|-------|---|
| T cells | 720 | 59.7% |
| NK cells | 250 | 20.7% |
| B cells | 117 | 9.7% |
| Monocytes | 109 | 9.0% |
| Platelets | 11 | 0.9% |

> Note: SingleR labels are used until `CLUSTER_CELLTYPE_MAP` is filled in `config.R`.
> T cells currently includes both CD4+ and CD8+ — manual annotation will split these.

## Environment

- **Conda env:** `scrna_seurat`
- **R:** 4.3.3
- **Seurat:** v5.1.0
- **Key packages:** harmony, scDblFinder, SingleR, celldex, ggplot2, patchwork
