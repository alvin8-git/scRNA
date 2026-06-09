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
- **9-step modular pipeline** — each step saves intermediate `.rds` objects so you can re-run from any point
- **Per-sample RDS cache** — steps 01–03 write results to `sample_cache/<name>/`; subsequent runs with the same sample in a different integration combination skip reprocessing automatically
- **Parallel execution** — `future` multicore (16 GB worker budget) for Seurat steps, `BiocParallel` for scDblFinder and SingleR
- **Automated cell type annotation** — SingleR majority vote per cluster + scType marker-based scoring (no reference bias, gene-frequency weighted); suggested `CLUSTER_CELLTYPE_MAP` printed to log; partial manual overrides with automatic SingleR fallback for unmapped clusters
- **Contamination cell type detection** — rare/contaminating populations (Neutrophil, RBC, HSPC, Platelet, Basophil, Eosinophil, Mast cell) override cluster majority-vote labels via per-cell SingleR labels so they always appear on UMAP regardless of how many cells are in the cluster
- **Extended label normalisation** — 30-entry `SINGLER_NORM` maps all raw HumanPrimaryCellAtlasData labels (including granulocyte/erythroid lineages) to canonical pipeline names
- **Differential expression** — step 06b runs `FindMarkers` per cell type between samples; volcano plots, DE tables, and inflammatory/interferon/T-activation module score plots saved automatically (multi-sample only)
- **Combined PDF reports** — five numbered reports (`01-QC_report.pdf` … `05-Integrated_report.pdf`) plus an A4-normalised `Overall_report.pdf` (curated cross-stage summary); split-by-sample UMAP paginated at 2 panels per page with a shared cell-type legend
- **Comparison report** — step 08 (`08_comparison_report.R`) generates a standalone `Comparison_report.pdf` addressing sample quality, doublet rates, composition, integration, and DE in one document
- **Dynamic results directory** — named after input samples (e.g. `results_H1andH2/`)

---

## Documentation

| Document | Type | Description |
|----------|------|-------------|
| [Getting Started Tutorial](docs/tutorial-first-run.md) | Tutorial | Install, run, annotate — from raw matrices to a full report |
| [How to Re-run from a Specific Step](docs/howto-rerun-steps.md) | How-to | Restart the pipeline from step 03, 05, or 07 without reprocessing |
| [How to Override Cell Type Annotations](docs/howto-override-annotations.md) | How-to | Fill and correct `CLUSTER_CELLTYPE_MAP` after first run |
| [How to Run on Bat Whole Blood](docs/howto-bat-whole-blood.md) | How-to | Bat-specific overrides, expected proportions, common mislabellings |
| [How to Add or Replace a Sample](docs/howto-add-sample.md) | How-to | Add a new sample to an existing cohort using the sample cache |
| [Configuration Reference](docs/reference-config.md) | Reference | Every `config.R` parameter with type, default, and guidance |
| [Pipeline Steps Reference](docs/reference-pipeline-steps.md) | Reference | Inputs, outputs, and config keys for each of the 10 steps |
| [Output Files Reference](docs/reference-outputs.md) | Reference | Full directory layout and key file descriptions |
| [Pipeline Architecture](docs/explanation-architecture.md) | Explanation | Why the modular design, sample cache, and Harmony were chosen |
| [Annotation Strategy](docs/explanation-annotation.md) | Explanation | How SingleR, contamination overrides, and sub-type refinement interact |

---

## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 01 | `01_load_qc.R` | Load 10x matrices, compute QC metrics, filter low-quality cells |
| 02 | `02_doublets.R` | Doublet detection with scDblFinder |
| 03 | `03_individual.R` | Per-sample: normalize, HVG, PCA, UMAP, cluster, markers |
| 04 | `04_integrate.R` | Merge samples + Harmony batch correction (or single-sample UMAP) |
| 05 | `05_annotate.R` | SingleR automated annotation + sub-type refinement + canonical marker plots |
| 06 | `06_visualize.R` | Publication-quality UMAP, dot plots, heatmaps, violin plots |
| 06b | `06b_differential.R` | DE between samples per cell type; volcano plots; module scores (multi-sample only) |
| 07 | `07_finalize_reports.R` | Merge all PDFs into 5 final combined reports + Overall_report.pdf |
| 08 | `08_comparison_report.R` | Standalone cross-sample comparison report (quality, composition, integration, DE) |
| 09 | `09_bootstrap_proportions.R` | Bootstrap-normalised cell type proportions across samples; multinomial 95% CI; pairwise chi-squared tests |
| 10 | `10_rarefaction.R` | Rarefaction analysis — minimum capture depth for stable proportion estimates |

---

## Requirements

- **Linux** (tested on Ubuntu 20.04)
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://github.com/mamba-org/mamba)
- ~8 GB RAM minimum; 32–64 GB recommended for multi-sample runs with large cell counts (>20 K cells); `future_mem_gb` in `config.R` defaults to 16 GB per worker
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
`Seurat`, `harmony`, `scDblFinder`, `SingleR`, `celldex`, `BiocParallel`, `ggplot2`, `patchwork`, `cowplot`, `pheatmap`, `viridis`, `dplyr`, `pdftools`, `magick`

> `pdftools` and `magick` are required for the A4-normalised `Overall_report.pdf`. Without them the pipeline falls back to a simple PDF merge.

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

# Species flag (default: human) — applies bat-specific marker/QC/reference overrides
bash pipeline/run_pipeline.sh bat /path/to/ES03 /path/to/ES12
bash pipeline/run_pipeline.sh human /path/to/H1 /path/to/H2

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
├── qc/                          # QC violin & scatter plots, cell_fate.csv
├── doublets/                    # scDblFinder score plots
├── individual/                  # per-sample UMAP, markers, dot plots
├── integrated/                  # Harmony UMAP, heatmap, violin plots
├── annotation/                  # SingleR scores, marker feature plots
├── reports/
│   ├── 01-QC_report.pdf
│   ├── 02-Doublet_report.pdf
│   ├── 03-Individual_report.pdf
│   ├── 04-Annotation_report.pdf
│   ├── 05-Integrated_report.pdf
│   └── Overall_report.pdf       # A4-normalised curated summary
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
Scores every cell against a reference transcriptome. The reference is set by `SINGLER_REF` in `config.R`:

| Value | Reference | Best for |
|-------|-----------|---------|
| `"HumanPrimaryCellAtlas"` | HumanPrimaryCellAtlasData (default) | General human tissues |
| `"MonacoImmune"` | MonacoImmuneData | Blood/PBMC — resolves CD4/CD8/Treg/γδ T, monocyte subtypes, pDC/mDC (29 types) |

The bat species override sets `MonacoImmune` automatically. Cell types detectable include:

- T cells (CD4, CD8, Tregs), NK cells, NKT cells
- B cells, Plasma cells
- Monocytes (CD14+, CD16+/FCGR3A+), Macrophages
- Dendritic cells (mDC, pDC), Mast cells
- Neutrophils, Eosinophils, Basophils *(if present — rare in healthy PBMC)*
- Megakaryocytes/Platelets, Erythrocytes
- Endothelial, Epithelial, Fibroblasts *(tissue samples)*

> Rare populations like Neutrophils will be detected and labelled by SingleR even if they are not in your canonical `MARKERS` list.

**2. Auto-generated cluster map + optional manual override via `CLUSTER_CELLTYPE_MAP`**
When `CLUSTER_CELLTYPE_MAP = NULL` (the default), the pipeline computes SingleR majority vote per cluster and uses that as the cell type label. A copy-pasteable `CLUSTER_CELLTYPE_MAP` is printed to `logs/05_annotate.log` after every run so you can review and correct it. Partial entries are supported — listed clusters get your label; any cluster not in the map falls back to SingleR automatically.

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
  algorithm   = 1,                              # 1 = Louvain, 4 = Leiden
  compare_res = c(0.5, 0.6, 0.8)               # side-by-side UMAP comparison saved to integrated/
)
```

| Resolution | Typical clusters (PBMC ~1000 cells) | Use when |
|-----------|--------------------------------------|----------|
| 0.3 | 8–10 | Starting point; broad populations |
| 0.5 | 12–15 | **Default — good for most PBMC runs** |
| 0.8 | 18–22 | Needed to split CD4/CD8 T or mono subtypes |

> Step 04 saves `cluster_resolution_comparison.pdf` showing `compare_res` side-by-side so you can pick the best resolution before re-running annotation.

> Inspect the elbow plot in `03-Individual_report.pdf` and the cluster UMAP at multiple resolutions before finalising.

---

### T Cell Sub-clustering

```r
SUBCLUSTER <- list(
  enabled    = TRUE,
  t_patterns = "T[_ ]cell|T cell|CD4|CD8|Treg|cytotox",  # regex to identify T cell clusters
  resolution = 0.8,    # sub-cluster resolution (start lower: 0.6–0.8 to avoid over-splitting)
  min_cells  = 20      # skip sub-clustering if fewer cells
)
```

Step 05 (`05_annotate.R`) automatically identifies clusters whose majority `final_cell_type` matches `t_patterns` and runs `FindSubCluster()` on them. Outputs saved to `annotation/`:

| File | Description |
|------|-------------|
| `tcell_subclusters_umap.pdf` | UMAP coloured by T cell sub-cluster identity |
| `tcell_subclusters_dotplot.pdf` | CD4/CD8/Treg marker dot plot per sub-cluster |
| `tcell_subcluster_summary.csv` | Sub-cluster cell counts and parent cluster mapping |

> At `resolution = 1.2`, T cell clusters may split into 15–18 sub-clusters (over-split). Start with `resolution = 0.6–0.8` for 4–8 meaningful sub-types (CD4 naive, CD4 memory, CD8 effector, Treg).

---

### Annotation Workflow

**Step 1 — First run (automatic, no config needed)**

Leave `CLUSTER_CELLTYPE_MAP = NULL` in `config.R`. The pipeline automatically assigns SingleR majority labels per cluster. After step 05 completes, the log prints a suggested map:

```
logs/05_annotate.log:
  CLUSTER_CELLTYPE_MAP <- c(
    "0" = "CD4 T (naive)",
    "1" = "NK",
    ...
  )
```

**Step 2 — Review (optional)**

1. Open `results_*/annotation/canonical_markers_dotplot.pdf`
2. Cross-check suggested labels against marker expression
3. Correct any wrong entries (SingleR occasionally mislabels DCs as Monocytes, neutrophils as Pre-B, etc.)

**Step 3 — Apply corrections (only if needed)**

Copy the suggested map from the log into `config.R`, editing any wrong labels:

```r
CLUSTER_CELLTYPE_MAP <- c(
  "0" = "CD4 T (naive)",
  "3" = "NK",           # override: SingleR said T_cells
  "11" = "DC"           # override: SingleR said Monocyte
  # all other clusters fall back to SingleR automatically
)
```

**Step 4 — Re-run**

```bash
bash pipeline/run_pipeline.sh /path/to/S1 /path/to/S2 05 06 07
```

> **Important:** Cluster numbers change between datasets. Never copy a `CLUSTER_CELLTYPE_MAP` from one sample run to another — always start from the log of that specific run.

---

## Adapting for Other Species or Tissues

| Parameter | Location | Example |
|-----------|----------|---------|
| Species keyword | `run_pipeline.sh` first arg | `bat` or `human` (default) |
| MT gene pattern | `config.R` line ~34 | `"^MT-"` → `"^mt-"` (mouse) |
| `max_percent_mt` | `config.R` `QC` list | 20 → 35 (heart/brain tissue) |
| SingleR reference | `config.R` `SINGLER_REF` | `"HumanPrimaryCellAtlas"` (default) → `"MonacoImmune"` (blood) |
| Canonical markers | `config.R` `MARKERS` list | Replace with tissue-specific genes |
| Contamination types | `config.R` `CONTAMINATION_TYPES` | Remove RBC/Neutrophil for whole-blood samples |
| Cell type colours | `config.R` `CELLTYPE_COLORS` | Update names to match expected types |
| Harmony theta | `config.R` `HARMONY$theta` | 2 → 3–5 for strong batch effects |
| Clustering resolution | `config.R` `CLUSTER$default_res` | 0.5 (PBMC) → 1.0 (whole blood / complex tissue) |

The `bat` keyword applies all *Eonycteris spelaea* whole-blood overrides automatically (MonacoImmune reference, res=1.0, γδ T markers, adjusted contamination list, bat-specific subtype markers).

---

## License

MIT License — see [LICENSE](LICENSE) for details.
