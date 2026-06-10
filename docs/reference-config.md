# Pipeline Configuration Reference

All pipeline behaviour is controlled by `pipeline/config.R`. Every R script sources this file at startup. Edit it before running; all parameter values documented here are the shipped defaults.

---

## Config validation (`validate_config.R`)

`pipeline/validate_config.R` runs automatically before step 01 (called by `run_pipeline.sh`). It exits with code 1 and a clear error message on any of the following:

| Check | What it catches |
|-------|----------------|
| CELLTYPE_COLORS coverage | A type in `CLUSTER_CELLTYPE_MAP` has no colour entry — plots would silently get grey |
| CLUSTER_CELLTYPE_MAP key format | Keys that are bare integers instead of quoted strings (`0` vs `"0"`) — the map would silently not match |
| SAMPLE_PATHS existence | A path in `SAMPLE_PATHS` doesn't exist on disk |

You can also run it standalone: `Rscript pipeline/validate_config.R`. It uses a `commandArgs()` path fallback so it works whether sourced or invoked directly without editing hard-coded paths.

> **Pending:** Check 3 (SAMPLE_PATHS) fails if data lives on a NAS or drive that isn't mounted yet. A `--skip-paths` flag is tracked in TODO.md.

---

## Base paths

```r
BASE_DIR        <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
SAMPLE_CACHE_DIR <- file.path(BASE_DIR, "sample_cache")
```

`BASE_DIR` resolves to the repository root. Override it without touching `config.R` by setting the `SCRNA_BASE_DIR` env var:

```bash
SCRNA_BASE_DIR=/mnt/nas/project bash pipeline/run_pipeline.sh /path/to/SampleA
```

`SAMPLE_CACHE_DIR` is shared across all sample combinations — deleting a subdirectory forces that sample to be reprocessed from step 01.

---

## Sample resolution

Samples are set via environment variables (injected by `run_pipeline.sh`) or hardcoded fallbacks.

| Env var | Type | Description |
|---------|------|-------------|
| `SCRNA_SAMPLE1` … `SCRNA_SAMPLEN` | path | Absolute paths to sample folders in order |
| `SCRNA_SPECIES` | string | `human` (default) or `bat` or `bat_wing` |
| `SCRNA_CONDITION` | string | Comma-separated `name=label` pairs for DEG grouping |
| `SCRNA_BASE_DIR` | path | Overrides `BASE_DIR` — point outputs at a different root without editing `config.R` |

Example:
```bash
export SCRNA_SAMPLE1=/data/H1
export SCRNA_SAMPLE2=/data/H2
export SCRNA_SPECIES=human
export SCRNA_CONDITION="H1=control,H2=treated"
```

Harmony integration runs automatically when more than one sample is provided.

---

## QC

```r
QC <- list(
  min_features   = 200,   # min unique genes per cell
  max_features   = 5000,  # max unique genes per cell (doublet proxy)
  min_counts     = 500,   # min total UMIs per cell
  max_counts     = 25000, # max total UMIs per cell
  max_percent_mt = 20     # max mitochondrial gene %
)
```

**Guidance by tissue:**

| Parameter | Human PBMC | Whole blood | Non-PBMC tissue |
|-----------|-----------|-------------|-----------------|
| `min_features` | 300–500 | 200–300 | 200–500 |
| `max_features` | 4000–5000 | 4000–6000 | 6000–8000 |
| `min_counts` | 500–1000 | 500 | 500–1000 |
| `max_percent_mt` | **10–15%** | 20% | 20–40% |

The default `max_percent_mt = 20` is permissive for PBMC; lower to 10–15% for cleaner lymphocyte data.

---

## Doublet detection

```r
DOUBLET <- list(
  db_rate = NULL  # NULL = auto-estimate from cell count (recommended)
)
```

`NULL` lets scDblFinder estimate the expected doublet rate from the number of recovered cells (~0.8% per 1,000 cells). Set a numeric value (e.g. `0.05`) to fix the rate.

---

## Normalisation & HVG

```r
NORM <- list(
  method    = "LogNormalize",
  scale_fac = 10000,
  n_hvg     = 2000
)
```

`n_hvg`: number of highly variable genes used for PCA. 2,000 is appropriate for most PBMC runs; raise to 3,000–5,000 for complex tissues.

---

## Dimensionality reduction

```r
DIMS <- list(
  pca_dims = 50,   # PCs computed
  umap_dims = 1:20 # PCs fed into UMAP and clustering
)
```

Inspect the elbow plot in `03-Individual_report.pdf` to confirm `umap_dims` captures most variance. Typical PBMC: PCs 1–15; whole blood with granulocytes: PCs 1–20.

---

## Clustering

```r
CLUSTER <- list(
  resolutions = c(0.3, 0.4, 0.5, 0.6, 0.8),
  default_res = 0.5,
  compare_res = c(0.5, 0.6, 0.8),
  algorithm   = 1                             # 1 = Louvain, 4 = Leiden
)
```

All resolutions in `resolutions` are computed and stored in cell metadata. Only `default_res` is used downstream. `compare_res` controls the side-by-side comparison UMAP saved in `integrated/`.

| Resolution | Typical clusters (PBMC ~1 K cells) |
|-----------|--------------------------------------|
| 0.3 | 8–10 |
| 0.5 | 12–15 (default) |
| 0.8 | 18–22 |

---

## T cell sub-clustering

```r
SUBCLUSTER <- list(
  enabled    = TRUE,
  t_patterns = "T[_ ]cell|T cell|CD4|CD8|Treg|cytotox",
  resolution = 0.8,
  min_cells  = 20
)
```

`t_patterns`: regex matched against `final_cell_type`; clusters where the majority label matches are sub-clustered. Set `enabled = FALSE` to skip entirely. `min_cells`: clusters smaller than this are skipped.

---

## Harmony integration

```r
HARMONY <- list(
  group_by_vars = "sample",
  theta         = 2,
  lambda        = 1,
  nclust        = 50,
  max_iter      = 20,
  dims_use      = 1:20
)
```

| Parameter | Effect | When to change |
|-----------|--------|----------------|
| `theta` | Diversity penalty. Higher = stronger batch correction | Raise to 3–5 for strong batch effects (different sequencing runs) |
| `lambda` | Ridge regression penalty | Rarely changed |
| `nclust` | Harmony soft-clustering centroids | Raise for >10 samples |
| `dims_use` | PCs fed into Harmony | Match `DIMS$umap_dims` |

---

## SingleR reference

```r
SINGLER_REF <- "HumanPrimaryCellAtlas"
```

| Value | Reference | Best for |
|-------|-----------|---------|
| `"HumanPrimaryCellAtlas"` | HumanPrimaryCellAtlasData | Broad human tissues (default) |
| `"MonacoImmune"` | MonacoImmuneData | Blood/PBMC — resolves CD4/CD8/Treg/γδ T, monocyte subtypes, pDC/mDC (29 types) |

The `bat` species override automatically sets `MonacoImmune`.

---

## Canonical marker genes (`MARKERS`)

`MARKERS` is a named list of character vectors. Keys are cell-type names matching `CELLTYPE_COLORS`; values are HGNC gene symbols. Used by step 05 for marker dot plots and by step 06 for feature plots.

Example subset:
```r
MARKERS <- list(
  "T cell"     = c("CD3D", "CD3E", "CD3G"),
  "CD4 T"      = c("CD3D", "CD3E", "CD4", "IL7R"),
  "CD8 T"      = c("CD3D", "CD3E", "CD8A", "CD8B"),
  "NK"         = c("GNLY", "NKG7", "KLRD1", "GZMB", "NCAM1"),
  "B cell"     = c("MS4A1", "CD79A", "CD79B", "TCL1A"),
  "CD14+ Mono" = c("CD14", "LYZ", "S100A8", "S100A9"),
  "Neutrophil" = c("S100A12", "S100A9", "BST1", "G0S2", "FCGR3B")
)
```

For bat whole blood, the `bat` species keyword replaces these with *Eonycteris spelaea*-validated orthologues.

### `MARKERS$compute_integrated`

```r
MARKERS <- list(
  ...
  compute_integrated = FALSE   # set TRUE to run FindAllMarkers after step 04 integration
)
```

When `FALSE` (default), step 04 skips `FindAllMarkers` entirely. The full marker sweep on a 20 K-cell integrated dataset takes 20–30 minutes and isn't needed for annotation. Set to `TRUE` to write `integrated/integrated_cluster_markers.csv`. This flag was added because the previous default (always run) made step 04 the slowest step for no benefit on most runs.

---

## Sub-type refinement markers (`REFINEMENT_MARKERS`)

Used by step 05 to split coarse SingleR labels (e.g. `"CD4 T"`) into sub-types (`"CD4 T (naive)"`, `"CD4 T (memory)"`, `"CD4 T (effector)"`).

```r
REFINEMENT_MARKERS <- list(
  "CD4 T" = list(
    "CD4 T (naive)"    = c("CCR7", "TCF7", "LEF1", "SELL"),
    "CD4 T (memory)"   = c("IL7R", "AQP3", "GPR183"),
    "CD4 T (effector)" = c("GZMK", "TNFRSF4", "CCL5")
  ),
  "CD8 T" = list(
    "CD8 T (naive)"    = c("CCR7", "TCF7", "LEF1"),
    "CD8 T (effector)" = c("GZMB", "PRF1", "GNLY", "IFNG")
  )
)
```

Set to `NULL` to disable sub-type refinement and keep coarse labels.

---

## Manual cluster annotation (`CLUSTER_CELLTYPE_MAP`)

```r
CLUSTER_CELLTYPE_MAP <- NULL
```

`NULL` = auto-annotate using SingleR majority vote per cluster (recommended first run). After step 05 runs, the log file `logs/05_annotate.log` prints a copy-pasteable map. Edit wrong entries and paste into `config.R`, then re-run steps 05–07.

Partial maps are supported: listed clusters get your label; any cluster not in the map falls back to SingleR.

> Cluster numbers change between datasets. Never copy a `CLUSTER_CELLTYPE_MAP` from one sample combination to another.

---

## Contamination types (`CONTAMINATION_TYPES`)

```r
CONTAMINATION_TYPES <- c("Neutrophil", "RBC", "HSPC", "Platelet",
                          "Basophil", "Eosinophil", "Mast cell")
```

Cell types in this list bypass majority-vote cluster labelling — each cell retains its per-cell SingleR label. This ensures rare populations always appear on UMAP regardless of cluster size. Remove types that are biologically expected (e.g. remove `"Neutrophil"` for whole-blood samples where neutrophils are the majority).

---

## Color palettes

```r
SAMPLE_COLORS    <- c(H1 = "#E64B35", H2 = "#4DBBD5")
CELLTYPE_COLORS  <- c("CD4 T" = "#E64B35", "NK" = "#00A087", ...)
```

`SAMPLE_COLORS`: add an entry for each new sample name. The pipeline auto-assigns `hue_pal()` colours for samples not in this list.

`CELLTYPE_COLORS`: keys must match `final_cell_type` labels in the annotated Seurat object exactly.

---

## Parallelism

```r
WORKERS       <- min(parallel::detectCores() - 2, 8)
future_mem_gb <- 16
```

`WORKERS`: capped at 8 to avoid memory exhaustion. `future_mem_gb`: maximum RAM allocated per `future` worker (16 GB). Raise for large datasets (>50 K cells); lower on machines with <32 GB RAM.

---

## Plot defaults

```r
PLOT <- list(
  width      = 8,
  height     = 7,
  dpi        = 300,
  pt_size    = 0.8,
  label_size = 4
)
```

`pt_size`: point size in all UMAPs. Raise for sparse datasets (<5 K cells); lower for dense datasets (>100 K cells). `label_size`: cluster label font size on UMAPs.

---

## Species overrides

Set via the `SCRNA_SPECIES` env var (injected by `run_pipeline.sh`):

| Value | Effect |
|-------|--------|
| `"human"` (default) | Standard PBMC settings |
| `"bat"` | MonacoImmune reference, res=1.0, γδ T markers, bat-specific REFINEMENT_MARKERS, adjusted CONTAMINATION_TYPES |
| `"bat_wing"` | Bat species with wing-tissue QC thresholds (higher `max_percent_mt`, lower `min_counts`) |

---

## Related

- [Pipeline Steps Reference](reference-pipeline-steps.md) — what each step reads and writes
- [How to Override Annotations](howto-override-annotations.md) — step-by-step CLUSTER_CELLTYPE_MAP workflow
- [How to Run on Bat Whole Blood](howto-bat-whole-blood.md) — bat-specific parameter guidance
