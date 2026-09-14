# Pipeline Configuration Reference

All pipeline behaviour is controlled by `pipeline/config.R`. Every R script sources this file at startup. Edit it before running; all parameter values documented here are the shipped defaults.

---

## Config validation (`validate_config.R`)

`pipeline/validate_config.R` runs automatically before step 01 (called by `run_pipeline.sh`) and can be run standalone: `Rscript pipeline/validate_config.R`. It resolves `config.R` relative to its own location (`sys.frame`/`commandArgs` — no hard-coded path), so it works whether sourced or invoked directly.

| Check | What it catches | On failure |
|-------|----------------|-----------|
| CELLTYPE_COLORS coverage | A type in `CLUSTER_CELLTYPE_MAP` has no colour entry — plots would silently get grey | **error** (exit 1) |
| CLUSTER_CELLTYPE_MAP key format | Keys that are not cluster numbers — the map would silently not match | **error** (exit 1) |
| SAMPLE_PATHS existence | A path in `SAMPLE_PATHS` doesn't exist on disk | **warning** (continues) |

`SAMPLE_PATHS` existence is a **warning by default**, not a hard error — data often lives on a NAS or external drive that isn't mounted at validation time, and config-only checks should still pass. Pass `--strict-paths` to escalate missing paths to an error (e.g. in CI):

```bash
Rscript pipeline/validate_config.R --strict-paths
```

Source-level regression checks (e.g. "step 10 must not hardcode a ground-truth sample") live in `pipeline/tests/test_regressions.R`, not the validator — `validate_config.R` covers config invariants only.

---

## Base paths

```r
BASE_DIR         <- Sys.getenv("SCRNA_BASE_DIR", "/data/alvin/scRNA")
SAMPLE_CACHE_DIR <- file.path(BASE_DIR, "sample_cache")
```

`BASE_DIR` resolves to the repository root. Override it without touching `config.R` by setting the `SCRNA_BASE_DIR` env var:

```bash
SCRNA_BASE_DIR=/mnt/nas/project bash pipeline/run_pipeline.sh /path/to/SampleA
```

`SAMPLE_CACHE_DIR` is shared across all sample combinations — deleting a subdirectory forces that sample to be reprocessed from step 01.

---

## Cache invalidation engine

Steps 01–03 cache their per-sample output `.rds` under `sample_cache/<name>/`, each paired with a `.hash` sidecar. Before recomputing, a step compares the stored hash against `cache_hash(nm, step)`:

```
cache_hash(nm, step) = md5( cumulative_params(step) + resolved_input_path + matrix_fingerprint(path) )
```

- **cumulative_params** — per-step and nested: step 01 hashes `QC` + species; 02 adds `DOUBLET`; 03 adds `NORM`/`DIM`/`CLUSTER`. Because `digest` hashes the whole nested list, changing a threshold vector or array parameter (e.g. `QC$max_features`, `CLUSTER$resolutions`) changes the digest and invalidates **only** the steps that depend on it — a clustering change does not bust the QC (01) or doublet (02) caches.
- **resolved_input_path** — the absolute sample path, so two experiments that share a folder name (e.g. `PBMC/`) never collide in `sample_cache/`.
- **matrix_fingerprint** — `size + mtime` of the 10x matrix files, so editing the input matrix busts the cache even when config is byte-identical.

A mismatch logs `[CACHE STALE]` and recomputes once (linear; no retry loop or halt); a match logs `[CACHE HIT]` and copies the cached object. Parameters that only affect steps 04+ (`HARMONY`, `MARKERS`, `SINGLER_REF`, `CLUSTER_CELLTYPE_MAP`) are deliberately excluded from the key — those steps are not cached and recompute every run.

---

## Sample resolution

Samples are set via environment variables (injected by `run_pipeline.sh`) or hardcoded fallbacks.

| Env var | Type | Description |
|---------|------|-------------|
| `SCRNA_SAMPLE1` … `SCRNA_SAMPLEN` | path | Absolute paths to sample folders in order |
| `SCRNA_SPECIES` | string | `human` (default), `bat`, `bat_wing`, or `cm` |
| `SCRNA_CONDITION` | string | Comma-separated `name=label` pairs for DEG grouping |
| `SCRNA_BASE_DIR` | path | Overrides `BASE_DIR` — point outputs at a different root without editing `config.R` |
| `SCRNA_RESULTS_DIR` | path | Point a step at an existing run directory instead of deriving one from the sample list. Also suppresses creation of per-sample `individual/` subdirectories |
| `SCRNA_REFERENCE_MODEL` | path | Frozen SingleR model built by `build_reference.R`. Sets `REFERENCE_MODEL`; empty (default) makes steps 05r and 08c self-skip |
| `SCRNA_ANCHORS` | string | Comma-separated anchor sample names for the cross-run benchmark. Sets `ANCHOR_SAMPLES`, default `Aksh1,ES332` |
| `SCRNA_DRIFT_PP` | number | Anchor drift threshold in percentage points. Sets `DRIFT_FLAG_PP`, default `5` |

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

**`bat_wing` raises two of these.** The overlay sets `max_features = 8000` and `max_counts = 60000`,
because keratinocytes and fibroblasts are far larger and more transcriptionally complex than
lymphocytes — the blood caps discard 5–11% of genuine tissue cells. `max_percent_mt` stays at 20:
this platform runs a ~0.5–1% mitochondrial baseline (median 0.80% across 40 samples, human PBMC
included), so 20% is already a near-inert filter and raising it buys nothing.

---

## Doublet detection

```r
DOUBLET <- list(
  doublet_rate = list(H1 = 0.031, H2 = 0.077),  # NULL for unknown samples → auto
  PCs  = 1:15,
  sct  = FALSE
)
```

`doublet_rate` is a per-sample named list of expected doublet rates — the shipped default covers the
bundled `H1`/`H2` example data only, so add an entry per sample or set the whole field to `NULL`.
`NULL` lets scDblFinder estimate the rate from the number of recovered cells (~0.8% per 1,000 cells).
`PCs`: principal components fed to scDblFinder. `sct`: use SCTransform normalisation instead of
LogNormalize for the doublet call.

---

## Normalisation & HVG

```r
NORM <- list(
  method       = "LogNormalize",
  scale_factor = 10000,
  n_hvg        = 2000,
  hvg_method   = "vst"
)
```

`n_hvg`: number of highly variable genes used for PCA. 2,000 is appropriate for most PBMC runs; raise to 3,000–5,000 for complex tissues.

---

## Dimensionality reduction

```r
DIM <- list(
  npcs      = 30,    # PCs computed
  dims_use  = 1:20,  # PCs fed into UMAP and clustering
  umap_seed = 42     # seed for reproducible UMAP embeddings
)
```

Inspect the elbow plot in `03-Individual_report.pdf` to confirm `dims_use` captures most variance. Typical PBMC: PCs 1–15; whole blood with granulocytes: PCs 1–20.

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

`t_patterns`: regex matched against `cell_type`; clusters where the majority label matches are sub-clustered. Set `enabled = FALSE` to skip entirely. `min_cells`: clusters smaller than this are skipped.

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
| `dims_use` | PCs fed into Harmony | Match `DIM$dims_use` |

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

## Sub-type refinement markers (`SUBTYPE_MARKERS`)

Used by step 05 to split coarse SingleR labels (e.g. `"CD4 T"`) into sub-types (`"CD4 T (naive)"`, `"CD4 T (memory)"`, `"CD4 T (effector)"`).

```r
SUBTYPE_MARKERS <- list(
  "CD4 T" = list(
    "CD4 T (naive)"    = c("CCR7", "SELL", "TCF7", "LEF1"),
    "CD4 T (effector)" = c("GZMK", "GZMB", "TNFRSF4", "PRF1"),
    "CD4 T (memory)"   = c("IL7R", "AQP3", "GPR183", "S100A4")
  ),
  "B cell" = list(
    "B cell (naive)"   = c("IGHD", "IGHM", "TCL1A", "IL4R"),
    "B cell (memory)"  = c("IGHG1", "IGHG2", "IGHA1", "TNFRSF13B"),
    "Plasma"           = c("MZB1", "JCHAIN", "SDC1", "CD38", "XBP1", "PRDM1")
  ),
  "Monocyte" = list(
    "CD14+ Mono"       = c("CD14", "S100A8", "S100A9", "LYZ"),
    "FCGR3A+ Mono"     = c("FCGR3A", "CDKN1C", "MS4A7")
  )
)
```

Each top-level key must match a label SingleR can produce *after* `SINGLER_NORM` normalisation.
Scoring is the average expression of the listed genes per cluster; the highest-scoring sub-type wins.
Set to `NULL` to disable sub-type refinement and keep coarse labels.

---

## Manual cluster annotation (`CLUSTER_CELLTYPE_MAP`)

```r
CLUSTER_CELLTYPE_MAP <- NULL
```

`NULL` = auto-annotate using SingleR majority vote per cluster (recommended first run). After step 05 runs, the log file `logs/05_annotate.log` prints a copy-pasteable map. Edit wrong entries and paste into `config.R`, then re-run steps 05–07.

**Map every cluster, not just the wrong ones.** Clusters absent from the map fall back to the
*per-cell* SingleR label rather than the cluster majority (`05_annotate.R:379`), which fragments
each unmapped cluster into a mix of labels. A partial map is worse than no map.

**Guard a live map by sample set.** Cluster ids only mean something for the exact sample set and
integration that produced them, so wrap any committed map in a guard rather than assigning it at
top level:

```r
if (length(SAMPLE_NAMES) == 3 && setequal(SAMPLE_NAMES, c("T1", "T2", "T6"))) {
  CLUSTER_CELLTYPE_MAP <- c("0" = "Fibroblast", "1" = "RBC", ...)
}
```

Adding or removing a sample silently disables the map and reverts to auto-annotation, which is the
safe failure. Live maps: bat wing T1/T2/T6 in `config.R` (added because `HumanPrimaryCellAtlas`
swaps fibroblast and smooth muscle in wing tissue), and the 7-sample H1 cardio cohort in
`config_species_cm.R`; everything else is `NULL`.

**Guard on `QC` too when the floor varies.** The cardio maps additionally test
`isTRUE(QC$min_features == 1000)` (and `== 200` for the preserved earlier run), because the same
samples cluster differently at a different QC floor. A `QC`-guarded map must live in the species
overlay, not in `config.R`: `config.R` evaluates its map blocks before sourcing the overlays, so it
would read the base `QC` value (regression guard T8). Use `isTRUE(x == n)`, not
`identical(x, nL)` — `QC` values are doubles.

### Splitting one mixed cluster (`CLUSTER_SUBCLUSTER_MAP`)

```r
CLUSTER_SUBCLUSTER_MAP <- list(cluster = "1", resolution = 0.2,
  labels = c("0" = "Pluripotent", "1" = "Fibroblast", "2" = "Fibroblast", "3" = "Fibroblast", "4" = "Fibroblast"))
```

Optional; unset by default. Applied in step 05 right after `CLUSTER_CELLTYPE_MAP`: `FindSubCluster`
re-clusters `cluster` on the SNN graph at `resolution`, stores the result in the
`annot_subcluster` metadata column, and relabels that cluster's cells per subcluster. The cluster
still needs an entry in `CLUSTER_CELLTYPE_MAP` (a placeholder is fine). A subcluster missing from
`labels` stops step 05. Subcluster ids depend on the resolution and the graph, so re-curate after any
upstream change. Keep it inside the same guard as the map it refines.

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

`CELLTYPE_COLORS`: keys must match `cell_type` labels in the annotated Seurat object exactly.

---

## Parallelism

```r
PARALLEL <- list(
  workers       = ...,   # per-sample phase (01-03): cores-2, capped at 8
  merge_workers = ...,   # merged-object phase (04-06b): fewer (larger per-worker budget)
  future_mem_gb = 8L,    # per-worker RAM budget, per-sample phase
  merge_mem_gb  = 16L    # per-worker RAM budget, merged-object phase
)
```

Two worker budgets, both RAM-governed. `workers` drives the per-sample phase (steps 01–03), where objects are small (~8 GB/worker). `merge_workers` drives the merged-object phase (steps 04, 05, 06, 06b), where each `future`/`BiocParallel` worker can hold a copy of the full merged object, so it uses a larger per-worker budget (16 GB) and therefore fans out to **fewer** workers. A RAM governor reads `MemAvailable` and caps each count so `workers × budget` leaves ~20% headroom, preventing OOM in the fan-out phases. Both are also capped at `cores − 2` (max 8). OMP/OpenBLAS/MKL/BLAS threads are pinned to 1 per process so worker-spawned BLAS pools don't oversubscribe the CPU.

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
| `"bat"` | MonacoImmune reference, res=1.0, γδ T markers, bat-specific SUBTYPE_MARKERS, adjusted CONTAMINATION_TYPES |
| `"bat_wing"` | Wing-tissue markers (fibroblast, keratinocyte, endothelial, pericyte, macrophage, melanocyte), `HumanPrimaryCellAtlas` reference, res=0.5, tissue `SUBTYPE_MARKERS`, `WOUND_MODULES` for step 11, `CONTAMINATION_TYPES` reduced to RBC/HSPC/Platelet |
| `"cm"` | hESC/iPSC → cardiomyocyte differentiation (`config_species_cm.R`). `MARKERS` replaced by lineage panels (Cardiomyocyte, CM ventricular/atrial, Cardiac progenitor, Pluripotent, Fibroblast, Myofibroblast, Smooth Muscle, Pericyte, Endothelial, Epithelial, Epicardial, Hepatic/Endoderm, Proliferating); `HumanPrimaryCellAtlas` reference, which has **no cardiomyocyte label**, so a curated `CLUSTER_CELLTYPE_MAP` is expected; res 0.1–0.8, default 0.5; `CONTAMINATION_TYPES` Endothelial/Pericyte; `PARALLEL$merge_mem_gb` 40 |

The `bat` overlay does not change `QC`. The `bat_wing` overlay raises `QC$max_features` to 8000 and
`QC$max_counts` to 60000 (see [QC](#qc)). The `cm` overlay sets `QC$max_features` 9000,
`QC$max_counts` 70000 and `QC$min_features` 1000. The 1000 floor drops cells unevenly by day
(D0 20.6%, D11 5.6%, D20 25.7%, D30 2.0%), so compare proportions with that in mind.
`max_percent_mt` stays at the base value of 20 for every species. Because step 01's cache key hashes
`QC`, switching species between any two with different `QC` invalidates the step 01–03 cache for a
sample. cm's downstream trajectory is `pipeline/projects/cm/trajectory.R`; see
`docs/cm_differentiation_readiness.md`.

---

## Related

- [Pipeline Steps Reference](reference-pipeline-steps.md) — what each step reads and writes
- [How to Override Annotations](howto-override-annotations.md) — step-by-step CLUSTER_CELLTYPE_MAP workflow
- [How to Run on Bat Whole Blood](howto-bat-whole-blood.md) — bat-specific parameter guidance
