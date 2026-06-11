# =============================================================================
# config.R - Central configuration for PBMC scRNA-seq pipeline
# Source this file at the top of every pipeline script.
# =============================================================================

# --- Suppress noisy but harmless plot warnings ---
options(Seurat.warn.raster = FALSE)   # rasterizing >100k points is intentional
options(error = function() {
  msg <- tryCatch(conditionMessage(.GlobalEnv$.Last.error),
                  error = function(e) geterrmessage())
  message("\nERROR: ", msg)
  quit(status = 1, save = "no")
})
suppressWarnings(library(ggplot2))    # silence freetype/systemfonts version mismatch

# --- Shared PDF helpers (A4P, A4L, .combine_pdfs, .render, .build_page, .save_page) ---
source(file.path(dirname(sys.frame(1)$ofile %||% "."), "pdf_helpers.R"))

# --- Base paths ---
BASE_DIR          <- Sys.getenv("SCRNA_BASE_DIR", "/data/alvin/scRNA")
PIPELINE_DIR      <- file.path(BASE_DIR, "pipeline")
# Per-sample RDS cache — shared across all combo runs; delete a subfolder to force recompute
SAMPLE_CACHE_DIR  <- file.path(BASE_DIR, "sample_cache")

# --- Dynamic sample path resolution ---
# Override with environment variables (supports any number of samples):
#   SCRNA_SAMPLE1, SCRNA_SAMPLE2, SCRNA_SAMPLE3, ... = sample folder paths
# Integration (Harmony) runs automatically when >1 sample is provided.
# If no env vars set, falls back to hardcoded defaults below.

.FILTERED_DIRS <- c("filter_matrix", "filtered_feature_bc_matrix")
.RAW_DIRS      <- c("raw_matrix", "raw_feature_bc_matrix")
.MATRIX_DIRS   <- c(.FILTERED_DIRS, .RAW_DIRS)

.resolve_sample_path <- function(path) {
  matrix_files <- c("matrix.mtx.gz", "matrix.mtx", "barcodes.tsv.gz", "barcodes.tsv")
  if (any(file.exists(file.path(path, matrix_files)))) return(normalizePath(path))
  for (sub in .MATRIX_DIRS) {
    cand <- file.path(path, sub)
    if (dir.exists(cand)) return(normalizePath(cand))
  }
  normalizePath(path, mustWork = FALSE)
}

.resolve_sample_name <- function(path) {
  bn <- basename(normalizePath(path, mustWork = FALSE))
  if (bn %in% .MATRIX_DIRS)
    return(basename(dirname(normalizePath(path, mustWork = FALSE))))
  bn
}

.matrix_tag <- function(resolved_path) {
  if (basename(resolved_path) %in% .RAW_DIRS) "raw" else "filtered"
}

# Collect all SCRNA_SAMPLE{N} env vars in order
.env_paths <- character(0)
.i <- 1
repeat {
  .val <- Sys.getenv(paste0("SCRNA_SAMPLE", .i), unset = "")
  if (nchar(.val) == 0) break
  .env_paths <- c(.env_paths, .val)
  .i <- .i + 1
}

if (length(.env_paths) > 0) {
  SAMPLE_PATHS <- setNames(
    lapply(.env_paths, .resolve_sample_path),
    sapply(.env_paths, .resolve_sample_name)
  )
  SAMPLE_NAMES <- names(SAMPLE_PATHS)
} else {
  message("[config] No SCRNA_SAMPLE* env vars set — falling back to hardcoded H1/H2 defaults.")
  # ── HARDCODED DEFAULTS (edit here for manual runs) ─────────────────────────
  SAMPLE_PATHS <- list(
    H1 = file.path(BASE_DIR, "Samples", "H1", "filter_matrix"),
    H2 = file.path(BASE_DIR, "Samples", "H2", "filter_matrix")
  )
  SAMPLE_NAMES <- c("H1", "H2")
}

# Append _filtered or _raw so runs on the same sample stay in separate folders
.tags   <- sapply(unlist(SAMPLE_PATHS), .matrix_tag)
.mtag   <- if (any(.tags == "raw")) "raw" else "filtered"
rm(.env_paths, .i, .val, .FILTERED_DIRS, .RAW_DIRS, .MATRIX_DIRS,
   .resolve_sample_path, .resolve_sample_name, .matrix_tag, .tags)

SINGLE_SAMPLE <- length(SAMPLE_NAMES) == 1

# --- Condition labels (for DEG / CellChat / trajectory steps) ---
# Set via SCRNA_CONDITION env var: comma-separated name=label pairs
# Example: SCRNA_CONDITION="ADay0_healthy=healthy,BDay1_recovering=recovering"
.cond_raw <- Sys.getenv("SCRNA_CONDITION", unset = "")
if (nchar(.cond_raw) > 0) {
  .pairs <- strsplit(.cond_raw, ",")[[1]]
  .kv    <- strsplit(.pairs, "=")
  SAMPLE_CONDITIONS <- setNames(
    sapply(.kv, `[[`, 2),
    sapply(.kv, `[[`, 1)
  )
} else {
  SAMPLE_CONDITIONS <- setNames(rep("condition1", length(SAMPLE_NAMES)), SAMPLE_NAMES)
}
CONDITION_LEVELS <- unique(SAMPLE_CONDITIONS)
rm(.cond_raw)
if (exists(".pairs")) rm(.pairs)
if (exists(".kv"))    rm(.kv)

# --- Results directory — always under BASE_DIR/Results/ for organisation ---
# Single sample  : Results/results_H1_filtered/
# Multiple samples: Results/results_H1-H2-H3_filtered/
RESULTS_DIR <- file.path(BASE_DIR, "Results",
  paste0("results_", paste(SAMPLE_NAMES, collapse = "-"), "_", .mtag)
)
rm(.mtag)

# --- QC Thresholds ---
QC <- list(
  min_features   = 200,
  max_features   = 5000,
  min_counts     = 500,
  max_counts     = 25000,
  max_percent_mt = 20
)

# --- Doublet Detection ---
DOUBLET <- list(
  doublet_rate = list(H1 = 0.031, H2 = 0.077),  # NULL for unknown samples → auto
  PCs  = 1:15,
  sct  = FALSE
)

# --- Normalization & HVG ---
NORM <- list(
  method       = "LogNormalize",
  scale_factor = 10000,
  n_hvg        = 2000,
  hvg_method   = "vst"
)

# --- Dimensionality Reduction ---
DIM <- list(
  npcs      = 30,
  dims_use  = 1:20,
  umap_seed = 42
)

# --- Clustering ---
CLUSTER <- list(
  resolutions = c(0.3, 0.4, 0.5, 0.6, 0.8),
  default_res = 0.5,
  compare_res = c(0.5, 0.6, 0.8),   # resolutions shown side-by-side in comparison UMAP
  algorithm   = 1
)

# --- T Cell Sub-clustering ---
# Set enabled = FALSE to skip.  t_patterns is a regex matched against cell_type labels.
SUBCLUSTER <- list(
  enabled    = TRUE,
  t_patterns = "T[_ ]cell|T cell|CD4|CD8|Treg|cytotox",
  resolution = 0.8,    # higher res → more sub-clusters
  min_cells  = 20      # skip if fewer T cells than this
)

# --- Harmony Integration ---
HARMONY <- list(
  group_by_vars = "sample",
  theta         = 2,
  lambda        = 1,
  nclust        = 50,
  max_iter      = 20,
  dims_use      = 1:20
)

# --- Canonical PBMC Marker Genes ---
# Core immune populations
MARKERS <- list(
  T_pan       = c("CD3D", "CD3E"),
  CD4_T       = c("CD4", "IL7R", "CCR7"),
  CD8_T       = c("CD8A", "CD8B", "GZMK"),
  Treg        = c("FOXP3", "IL2RA", "CTLA4"),        # regulatory T cells
  NK          = c("NKG7", "GNLY", "KLRD1"),
  B_cell      = c("MS4A1", "CD79A", "CD19"),
  Plasma      = c("MZB1", "JCHAIN", "SDC1", "CD38", "XBP1", "PRDM1"),  # plasma / plasmablasts
  CD14_mono   = c("CD14", "LYZ", "CST3", "S100A8"),
  FCGR3A_mono = c("FCGR3A", "MS4A7"),
  Neutrophil  = c("FCGR3B", "CSF3R", "CXCR2", "CEACAM8"),  # inflammation marker
  DC          = c("FCER1A", "CLEC9A"),
  Platelet    = c("PPBP", "PF4"),
  Eosinophil  = c("SIGLEC8", "CCR3", "EPX"),           # granulocyte; rare in healthy blood
  Mast_cell   = c("TPSAB1", "CPA3", "MS4A2", "KIT"),   # tissue-resident; contamination if present
  # Contamination indicators — presence signals poor sample quality
  RBC         = c("HBB", "HBA1", "HBA2", "GYPA"),    # red blood cell contamination
  HSPC        = c("CD34", "GATA2", "AVP")             # haematopoietic progenitors
)

ALL_MARKERS <- unique(unlist(MARKERS))

# Set TRUE to run FindAllMarkers on the integrated object after Harmony.
# Saves to integrated_cluster_markers.csv. Slow: 20-40 min for large datasets.
# Not consumed by automated downstream steps; useful for manual exploration only.
MARKERS$compute_integrated <- FALSE

# --- SingleR Reference ---
# "HumanPrimaryCellAtlas" (default, broad) | "MonacoImmune" (blood-optimised, resolves CD4/CD8/γδ)
SINGLER_REF <- "HumanPrimaryCellAtlas"

# --- Sub-type Refinement Markers ---
# Used by 05_annotate.R to refine generic SingleR labels into biologically
# meaningful sub-types (e.g. "CD4 T" → "CD4 T (naive)" / "(memory)" / "(effector)").
# Each top-level key must match a label that SingleR (after SINGLER_NORM) can produce.
# Scoring: average expression of listed genes per cluster; highest score wins.
# Set to NULL to disable sub-type refinement entirely.
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

# --- Parallelism ---
# Two worker budgets, both capped at (cores - 2, max 8) and then by a RAM
# governor that leaves ~20% headroom:
#   workers        per-sample phase (steps 01-03): small objects, ~8 GB/worker.
#   merge_workers  merged-object phase (steps 04-06b): each future/BiocParallel
#                  worker can hold a copy of the full merged object, so it uses a
#                  larger per-worker budget and therefore fans out to fewer
#                  workers. Big-merge steps no longer spawn as wide as the
#                  per-sample steps (office-hours P2).
.future_mem_gb <- 8L    # per-worker budget, per-sample phase (GB)
.merge_mem_gb  <- 16L   # per-worker budget, merged-object phase (GB)
.cpu_workers   <- min(8L, max(1L, parallel::detectCores() - 2L))
.avail_gb <- tryCatch({
  kb <- as.numeric(sub("[^0-9]*([0-9]+).*", "\\1",
                       grep("MemAvailable", readLines("/proc/meminfo"), value = TRUE)))
  if (length(kb) == 1L && is.finite(kb)) kb / 1024^2 else NA_real_
}, error = function(e) NA_real_)
.gov <- function(mem_gb) if (is.finite(.avail_gb))
  max(1L, min(.cpu_workers, as.integer(floor(.avail_gb * 0.8 / mem_gb)))) else .cpu_workers
PARALLEL <- list(
  workers       = .gov(.future_mem_gb),
  merge_workers = .gov(.merge_mem_gb),
  future_mem_gb = .future_mem_gb,
  merge_mem_gb  = .merge_mem_gb
)
if (PARALLEL$workers < .cpu_workers || PARALLEL$merge_workers < PARALLEL$workers)
  message(sprintf("[config] RAM governor: %d workers per-sample / %d workers merged-phase (%.0f GB avail; %d CPU).",
                  PARALLEL$workers, PARALLEL$merge_workers, .avail_gb, .cpu_workers))
rm(.future_mem_gb, .merge_mem_gb, .cpu_workers, .avail_gb, .gov)

# --- Serialization ---
# Opt-in to qs (3-8x faster RDS I/O via LZ4/ZSTD). Falls back to base saveRDS/readRDS.
IO <- list(use_qs = requireNamespace("qs", quietly = TRUE))
.save_rds <- function(obj, path, ...) {
  if (isTRUE(IO$use_qs)) qs::qsave(obj, path, nthreads = min(4L, PARALLEL$workers), ...)
  else saveRDS(obj, path, ...)
}
.read_rds <- function(path, ...) {
  if (isTRUE(IO$use_qs)) qs::qread(path, nthreads = min(4L, PARALLEL$workers), ...)
  else readRDS(path, ...)
}

# --- Cache invalidation hash ---
# Per-step cumulative parameters. Each step's cache depends on its own params
# AND all upstream params (its input is the previous step's output), so the keys
# are nested. Splitting per step means a downstream knob (e.g. CLUSTER) no longer
# invalidates the QC (01) or doublet (02) caches.
.cache_params <- function(step) {
  species <- Sys.getenv("SCRNA_SPECIES", "human")
  p01 <- list(qc = QC, species = species)
  p02 <- c(p01, list(doublet = DOUBLET))
  p03 <- c(p02, list(norm = NORM, dim_red = DIM, cluster = CLUSTER))
  switch(step, "01" = p01, "02" = p02, "03" = p03,
         stop("cache_hash: unknown step '", step, "'"))
}

# Fingerprint the 10x matrix files (size + mtime). Lets a changed input bust the
# cache even when config is identical. Returns "" if no matrix files are found,
# so the key degrades gracefully to (params + path) for non-10x or moved inputs.
.matrix_fingerprint <- function(path) {
  if (is.null(path) || !nzchar(path)) return("")
  files <- c("matrix.mtx.gz", "matrix.mtx",
             "barcodes.tsv.gz", "barcodes.tsv",
             "features.tsv.gz", "features.tsv",
             "genes.tsv.gz", "genes.tsv")
  fp   <- file.path(path, files)
  info <- file.info(fp[file.exists(fp)])
  if (nrow(info) == 0) return("")
  paste(rownames(info), info$size, as.integer(info$mtime), collapse = ";")
}

# Per-sample, per-step cache key: cumulative step params + the sample's resolved
# input PATH + a fingerprint of its matrix files. Keying on the absolute path
# (not just the sample name `nm`) stops two different experiments that share a
# folder name from colliding in sample_cache/; the fingerprint busts the cache
# when the source matrix changes under an unchanged config.
cache_hash <- function(nm, step) {
  path <- tryCatch(SAMPLE_PATHS[[nm]], error = function(e) NULL)
  digest::digest(
    list(params      = .cache_params(step),
         path        = if (is.null(path)) nm else normalizePath(path, mustWork = FALSE),
         fingerprint = .matrix_fingerprint(path)),
    algo = "md5"
  )
}

# --- Manual Cluster → Cell Type Map ---
# NULL  → auto-annotate using SingleR majority vote per cluster (recommended for first run).
#         The log (logs/05_annotate.log) will print a copy-pasteable map you can paste below.
# c(...) → use the map below; any cluster NOT listed falls back to SingleR automatically.
#
# NOTE: cluster numbers change between datasets — do not copy a map from one sample to another.
#
# Verified map for DemoScRNA_filtered (13 clusters, res=0.5):
#   "0"  = "CD4 T (naive)"    CCR7, TCF7, LEF1, IL7R
#   "1"  = "CD4 T (effector)" CD3D/E, GZMK, TNFRSF4 (OX40)
#   "2"  = "CD4 T (memory)"   CD3D/E, lower CCR7/LEF1
#   "3"  = "NK"               GNLY, NKG7, KLRD1, GZMB
#   "4"  = "B cell (naive)"   TCL1A, IGHD, IGHM, CD79A
#   "5"  = "FCGR3A+ Mono"     CDKN1C, FCGR3A, MS4A7
#   "6"  = "CD14+ Mono"       S100A8, LYZ, CD14, FCAR
#   "7"  = "Neutrophil"       S100A12, S100A9, BST1, G0S2
#   "8"  = "B cell (memory)"  IGHG1/G2, IGHA1, TNFRSF13B
#   "9"  = "γδ T"             TRGC2, TRGC1
#   "10" = "CD8 T"            CD8A, CD8B, CCR7, LEF1
#   "11" = "DC"               FCER1A, CLEC10A, FLT3
#   "12" = "Platelet"         PPBP, PF4, ITGB3, GP9
CLUSTER_CELLTYPE_MAP <- NULL  # reset for each new analysis

# --- Color Palettes ---
SAMPLE_COLORS <- c(H1 = "#E64B35", H2 = "#4DBBD5")

# Cell types always preserved at per-cell level — never overridden by cluster majority vote.
# These are contamination or rare types expected in low-quality or sorted samples.
# Add / remove types to control which populations bypass the majority-vote labelling.
CONTAMINATION_TYPES <- c("Neutrophil", "RBC", "HSPC", "Platelet",
                          "Basophil", "Eosinophil", "Mast cell")

CELLTYPE_COLORS <- c(
  # CD4 T subtypes — red family
  "CD4 T"              = "#E64B35",
  "CD4 T (naive)"      = "#E64B35",
  "CD4 T (memory)"     = "#FF7043",
  "CD4 T (effector)"   = "#FF8A65",
  # CD8 T subtypes — blue family
  "CD8 T"              = "#4DBBD5",
  "CD8 T (naive)"      = "#6ACDE6",
  "CD8 T (memory)"     = "#3A8BA5",
  "CD8 T (effector)"   = "#1D6680",
  # Other T
  "Treg"               = "#FF7F0E",
  "γδ T"               = "#FFC107",
  "NKT"                = "#17BECF",
  # NK — teal
  "NK"                 = "#00A087",
  # B cell subtypes — dark blue family
  "B cell"             = "#3C5488",
  "B cell (naive)"     = "#3C5488",
  "B cell (memory)"    = "#5C74A8",
  "Plasma"             = "#7B4F9E",
  # Monocyte subtypes — salmon/purple family
  "Monocyte"           = "#F39B7F",
  "CD14+ Mono"         = "#F39B7F",
  "FCGR3A+ Mono"       = "#8491B4",
  # Myeloid
  "Neutrophil"         = "#E377C2",   # pink — distinct from NK teal
  "DC"                 = "#91D1C2",
  "cDC1"               = "#70BFB0",
  "cDC2"               = "#A8E6D8",
  "pDC"                = "#4FA090",
  "Platelet"           = "#DC0000",
  # Contamination / rare populations
  "RBC"                = "#A52A2A",
  "HSPC"               = "#8C564B",   # brown — distinct from Monocyte salmon
  "Basophil"           = "#B5B000",
  "Eosinophil"         = "#F4A460",
  "Mast cell"          = "#9400D3",
  # Non-immune / stromal (tissue contamination)
  "Endothelial"        = "#636363",
  "Epithelial"         = "#969696",
  "Fibroblast"         = "#BDBDBD",
  "Smooth Muscle"      = "#D9D9D9",
  "Unknown"            = "#B09C85"
)


# =============================================================================
# Species overrides — read SCRNA_SPECIES env var (set by run_pipeline.sh)
# Supported: "human" (default) | "bat" (whole blood) | "bat_wing" (wing tissue)
# =============================================================================
.species <- Sys.getenv("SCRNA_SPECIES", unset = "")
if (!nzchar(.species)) {
  message("[config] SCRNA_SPECIES not set — defaulting to 'human' (set SCRNA_SPECIES=bat for bat data).")
  .species <- "human"
}

source(file.path(PIPELINE_DIR, "config_species_bat.R"))  # bat / bat_wing overrides (extracted)

if (!exists("WOUND_MODULES")) WOUND_MODULES <- list()

rm(.species)

# --- Plot Defaults ---
PLOT <- list(
  width      = 8,
  height     = 7,
  dpi        = 300,
  pt_size    = 0.8,
  label_size = 4
)

# --- Utility ---
`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0) x else y

# --- Output Subdirectories ---
DIRS <- list(
  qc           = file.path(RESULTS_DIR, "qc"),
  doublets     = file.path(RESULTS_DIR, "doublets"),
  individual   = file.path(RESULTS_DIR, "individual"),
  integrated   = file.path(RESULTS_DIR, "integrated"),
  annotation   = file.path(RESULTS_DIR, "annotation"),
  differential = file.path(RESULTS_DIR, "differential"),
  pathways     = file.path(RESULTS_DIR, "pathways"),
  cellchat     = file.path(RESULTS_DIR, "cellchat"),
  trajectory   = file.path(RESULTS_DIR, "trajectory"),
  logs         = file.path(RESULTS_DIR, "logs"),
  reports      = RESULTS_DIR
)

invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))
invisible(lapply(file.path(DIRS$individual, SAMPLE_NAMES),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

# =============================================================================
# PDF Report Helpers
# =============================================================================

# Tag a plot with page dimensions and whether it's "small" (2-per-page eligible)
set_page <- function(p, pw = 8.5, ph = 7.5, small = FALSE) {
  attr(p, "pw")    <- pw
  attr(p, "ph")    <- ph
  attr(p, "small") <- small
  p
}
mark_small <- function(p) set_page(p, pw = 8.0, ph = 5.0, small = TRUE)

# --- Per-plot caption lookup (pattern → "description\nGood: ... | Bad: ...") ---
PLOT_CAPTIONS <- list(
  "QC Violin" =
    "Gene count, UMI count, and mitochondrial % distributions per cell.\nGood: Unimodal peaks; %MT median < 10%.  Bad: Bimodal nFeature or high %MT peak (dead/lysed cells).",
  "QC Scatter" =
    "UMI vs genes (left) and UMI vs %MT (right). Red dashed lines = filter thresholds.\nGood: Tight diagonal band; cells inside threshold box.  Bad: Cloud above %MT line or below gene line (empty drops / dying cells).",
  "Doublet.*UMAP" =
    "UMAP coloured by doublet (red) vs singlet (grey) classification from scDblFinder.\nGood: Doublets at cluster boundaries (mixed-identity droplets).  Bad: Doublets spread uniformly = miscalibrated rate or poor library.",
  "Doublet.*Hist|Score Dist" =
    "Doublet probability score distribution (0 = singlet, 1 = doublet).\nGood: Bimodal — large peak near 0, small peak near 1.  Bad: Single broad peak = doublets indistinguishable from singlets.",
  "Variable Genes|Highly Variable" =
    "All genes by average expression vs variance; top 2000 HVGs (red) drive PCA.\nGood: Smooth curve; known markers (CD3D, CD14) in top HVGs.  Bad: Only ribo/MT genes in top 10 = normalisation issue.",
  "Elbow" =
    "Variance explained per PC; red dashed line = PCs used for UMAP and clustering.\nGood: Clear elbow with cutoff just past it (PC 10-20 for PBMC).  Bad: No elbow = no structure; cutoff before elbow = noisy clusters.",
  "PC Heatmap" =
    "Top/bottom loading genes per PC across cells ordered by PC score.\nGood: Distinct gradient with known markers per PC.  Bad: Ribo (RPL/RPS) or stress genes dominate = technical, not biological, PCs.",
  "UMAP.*Cluster|Cluster.*UMAP" =
    "UMAP embedding coloured by Seurat cluster identity.\nGood: Compact, well-separated clusters.  Bad: One blob = too few PCs; many tiny clusters = resolution too high.",
  "Top.*Marker.*Dot|Dot.*Top.*Marker" =
    "Top 5 DE genes per cluster: dot size = % expressing, colour = avg expression.\nGood: Each cluster has unique markers.  Bad: Shared markers or low fold-change = unresolved populations.",
  "Canonical Marker.*Feature|Feature.*Canonical" =
    "Canonical PBMC marker expression overlaid on UMAP (grey = low, red = high).\nGood: Each marker enriched in one UMAP region.  Bad: Diffuse uniform expression = normalisation failure or very noisy data.",
  "SingleR.*Score Heatmap|Score Heatmap" =
    "SingleR reference scores per cell (columns) vs reference cell types (rows).\nGood: One bright score per cell, rest dark = confident annotation.  Bad: Uniformly moderate scores = ambiguous cell type.",
  "SingleR.*Label|Labels.*UMAP" =
    "UMAP coloured by pruned SingleR automated cell type labels.\nGood: Labels spatially coherent and match expected tissue types.  Bad: Scattered labels = wrong reference; many 'Unassigned' = reference gap.",
  "SingleR.*Delta|Delta Score" =
    "SingleR annotation confidence (top score minus second-best) on UMAP.\nGood: High delta (dark blue) in cluster centres.  Bad: All low delta (<0.1) = reference cannot discriminate cell types.",
  "Canonical PBMC Markers Dot" =
    "Canonical PBMC markers vs clusters — use this to fill CLUSTER_CELLTYPE_MAP in config.R.\nGood: Each cluster expresses one lineage's markers clearly.  Bad: Mixed or absent signal = unresolved clusters or wrong tissue.",
  "Feature.*T Cell|T Cell.*Marker" =
    "T cell marker expression (CD3D, CD4, CD8A, CCR7, GZMK) on UMAP.\nGood: CD3D marks all T cells; CD4/CD8A localise to separate non-overlapping regions.  Bad: Full overlap = CD4/CD8 not resolved.",
  "Feature.*NK|NK.*Marker" =
    "NK cell markers (NKG7, GNLY, KLRD1) on UMAP.\nGood: Co-localised region distinct from CD3D+ T cells.  Bad: Same region as T cells = NK/cytotoxic T not separated.",
  "Feature.*B Cell|B Cell.*Marker" =
    "B cell markers (MS4A1, CD79A, CD19) on UMAP.\nGood: Tight co-expression in an isolated cluster.  Bad: MS4A1 in monocyte region = likely annotation error.",
  "Feature.*Mono|Mono.*Marker" =
    "Monocyte markers (CD14, LYZ, FCGR3A, MS4A7) on UMAP.\nGood: CD14 and FCGR3A mark adjacent but distinct sub-clusters.  Bad: Full overlap = CD14/FCGR3A monocytes unresolved.",
  "Feature.*DC|DC.*Platelet" =
    "DC (FCER1A, CLEC9A) and platelet (PPBP, PF4) markers — typically rare populations.\nGood: Small but distinct expression clusters.  Bad: Absent = rare types filtered out or too few cells to cluster.",
  "Cell Type.*UMAP|UMAP.*Cell Type" =
    "UMAP coloured by final annotated cell type labels.\nGood: Spatially coherent, non-overlapping regions per type.  Bad: Intermixed types = annotation errors or insufficient cluster resolution.",
  "Harmony|Before.*After|After.*Harmony" =
    "UMAP by sample before (left) and after (right) Harmony batch correction.\nGood: Samples separate before; fully intermixed by cell type after.  Bad: Still separated after = increase HARMONY\\$theta (try 3-5).",
  "Split.*Sample|UMAP.*Split" =
    "Integrated UMAP split per sample — same embedding, one panel each.\nGood: Similar cluster layout and proportions across samples.  Bad: Missing clusters in one sample = absent cell type or sample QC failure.",
  "Canonical.*Cell Type|Dot.*Cell Type" =
    "Canonical markers vs annotated cell types: size = % expressing, colour = avg expression.\nGood: Diagonal specificity — each type expresses its expected markers only.  Bad: Cross-lineage expression = annotation needs refinement.",
  "Heatmap.*Marker|Top.*Markers.*Cluster" =
    "Top 3 DE genes per cluster across downsampled cells, grouped by annotated cell type.\nGood: Clear colour blocks per cell type.  Bad: Markers bleed across types = similar populations or incorrect annotation.",
  "Composition|Proportion|% of Sample|Number of Cells" =
    "Cell type proportions (%) and absolute counts per sample.\nGood: Consistent proportions across samples for stable populations.  Bad: Large differences = biological variation or uncorrected batch effect.",
  "Violin.*Marker|Key.*Lineage|Key Marker" =
    "Expression of 8 key lineage marker genes across annotated cell types.\nGood: Each marker high in its expected type and near-zero in all others.  Bad: Broad expression = cell type labels need revision."
)

.get_caption <- function(nm) {
  if (is.null(nm) || is.na(nm) || nchar(nm) == 0) return(NULL)
  for (pat in names(PLOT_CAPTIONS))
    if (grepl(pat, nm, ignore.case = TRUE, perl = TRUE)) return(PLOT_CAPTIONS[[pat]])
  NULL
}

# Save a named list of ggplot/patchwork objects to a multi-page PDF.
# Items may also be character file paths to existing PDFs (e.g. base-R plots).
# Small plots (attr "small" == TRUE) are paired 2-per-portrait-page.
# Names become bold title banners above each figure.
save_report_pdf <- function(plots, filepath) {
  suppressPackageStartupMessages({ library(cowplot); library(ggplot2) })

  wrap_gg <- function(p, nm) {
    # 1. Title banner
    if (!is.null(nm) && !is.na(nm) && nchar(nm) > 0) {
      if (inherits(p, "patchwork")) {
        p <- p + patchwork::plot_annotation(
          title = nm,
          theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0))
        )
      } else {
        title_row <- ggdraw() +
          draw_label(nm, fontface = "bold", size = 13, x = 0.04, y = 0.5, hjust = 0)
        p <- plot_grid(title_row, p, ncol = 1, rel_heights = c(0.06, 0.94))
      }
    }
    # 2. Caption row (description + Good/Bad guide)
    cap <- .get_caption(nm)
    if (!is.null(cap)) {
      cap_row <- ggdraw() +
        draw_label(cap, size = 7.5, color = "#444444",
                   x = 0.02, y = 0.55, hjust = 0, vjust = 1, lineheight = 1.3)
      p <- plot_grid(p, cap_row, ncol = 1, rel_heights = c(1, 0.14))
    }
    p
  }

  `%||%` <- function(x, y) if (!is.null(x) && length(x) > 0) x else y

  tmp_files <- character(0)
  i <- 1; n <- length(plots)

  while (i <= n) {
    p  <- plots[[i]]
    nm <- names(plots)[i]

    # Existing PDF file path — include directly
    if (is.character(p)) {
      if (file.exists(p)) tmp_files <- c(tmp_files, p)
      i <- i + 1; next
    }

    pw <- attr(p, "pw") %||% 8.5
    ph <- attr(p, "ph") %||% 7.5
    sm <- isTRUE(attr(p, "small"))

    # Try to pair two consecutive small plots on one portrait page
    next_is_small <- (i + 1) <= n &&
      !is.character(plots[[i + 1]]) &&
      isTRUE(attr(plots[[i + 1]], "small"))

    tf <- tempfile(fileext = ".pdf")
    if (sm && next_is_small) {
      p2 <- plots[[i + 1]]; nm2 <- names(plots)[i + 1]
      pdf(tf, width = 8.5, height = 11)
      print(plot_grid(wrap_gg(p, nm), wrap_gg(p2, nm2), ncol = 1, nrow = 2))
      dev.off()
      tmp_files <- c(tmp_files, tf)
      i <- i + 2
    } else {
      pdf(tf, width = pw, height = ph)
      print(wrap_gg(p, nm))
      dev.off()
      tmp_files <- c(tmp_files, tf)
      i <- i + 1
    }
  }

  .combine_pdfs(tmp_files, filepath)
  # Remove only temp files (not user-supplied paths)
  unlink(tmp_files[startsWith(tmp_files, tempdir())])
  message("  Report saved: ", filepath)
  invisible(filepath)
}

# =============================================================================
# ARCHIVED CLUSTER_CELLTYPE_MAPs — reference only, do NOT activate.
# Each entry includes the dataset name and run date.
# Current active map is set in the Manual Cluster -> Cell Type Map section above.
# =============================================================================

# Verified map for ES03_newkit + ES12_newkit (bat, 17 clusters, res=1.0, MonacoImmune)
# NK clusters 4/7/11/13 confirmed CD3+ cytotoxic T cells by marker check (CD3D 50-70%, CD3E 85-89%)
# Cluster 16 confirmed Platelet by PF4 89%, ITGA2B 72%, GP9 50%
# All other clusters fall back to SingleR (Monaco) automatically
# NOTE: set to NULL for fresh multi-sample runs — cluster numbers change between datasets
# 4-sample run (Sample6, Sample7, 10, ES03_newkit): clusters 3/11/15 confirmed CD8 T (CTL)
# by CD3E 90.7%, NCAM1 0.5% — SingleR incorrectly calls these NK
# ARCHIVED — do not reuse; cluster numbers change between datasets:
#   c("3"="CD8 T", "11"="CD8 T", "15"="CD8 T", "13"="DC", "18"="DC", "17"="Neutrophil")
