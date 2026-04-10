# =============================================================================
# config.R - Central configuration for PBMC scRNA-seq pipeline
# Source this file at the top of every pipeline script.
# =============================================================================

# --- Base paths ---
BASE_DIR     <- "/data/alvin/scRNA"
PIPELINE_DIR <- file.path(BASE_DIR, "pipeline")

# --- Dynamic sample path resolution ---
# Override with environment variables (supports any number of samples):
#   SCRNA_SAMPLE1, SCRNA_SAMPLE2, SCRNA_SAMPLE3, ... = sample folder paths
# Integration (Harmony) runs automatically when >1 sample is provided.
# If no env vars set, falls back to hardcoded defaults below.

.resolve_sample_path <- function(path) {
  matrix_files <- c("matrix.mtx.gz", "matrix.mtx", "barcodes.tsv.gz", "barcodes.tsv")
  if (any(file.exists(file.path(path, matrix_files)))) return(normalizePath(path))
  # Check known subfolder names in priority order
  for (sub in c("filter_matrix", "filtered_feature_bc_matrix",
                "raw_matrix", "raw_feature_bc_matrix")) {
    cand <- file.path(path, sub)
    if (dir.exists(cand)) return(normalizePath(cand))
  }
  normalizePath(path, mustWork = FALSE)
}

.resolve_sample_name <- function(path) {
  bn <- basename(normalizePath(path, mustWork = FALSE))
  if (bn %in% c("filter_matrix", "filtered_feature_bc_matrix",
                "raw_matrix", "raw_feature_bc_matrix"))
    return(basename(dirname(normalizePath(path, mustWork = FALSE))))
  bn
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
  # ── HARDCODED DEFAULTS (edit here for manual runs) ─────────────────────────
  SAMPLE_PATHS <- list(
    H1 = file.path(BASE_DIR, "H1", "filter_matrix"),
    H2 = file.path(BASE_DIR, "H2", "filter_matrix")
  )
  SAMPLE_NAMES <- c("H1", "H2")
}
rm(.env_paths, .i, .val, .resolve_sample_path, .resolve_sample_name)

SINGLE_SAMPLE <- length(SAMPLE_NAMES) == 1

# --- Results directory — named after samples for easy identification ---
# Single sample  : results_H1/
# Multiple samples: results_H1andH2/  or  results_H1andH2andH3/
RESULTS_DIR <- file.path(BASE_DIR,
  paste0("results_", paste(SAMPLE_NAMES, collapse = "and"))
)

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
  algorithm   = 1
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

ALL_MARKERS <- unique(unlist(MARKERS))

# --- Parallelism ---
# Workers: all cores minus 2 (keep system responsive), capped at 8
PARALLEL <- list(
  workers          = min(8L, max(1L, parallel::detectCores() - 2L)),
  # Memory per future worker — reduce if OOM errors occur
  future_mem_gb    = 4L
)

# --- Manual Cluster → Cell Type Map ---
# Fill after inspecting 05_annotate.R outputs.
# Set NULL to use SingleR majority labels as fallback.
CLUSTER_CELLTYPE_MAP <- NULL

# --- Color Palettes ---
SAMPLE_COLORS <- c(H1 = "#E64B35", H2 = "#4DBBD5")

CELLTYPE_COLORS <- c(
  "CD4 T"        = "#E64B35",
  "CD8 T"        = "#4DBBD5",
  "NK"           = "#00A087",
  "B cell"       = "#3C5488",
  "CD14+ Mono"   = "#F39B7F",
  "FCGR3A+ Mono" = "#8491B4",
  "DC"           = "#91D1C2",
  "Platelet"     = "#DC0000",
  "Unknown"      = "#B09C85"
)

# --- Plot Defaults ---
PLOT <- list(
  width      = 8,
  height     = 7,
  dpi        = 300,
  pt_size    = 0.8,
  label_size = 4
)

# --- Output Subdirectories ---
DIRS <- list(
  qc         = file.path(RESULTS_DIR, "qc"),
  doublets   = file.path(RESULTS_DIR, "doublets"),
  individual = file.path(RESULTS_DIR, "individual"),
  integrated = file.path(RESULTS_DIR, "integrated"),
  annotation = file.path(RESULTS_DIR, "annotation"),
  logs       = file.path(RESULTS_DIR, "logs"),
  reports    = RESULTS_DIR
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

# Combine a vector of PDF file paths into one PDF
.combine_pdfs <- function(files, output) {
  files <- files[file.exists(files)]
  if (length(files) == 0) { message("  No PDFs to combine for: ", output); return(invisible(NULL)) }
  if (length(files) == 1) { file.copy(files, output, overwrite = TRUE); return(invisible(output)) }
  if (requireNamespace("qpdf", quietly = TRUE)) {
    qpdf::pdf_combine(files, output)
  } else if (nchar(Sys.which("pdfunite")) > 0) {
    system2("pdfunite", c(shQuote(files), shQuote(output)))
  } else if (nchar(Sys.which("gs")) > 0) {
    system2("gs", c("-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=pdfwrite",
                    paste0("-sOutputFile=", shQuote(output)), shQuote(files)))
  } else {
    warning("Cannot combine PDFs — install R package 'qpdf' or system 'pdfunite'. Copying first file.")
    file.copy(files[1], output, overwrite = TRUE)
  }
  invisible(output)
}

# Save a named list of ggplot/patchwork objects to a multi-page PDF.
# Items may also be character file paths to existing PDFs (e.g. base-R plots).
# Small plots (attr "small" == TRUE) are paired 2-per-portrait-page.
# Names become bold title banners above each figure.
save_report_pdf <- function(plots, filepath) {
  suppressPackageStartupMessages({ library(cowplot); library(ggplot2) })

  wrap_gg <- function(p, nm) {
    if (is.null(nm) || is.na(nm) || nchar(nm) == 0) return(p)
    if (inherits(p, "patchwork")) {
      return(p + patchwork::plot_annotation(
        title = nm,
        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0))
      ))
    }
    title_row <- ggdraw() +
      draw_label(nm, fontface = "bold", size = 13, x = 0.04, y = 0.5, hjust = 0)
    plot_grid(title_row, p, ncol = 1, rel_heights = c(0.06, 0.94))
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
