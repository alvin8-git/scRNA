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
  # ── HARDCODED DEFAULTS (edit here for manual runs) ─────────────────────────
  SAMPLE_PATHS <- list(
    H1 = file.path(BASE_DIR, "H1", "filter_matrix"),
    H2 = file.path(BASE_DIR, "H2", "filter_matrix")
  )
  SAMPLE_NAMES <- c("H1", "H2")
}

# Append _filtered or _raw so runs on the same sample stay in separate folders
.tags   <- sapply(unlist(SAMPLE_PATHS), .matrix_tag)
.mtag   <- if (any(.tags == "raw")) "raw" else "filtered"
rm(.env_paths, .i, .val, .FILTERED_DIRS, .RAW_DIRS, .MATRIX_DIRS,
   .resolve_sample_path, .resolve_sample_name, .matrix_tag, .tags)

SINGLE_SAMPLE <- length(SAMPLE_NAMES) == 1

# --- Results directory — named after samples for easy identification ---
# Single sample  : results_H1/
# Multiple samples: results_H1andH2/  or  results_H1andH2andH3/
RESULTS_DIR <- file.path(BASE_DIR,
  paste0("results_", paste(SAMPLE_NAMES, collapse = "and"), "_", .mtag)
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
# Core immune populations
MARKERS <- list(
  T_pan       = c("CD3D", "CD3E"),
  CD4_T       = c("CD4", "IL7R", "CCR7"),
  CD8_T       = c("CD8A", "CD8B", "GZMK"),
  Treg        = c("FOXP3", "IL2RA", "CTLA4"),        # regulatory T cells
  NK          = c("NKG7", "GNLY", "KLRD1"),
  B_cell      = c("MS4A1", "CD79A", "CD19"),
  Plasma      = c("MZB1", "JCHAIN", "SDC1"),          # plasma / plasmablasts
  CD14_mono   = c("CD14", "LYZ", "CST3", "S100A8"),
  FCGR3A_mono = c("FCGR3A", "MS4A7"),
  Neutrophil  = c("FCGR3B", "CSF3R", "CXCR2", "CEACAM8"),  # inflammation marker
  DC          = c("FCER1A", "CLEC9A"),
  Platelet    = c("PPBP", "PF4"),
  # Contamination indicators — presence signals poor sample quality
  RBC         = c("HBB", "HBA1", "HBA2", "GYPA"),    # red blood cell contamination
  HSPC        = c("CD34", "GATA2", "AVP")             # haematopoietic progenitors
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
  "Treg"         = "#FF7F0E",
  "NK"           = "#00A087",
  "B cell"       = "#3C5488",
  "Plasma"       = "#7B4F9E",
  "CD14+ Mono"   = "#F39B7F",
  "FCGR3A+ Mono" = "#8491B4",
  "Neutrophil"   = "#E377C2",   # pink — distinct from NK teal
  "DC"           = "#91D1C2",
  "Platelet"     = "#DC0000",
  # Contamination / rare populations
  "RBC"          = "#A52A2A",
  "HSPC"         = "#8C564B",   # brown — distinct from Monocyte salmon
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
