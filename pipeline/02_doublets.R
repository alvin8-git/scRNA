# =============================================================================
# 02_doublets.R - Doublet detection using scDblFinder (Bioconductor).
#                 Outputs: per-sample *_singlets.rds + doublets_report.pdf
# =============================================================================
source("/data/alvin/scRNA/pipeline/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(BiocParallel)
  library(ggplot2)
  library(dplyr)
})

bp_param <- MulticoreParam(workers = PARALLEL$workers)
message("Parallelism: ", PARALLEL$workers, " cores (BiocParallel)")

# =============================================================================
# FUNCTION  -  returns list(seu = singlets_seurat, plots = named_list)
# =============================================================================

run_doubletfinder <- function(seu, sample_name) {
  message("\n=== Doublet detection: ", sample_name,
          " (", ncol(seu), " cells) ===")
  n_cells <- ncol(seu)
  plots   <- list()

  sce <- as.SingleCellExperiment(seu)
  set.seed(DIM$umap_seed)
  sce <- scDblFinder(
    sce,
    dbr       = DOUBLET$doublet_rate[[sample_name]],  # NULL → auto
    nfeatures = 1000,
    BPPARAM   = bp_param,
    verbose   = FALSE
  )

  seu$doublet_score <- sce$scDblFinder.score
  seu$doublet_class <- sce$scDblFinder.class

  # Compute doublet stats before plotting so we can embed them in the title
  n_doublets <- sum(seu$doublet_class == "doublet")
  pct_dbl    <- round(n_doublets / n_cells * 100, 1)

  # Quick UMAP for visualisation only
  seu_vis <- NormalizeData(seu, verbose = FALSE)
  seu_vis <- FindVariableFeatures(seu_vis, nfeatures = 2000, verbose = FALSE)
  seu_vis <- ScaleData(seu_vis, verbose = FALSE)
  seu_vis <- RunPCA(seu_vis, npcs = 20, verbose = FALSE)
  seu_vis <- RunUMAP(seu_vis, dims = 1:15, seed.use = DIM$umap_seed,
                     umap.method = "uwot", verbose = FALSE)

  p_umap <- DimPlot(seu_vis,
    group.by = "doublet_class",
    cols     = c("singlet" = "#CCCCCC", "doublet" = "#E64B35"),
    pt.size  = 1.0
  ) + labs(title = paste0(sample_name, "  -  Doublet UMAP  (",
                          n_doublets, " in ", n_cells, " cells  (", pct_dbl, "%))")) +
    theme_classic()
  ggsave(file.path(DIRS$doublets, paste0(sample_name, "_doublet_umap.pdf")),
         p_umap, width = 7, height = 6)
  plots[[paste0(sample_name, "  -  Doublet UMAP")]] <- mark_small(p_umap)

  p_hist <- ggplot(seu@meta.data, aes(x = doublet_score, fill = doublet_class)) +
    geom_histogram(bins = 40, alpha = 0.8, position = "identity") +
    scale_fill_manual(values = c("singlet" = "#CCCCCC", "doublet" = "#E64B35")) +
    labs(title = paste0(sample_name, "  -  Doublet Score Distribution"),
         x = "scDblFinder Score", y = "Count") +
    theme_classic()
  ggsave(file.path(DIRS$doublets, paste0(sample_name, "_doublet_score_hist.pdf")),
         p_hist, width = 6, height = 4)
  plots[[paste0(sample_name, "  -  Doublet Score Histogram")]] <- mark_small(p_hist)

  message("  Doublets: ", n_doublets,
          " (", pct_dbl, "%)",
          " | Retaining: ", n_cells - n_doublets, " singlets")

  seu_singlets <- subset(seu, subset = doublet_class == "singlet")
  list(seu = seu_singlets, plots = plots)
}

# =============================================================================
# MAIN
# =============================================================================

report_plots <- list()

for (nm in SAMPLE_NAMES) {
  out_path   <- file.path(DIRS$individual, nm, paste0(nm, "_singlets.rds"))
  cache_path <- file.path(SAMPLE_CACHE_DIR, nm, paste0(nm, "_singlets.rds"))
  hash_path  <- sub("\\.rds$", ".hash", cache_path)
  if (file.exists(cache_path)) {
    stored <- if (file.exists(hash_path)) readLines(hash_path) else ""
    if (stored == cache_hash()) {
      message("  [CACHE HIT] ", nm, ": loading _singlets.rds from sample_cache/ (skipping doublet detection)")
      file.copy(cache_path, out_path, overwrite = TRUE)
    } else {
      message("  [CACHE STALE] ", nm, ": config changed — recomputing")
      seu <- readRDS(file.path(DIRS$individual, nm, paste0(nm, "_filtered.rds")))
      res <- run_doubletfinder(seu, nm)
      report_plots <- c(report_plots, res$plots)
      saveRDS(res$seu, out_path)
      message("  Saved singlets: ", ncol(res$seu), " cells")
      dir.create(file.path(SAMPLE_CACHE_DIR, nm), recursive = TRUE, showWarnings = FALSE)
      file.copy(out_path, cache_path, overwrite = TRUE)
      writeLines(cache_hash(), hash_path)
      message("  [CACHE] Saved ", nm, "_singlets.rds to sample_cache/")
    }
  } else {
    seu <- readRDS(file.path(DIRS$individual, nm, paste0(nm, "_filtered.rds")))
    res <- run_doubletfinder(seu, nm)
    report_plots <- c(report_plots, res$plots)
    saveRDS(res$seu, out_path)
    message("  Saved singlets: ", ncol(res$seu), " cells")
    dir.create(file.path(SAMPLE_CACHE_DIR, nm), recursive = TRUE, showWarnings = FALSE)
    file.copy(out_path, cache_path, overwrite = TRUE)
    writeLines(cache_hash(), hash_path)
    message("  [CACHE] Saved ", nm, "_singlets.rds to sample_cache/")
  }
}

save_report_pdf(report_plots, file.path(DIRS$doublets, "doublets_report.pdf"))

# --- Update cell_fate.csv with doublet removal counts ---
fate_file <- file.path(DIRS$qc, "cell_fate.csv")
if (file.exists(fate_file)) {
  cell_fate <- read.csv(fate_file, stringsAsFactors = FALSE)
  for (nm in SAMPLE_NAMES) {
    seu_post <- readRDS(file.path(DIRS$individual, nm, paste0(nm, "_singlets.rds")))
    idx <- match(nm, cell_fate$Sample)
    if (!is.na(idx)) {
      n_before_dbl   <- cell_fate$After_QC[idx]
      n_after_dbl    <- ncol(seu_post)
      cell_fate$After_doublet_removal[idx] <- n_after_dbl
      cell_fate$Removed_doublets[idx]      <- n_before_dbl - n_after_dbl
      cell_fate$Total_removed[idx]         <- cell_fate$GEM_barcodes[idx] - n_after_dbl
      cell_fate$Pct_retained[idx]          <- round(n_after_dbl /
                                               cell_fate$GEM_barcodes[idx] * 100, 1)
    }
  }
  write.csv(cell_fate, fate_file, row.names = FALSE)
  message("\n--- Cell Fate (load + QC + doublets) ---")
  print(cell_fate[, c("Sample", "GEM_barcodes", "Removed_load_low_genes",
                      "Removed_QC", "Removed_doublets",
                      "After_doublet_removal", "Pct_retained")])
}

message("\n02_doublets.R complete.")
