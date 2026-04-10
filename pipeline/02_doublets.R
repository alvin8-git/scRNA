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
  ) + labs(title = paste0(sample_name, "  -  scDblFinder")) + theme_classic()
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

  n_doublets <- sum(seu$doublet_class == "doublet")
  message("  Doublets: ", n_doublets,
          " (", round(n_doublets / n_cells * 100, 1), "%)",
          " | Retaining: ", n_cells - n_doublets, " singlets")

  seu_singlets <- subset(seu, subset = doublet_class == "singlet")
  list(seu = seu_singlets, plots = plots)
}

# =============================================================================
# MAIN
# =============================================================================

report_plots <- list()

for (nm in SAMPLE_NAMES) {
  seu <- readRDS(file.path(DIRS$individual, nm, paste0(nm, "_filtered.rds")))
  res <- run_doubletfinder(seu, nm)
  report_plots <- c(report_plots, res$plots)
  saveRDS(res$seu, file.path(DIRS$individual, nm, paste0(nm, "_singlets.rds")))
  message("  Saved singlets: ", ncol(res$seu), " cells")
}

save_report_pdf(report_plots, file.path(DIRS$doublets, "doublets_report.pdf"))
message("\n02_doublets.R complete.")
