# =============================================================================
# 03_individual.R - Per-sample: Normalize → HVG → Scale → PCA → UMAP →
#                   Cluster → Markers
#                   Outputs: per-sample *_seurat.rds + individual_report.pdf
# =============================================================================
source("/data/alvin/scRNA/pipeline/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(future.apply)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

plan("multicore", workers = PARALLEL$workers)
options(future.globals.maxSize = PARALLEL$future_mem_gb * 1024^3)
message("Parallelism: ", PARALLEL$workers, " cores (future multicore)")

# =============================================================================
# FUNCTION  -  returns list(seu = processed_seurat, plots = named_list)
# =============================================================================

process_individual <- function(seu, sample_name) {
  message("\n=== Individual analysis: ", sample_name,
          " (", ncol(seu), " cells) ===")
  outdir <- file.path(DIRS$individual, sample_name)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  plots <- list()

  # --- Normalize ---
  seu <- NormalizeData(seu, normalization.method = NORM$method,
                       scale.factor = NORM$scale_factor, verbose = FALSE)

  # --- HVG ---
  seu <- FindVariableFeatures(seu, selection.method = NORM$hvg_method,
                               nfeatures = NORM$n_hvg, verbose = FALSE)
  top10   <- head(VariableFeatures(seu), 10)
  p_hvg   <- VariableFeaturePlot(seu)
  p_hvg   <- LabelPoints(p_hvg, points = top10, repel = TRUE, xnudge = 0)
  p_hvg   <- p_hvg + labs(title = paste0(sample_name, "  -  Variable Genes"))
  ggsave(file.path(outdir, paste0(sample_name, "_hvg.pdf")), p_hvg, width = 8, height = 5)
  plots[[paste0(sample_name, "  -  Highly Variable Genes")]] <- mark_small(p_hvg)

  # --- Scale + PCA ---
  seu <- ScaleData(seu, vars.to.regress = "percent.mt", verbose = FALSE)
  seu <- RunPCA(seu, npcs = DIM$npcs, verbose = FALSE)

  p_elbow <- ElbowPlot(seu, ndims = DIM$npcs) +
    geom_vline(xintercept = max(DIM$dims_use), linetype = "dashed", color = "red") +
    labs(title = paste0(sample_name, "  -  Elbow Plot")) + theme_classic()
  ggsave(file.path(outdir, paste0(sample_name, "_elbow.pdf")), p_elbow, width = 6, height = 4)
  plots[[paste0(sample_name, "  -  PCA Elbow Plot")]] <- mark_small(p_elbow)

  # PC heatmaps (base R  -  save to file, include path for PDF combiner)
  pc_heatmap_path <- file.path(outdir, paste0(sample_name, "_pc_heatmaps.pdf"))
  pdf(pc_heatmap_path, width = 12, height = 15)
  DimHeatmap(seu, dims = 1:9, cells = min(300, ncol(seu)), balanced = TRUE)
  dev.off()
  plots[[paste0(sample_name, "  -  PC Heatmaps (dims 1-9)")]] <- pc_heatmap_path

  # --- UMAP + Clustering ---
  seu <- RunUMAP(seu, dims = DIM$dims_use, seed.use = DIM$umap_seed,
                 umap.method = "uwot", verbose = FALSE)
  seu <- FindNeighbors(seu, dims = DIM$dims_use, verbose = FALSE)

  for (res in CLUSTER$resolutions) {
    seu <- FindClusters(seu, resolution = res, algorithm = CLUSTER$algorithm, verbose = FALSE)
    message("  res=", res, ": ", length(unique(seu$seurat_clusters)), " clusters")
  }
  default_col <- paste0("RNA_snn_res.", CLUSTER$default_res)
  Idents(seu) <- default_col
  seu$seurat_clusters <- seu@meta.data[[default_col]]

  p_cluster <- DimPlot(seu, reduction = "umap", label = TRUE,
                        label.size = PLOT$label_size, pt.size = PLOT$pt_size) +
    labs(title = paste0(sample_name, "  -  Clusters (res=", CLUSTER$default_res, ")")) +
    theme_classic()
  ggsave(file.path(outdir, paste0(sample_name, "_umap_cluster.pdf")),
         p_cluster, width = PLOT$width, height = PLOT$height)
  plots[[paste0(sample_name, "  -  UMAP Clusters")]] <-
    set_page(p_cluster, pw = 8.5, ph = 7.5)

  p_sample <- DimPlot(seu, group.by = "sample",
                       cols = SAMPLE_COLORS[sample_name], pt.size = PLOT$pt_size) +
    labs(title = paste0(sample_name, "  -  Sample")) + theme_classic()
  ggsave(file.path(outdir, paste0(sample_name, "_umap_sample.pdf")),
         p_sample, width = PLOT$width, height = PLOT$height)
  plots[[paste0(sample_name, "  -  UMAP by Sample")]] <- mark_small(p_sample)

  # --- Cluster markers ---
  message("  Finding cluster markers...")
  all_markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1,
                                 logfc.threshold = 0.25, test.use = "wilcox",
                                 verbose = FALSE)
  write.csv(all_markers,
            file.path(outdir, paste0(sample_name, "_cluster_markers.csv")),
            row.names = FALSE)

  top5 <- all_markers %>%
    group_by(cluster) %>% slice_max(avg_log2FC, n = 5) %>%
    pull(gene) %>% unique()

  p_dot <- DotPlot(seu, features = top5) + RotatedAxis() +
    labs(title = paste0(sample_name, "  -  Top 5 Markers per Cluster"))
  ggsave(file.path(outdir, paste0(sample_name, "_dotplot_markers.pdf")),
         p_dot, width = max(12, length(top5) * 0.4), height = 7)
  plots[[paste0(sample_name, "  -  Top Marker Dot Plot")]] <-
    set_page(p_dot, pw = 11, ph = 7)

  # --- Canonical marker feature plots ---
  markers_present <- ALL_MARKERS[ALL_MARKERS %in% rownames(seu)]
  p_feat <- FeaturePlot(seu, features = markers_present, ncol = 4,
                         pt.size = 0.5, order = TRUE,
                         cols = c("lightgrey", "#E64B35")) &
    theme(plot.title = element_text(size = 9),
          axis.text = element_blank(), axis.ticks = element_blank())
  feat_h <- ceiling(length(markers_present) / 4) * 4
  ggsave(file.path(outdir, paste0(sample_name, "_umap_markers.pdf")),
         p_feat, width = 16, height = feat_h)
  plots[[paste0(sample_name, "  -  Canonical Marker Feature Plots")]] <-
    set_page(p_feat, pw = 11, ph = min(feat_h * 11 / 16, 10))

  message("  Done. Outputs in: ", outdir)
  list(seu = seu, plots = plots)
}

# =============================================================================
# MAIN
# =============================================================================

report_plots <- list()

for (nm in SAMPLE_NAMES) {
  seu <- readRDS(file.path(DIRS$individual, nm, paste0(nm, "_singlets.rds")))
  res <- process_individual(seu, nm)
  report_plots <- c(report_plots, res$plots)
  saveRDS(res$seu, file.path(DIRS$individual, nm, paste0(nm, "_seurat.rds")))
}

save_report_pdf(report_plots, file.path(DIRS$individual, "individual_report.pdf"))
message("\n03_individual.R complete.")
