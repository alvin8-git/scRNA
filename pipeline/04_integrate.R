# =============================================================================
# 04_integrate.R - Merge samples, run Harmony batch correction (multi-sample)
#                  or direct UMAP/cluster (single-sample), generate integration QC.
#                  Outputs: integrated_seurat.rds + integration_report.pdf
# =============================================================================
source("/data/alvin/scRNA/pipeline/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(future.apply)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  if (!SINGLE_SAMPLE) library(harmony)
})

plan("multicore", workers = PARALLEL$workers)
options(future.globals.maxSize = PARALLEL$future_mem_gb * 1024^3)
message("Parallelism: ", PARALLEL$workers, " cores (future multicore)")

report_plots <- list()

# Load all per-sample processed objects
seu_list <- setNames(
  lapply(SAMPLE_NAMES, function(nm)
    readRDS(file.path(DIRS$individual, nm, paste0(nm, "_seurat.rds")))),
  SAMPLE_NAMES
)

# =============================================================================
# SINGLE-SAMPLE PATH  -  no batch correction needed
# =============================================================================
if (SINGLE_SAMPLE) {
  nm     <- SAMPLE_NAMES[1]
  merged <- seu_list[[1]]
  message("Single-sample mode: ", nm, " (", ncol(merged), " cells)")

  merged <- NormalizeData(merged, verbose = FALSE)
  merged <- FindVariableFeatures(merged, nfeatures = NORM$n_hvg, verbose = FALSE)
  merged <- ScaleData(merged, vars.to.regress = "percent.mt", verbose = FALSE)
  merged <- RunPCA(merged, npcs = DIM$npcs, verbose = FALSE)
  merged <- RunUMAP(merged, dims = DIM$dims_use, seed.use = DIM$umap_seed,
                    umap.method = "uwot", verbose = FALSE)
  merged <- FindNeighbors(merged, dims = DIM$dims_use, verbose = FALSE)

  for (res in CLUSTER$resolutions) {
    merged <- FindClusters(merged, resolution = res,
                           algorithm = CLUSTER$algorithm, verbose = FALSE)
    message("  res=", res, ": ", length(unique(merged$seurat_clusters)), " clusters")
  }
  default_col <- paste0("RNA_snn_res.", CLUSTER$default_res)
  Idents(merged) <- default_col
  merged$seurat_clusters <- merged@meta.data[[default_col]]

  p_cluster <- DimPlot(merged, reduction = "umap", label = TRUE,
                        label.size = PLOT$label_size, pt.size = PLOT$pt_size) +
    labs(title = paste0(nm, "  -  Clusters (res=", CLUSTER$default_res, ")")) +
    theme_classic()
  ggsave(file.path(DIRS$integrated, "integrated_umap_cluster.pdf"),
         p_cluster, width = PLOT$width, height = PLOT$height, dpi = PLOT$dpi)
  report_plots[["Clusters"]] <- set_page(p_cluster, 8.5, 7.5)

# =============================================================================
# MULTI-SAMPLE PATH  -  Harmony batch correction
# =============================================================================
} else {
  message(paste(SAMPLE_NAMES, collapse = " + "), " | ",
          paste(sapply(seu_list, ncol), collapse = " + "), " cells")

  # Step 1: Merge
  merged <- suppressWarnings(
    merge(seu_list[[1]],
          y            = seu_list[-1],
          add.cell.ids = NULL,
          merge.data   = FALSE,
          project      = "integrated")
  )
  merged <- suppressWarnings(JoinLayers(merged))
  message("Merged: ", ncol(merged), " total cells")

  # Step 2: Normalize, HVG, Scale, PCA
  merged <- NormalizeData(merged, verbose = FALSE)
  merged <- FindVariableFeatures(merged, selection.method = NORM$hvg_method,
                                  nfeatures = NORM$n_hvg, verbose = FALSE)
  merged <- ScaleData(merged, vars.to.regress = "percent.mt", verbose = FALSE)
  merged <- RunPCA(merged, npcs = DIM$npcs, verbose = FALSE)

  # Step 3: UMAP before Harmony
  merged <- RunUMAP(merged, dims = DIM$dims_use,
                    reduction.name = "umap_uncorrected", reduction.key = "UMAPunc_",
                    seed.use = DIM$umap_seed, umap.method = "uwot", verbose = FALSE)

  sample_cols <- SAMPLE_COLORS[names(SAMPLE_COLORS) %in% SAMPLE_NAMES]
  if (length(sample_cols) < length(SAMPLE_NAMES)) {
    missing <- setdiff(SAMPLE_NAMES, names(SAMPLE_COLORS))
    sample_cols <- c(sample_cols, setNames(scales::hue_pal()(length(missing)), missing))
  }

  p_before <- DimPlot(merged, reduction = "umap_uncorrected",
                       group.by = "sample", cols = sample_cols, pt.size = PLOT$pt_size) +
    labs(title = "Before Harmony  -  by sample") + theme_classic()

  # Step 4: Harmony
  message("Running Harmony integration...")
  merged <- RunHarmony(merged,
    group.by.vars    = HARMONY$group_by_vars,
    dims.use         = DIM$dims_use,
    theta            = HARMONY$theta,
    lambda           = HARMONY$lambda,
    nclust           = HARMONY$nclust,
    max.iter.harmony = HARMONY$max_iter,
    reduction        = "pca", reduction.save = "harmony",
    plot_convergence = FALSE, verbose = FALSE
  )

  # Step 5: UMAP on Harmony embedding
  merged <- RunUMAP(merged, reduction = "harmony", dims = HARMONY$dims_use,
                    reduction.name = "umap", reduction.key = "UMAP_",
                    seed.use = DIM$umap_seed, umap.method = "uwot", verbose = FALSE)

  p_after <- DimPlot(merged, reduction = "umap",
                      group.by = "sample", cols = sample_cols, pt.size = PLOT$pt_size) +
    labs(title = "After Harmony  -  by sample") + theme_classic()

  p_compare <- p_before | p_after
  ggsave(file.path(DIRS$integrated, "harmony_before_after.pdf"),
         p_compare, width = 16, height = 7, dpi = PLOT$dpi)
  report_plots[["Harmony Integration  -  Before vs After (by Sample)"]] <-
    set_page(p_compare, pw = 11, ph = 5)

  # Step 6: Cluster on Harmony embedding
  merged <- FindNeighbors(merged, reduction = "harmony",
                           dims = HARMONY$dims_use, verbose = FALSE)
  for (res in CLUSTER$resolutions) {
    merged <- FindClusters(merged, resolution = res,
                           algorithm = CLUSTER$algorithm, verbose = FALSE)
    message("  res=", res, ": ", length(unique(merged$seurat_clusters)), " clusters")
  }
  default_col <- paste0("RNA_snn_res.", CLUSTER$default_res)
  Idents(merged) <- default_col
  merged$seurat_clusters <- merged@meta.data[[default_col]]

  # Step 7: Core UMAP plots
  p_cluster <- DimPlot(merged, reduction = "umap", label = TRUE,
                        label.size = PLOT$label_size, pt.size = PLOT$pt_size) +
    labs(title = paste0("Integrated  -  Clusters (res=", CLUSTER$default_res, ")")) +
    theme_classic()
  ggsave(file.path(DIRS$integrated, "integrated_umap_cluster.pdf"),
         p_cluster, width = PLOT$width, height = PLOT$height, dpi = PLOT$dpi)
  report_plots[["Integrated  -  UMAP Clusters"]] <- set_page(p_cluster, 8.5, 7.5)

  p_sample <- DimPlot(merged, reduction = "umap",
                       group.by = "sample", cols = sample_cols, pt.size = PLOT$pt_size) +
    labs(title = "Integrated  -  by Sample") + theme_classic()
  ggsave(file.path(DIRS$integrated, "integrated_umap_sample.pdf"),
         p_sample, width = PLOT$width, height = PLOT$height, dpi = PLOT$dpi)
  report_plots[["Integrated  -  UMAP by Sample"]] <- mark_small(p_sample)

  p_split <- DimPlot(merged, reduction = "umap", split.by = "sample",
                      pt.size = PLOT$pt_size, ncol = 2) +
    labs(title = "Integrated  -  Split by Sample") + theme_classic()
  ggsave(file.path(DIRS$integrated, "integrated_umap_split_sample.pdf"),
         p_split, width = 14, height = 7, dpi = PLOT$dpi)
  report_plots[["Integrated  -  UMAP Split by Sample"]] <- set_page(p_split, pw = 11, ph = 5.5)
}

# =============================================================================
# Cluster markers (both paths)
# =============================================================================
message("Finding integrated cluster markers...")
integrated_markers <- FindAllMarkers(merged, only.pos = TRUE,
                                      min.pct = 0.1, logfc.threshold = 0.25,
                                      verbose = FALSE)
write.csv(integrated_markers,
          file.path(DIRS$integrated, "integrated_cluster_markers.csv"),
          row.names = FALSE)

saveRDS(merged, file.path(DIRS$integrated, "integrated_seurat.rds"))

save_report_pdf(report_plots, file.path(DIRS$integrated, "integration_report.pdf"))

message("\n04_integrate.R complete.")
message("Object: ", ncol(merged), " cells, ",
        length(unique(merged$seurat_clusters)), " clusters")
