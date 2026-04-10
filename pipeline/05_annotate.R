# =============================================================================
# 05_annotate.R - Cell type annotation:
#   Part 1: SingleR automated annotation
#   Part 2: Canonical marker dot plot for manual annotation
#   Part 3: Assign final cell_type labels
#   Outputs: integrated_annotated.rds + annotation_report.pdf
# =============================================================================
source("/data/alvin/scRNA/pipeline/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(BiocParallel)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(dplyr)
  library(scales)
})

bp_param <- MulticoreParam(workers = PARALLEL$workers)
message("Parallelism: ", PARALLEL$workers, " cores (BiocParallel)")

merged <- readRDS(file.path(DIRS$integrated, "integrated_seurat.rds"))
message("Loaded: ", ncol(merged), " cells, ",
        length(unique(merged$seurat_clusters)), " clusters")

report_plots <- list()

# =============================================================================
# PART 1: SingleR automated annotation
# =============================================================================
message("\n--- Running SingleR ---")

ref            <- HumanPrimaryCellAtlasData()
expr_mat       <- GetAssayData(merged, assay = "RNA", layer = "data")
singler_result <- SingleR(test = expr_mat, ref = ref, labels = ref$label.main,
                           fine.tune = TRUE, prune = TRUE, BPPARAM = bp_param)

merged$singler_label  <- singler_result$labels
merged$singler_pruned <- singler_result$pruned.labels
merged$singler_delta  <- apply(singler_result$scores, 1,
                                function(x) { s <- sort(x, decreasing = TRUE); s[1] - s[2] })
merged$singler_label_clean <- ifelse(is.na(merged$singler_pruned),
                                      "Unassigned", merged$singler_pruned)

# SingleR score heatmap (base R via pheatmap  -  save file, include in report)
singler_heatmap_path <- file.path(DIRS$annotation, "singler_scores_heatmap.pdf")
pdf(singler_heatmap_path, width = 12, height = 8)
plotScoreHeatmap(singler_result, fontsize_row = 6,
  annotation_col = data.frame(Cluster = as.character(merged$seurat_clusters),
                               row.names = colnames(merged)))
dev.off()
report_plots[["SingleR  -  Score Heatmap (per-cell annotation confidence)"]] <-
  singler_heatmap_path

# SingleR labels UMAP
p_singler <- DimPlot(merged, group.by = "singler_label_clean", reduction = "umap",
                      label = TRUE, label.size = 3, pt.size = PLOT$pt_size, repel = TRUE) +
  labs(title = "SingleR  -  HumanPrimaryCellAtlasData") + theme_classic() +
  theme(legend.position = "right", legend.text = element_text(size = 8))
ggsave(file.path(DIRS$annotation, "singler_labels_umap.pdf"),
       p_singler, width = 10, height = 7, dpi = PLOT$dpi)
report_plots[["SingleR  -  Labels on UMAP"]] <- set_page(p_singler, 8.5, 7.5)

# Delta score UMAP (confidence)
p_delta <- FeaturePlot(merged, features = "singler_delta", reduction = "umap",
                        pt.size = PLOT$pt_size, cols = c("lightgrey", "#3C5488")) +
  labs(title = "SingleR  -  Annotation Confidence (Delta Score)") + theme_classic()
ggsave(file.path(DIRS$annotation, "singler_delta_umap.pdf"),
       p_delta, width = 8, height = 7, dpi = PLOT$dpi)
report_plots[["SingleR  -  Confidence Delta Score on UMAP"]] <- mark_small(p_delta)

# Cluster-level majority vote
cluster_singler <- merged@meta.data %>%
  group_by(seurat_clusters, singler_label_clean) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  rename(majority_singler = singler_label_clean)
message("  SingleR majority labels per cluster:")
print(as.data.frame(cluster_singler))

# =============================================================================
# PART 2: Canonical marker plots for manual annotation
# =============================================================================
message("\n--- Canonical marker plots ---")

markers_present <- lapply(MARKERS, function(g) g[g %in% rownames(merged)])

p_dot_manual <- DotPlot(merged, features = unique(unlist(markers_present)),
                         group.by = "seurat_clusters", dot.scale = 8, dot.min = 0.01) +
  scale_color_viridis_c(option = "plasma") + RotatedAxis() +
  labs(title = "Canonical PBMC Markers  -  use this to fill CLUSTER_CELLTYPE_MAP",
       x = NULL, y = "Cluster")
ggsave(file.path(DIRS$annotation, "canonical_markers_dotplot.pdf"),
       p_dot_manual, width = 16, height = 8, dpi = PLOT$dpi)
report_plots[["Canonical PBMC Markers Dot Plot (for manual cluster annotation)"]] <-
  set_page(p_dot_manual, pw = 11, ph = 7)

plot_marker_group <- function(genes, title, filename) {
  genes <- genes[genes %in% rownames(merged)]
  if (length(genes) == 0) return(invisible(NULL))
  ncols <- min(3, length(genes))
  nrows <- ceiling(length(genes) / ncols)
  p <- FeaturePlot(merged, features = genes, reduction = "umap", ncol = ncols,
                    pt.size = 0.4, order = TRUE, cols = c("lightgrey", "#E64B35")) &
    theme(plot.title = element_text(size = 9),
          axis.text = element_blank(), axis.ticks = element_blank())
  ggsave(file.path(DIRS$annotation, filename), p,
         width = ncols * 5, height = nrows * 4, dpi = PLOT$dpi)
  report_plots[[title]] <<- set_page(p, pw = min(ncols * 5, 11),
                                      ph = min(nrows * 4, 9))
}

plot_marker_group(c(MARKERS$T_pan, MARKERS$CD4_T, MARKERS$CD8_T),
                  "Feature Plot  -  T Cell Markers (CD3D, CD4, CD8A, CCR7, GZMK…)",
                  "feature_T_cells.pdf")
plot_marker_group(MARKERS$NK,
                  "Feature Plot  -  NK Cell Markers (NKG7, GNLY, KLRD1)",
                  "feature_NK.pdf")
plot_marker_group(MARKERS$B_cell,
                  "Feature Plot  -  B Cell Markers (MS4A1, CD79A, CD19)",
                  "feature_B_cells.pdf")
plot_marker_group(c(MARKERS$CD14_mono, MARKERS$FCGR3A_mono),
                  "Feature Plot  -  Monocyte Markers (CD14, LYZ, FCGR3A, MS4A7…)",
                  "feature_monocytes.pdf")
plot_marker_group(c(MARKERS$DC, MARKERS$Platelet),
                  "Feature Plot  -  DC & Platelet Markers (FCER1A, CLEC9A, PPBP, PF4)",
                  "feature_DC_platelet.pdf")

# =============================================================================
# PART 3: Assign final cell_type labels
# =============================================================================
message("\n--- Assigning cell type labels ---")

if (!is.null(CLUSTER_CELLTYPE_MAP)) {
  merged$cell_type <- unname(CLUSTER_CELLTYPE_MAP[as.character(merged$seurat_clusters)])
  merged$cell_type[is.na(merged$cell_type)] <- "Unknown"
  message("  Using manual CLUSTER_CELLTYPE_MAP")
} else {
  cluster_map <- setNames(cluster_singler$majority_singler,
                           as.character(cluster_singler$seurat_clusters))
  merged$cell_type <- unname(cluster_map[as.character(merged$seurat_clusters)])
  merged$cell_type[is.na(merged$cell_type)] <- "Unassigned"
  message("  Using SingleR majority labels (set CLUSTER_CELLTYPE_MAP in config.R for manual labels)")
}

annot_table <- merged@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_cells         = n(),
    final_cell_type = names(which.max(table(cell_type))),
    singler_majority = names(which.max(table(singler_label_clean))),
    .groups = "drop"
  )
write.csv(annot_table, file.path(DIRS$annotation, "cluster_annotation_table.csv"),
          row.names = FALSE)
message("\nCluster annotation summary:")
print(as.data.frame(annot_table))

# Cell type UMAP
ct_colors_used <- CELLTYPE_COLORS[names(CELLTYPE_COLORS) %in% unique(merged$cell_type)]
extra <- setdiff(unique(merged$cell_type), names(CELLTYPE_COLORS))
if (length(extra) > 0)
  ct_colors_used <- c(ct_colors_used, setNames(hue_pal()(length(extra)), extra))

p_ct <- DimPlot(merged, group.by = "cell_type", reduction = "umap",
                 cols = ct_colors_used, label = TRUE,
                 label.size = PLOT$label_size, repel = TRUE, pt.size = PLOT$pt_size) +
  labs(title = "Integrated  -  Cell Types") + theme_classic()
ggsave(file.path(DIRS$annotation, "celltype_umap.pdf"),
       p_ct, width = PLOT$width + 2, height = PLOT$height, dpi = PLOT$dpi)
report_plots[["Integrated  -  Cell Type UMAP"]] <- set_page(p_ct, 8.5, 7.5)

saveRDS(merged, file.path(DIRS$integrated, "integrated_annotated.rds"))

save_report_pdf(report_plots, file.path(DIRS$annotation, "annotation_report.pdf"))

message("\n05_annotate.R complete.")
message("Key output: ", file.path(DIRS$annotation, "canonical_markers_dotplot.pdf"))
