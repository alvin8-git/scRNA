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

# SingleR labels UMAP  -  Unassigned forced to grey, known types use CELLTYPE_COLORS
singler_labels <- unique(merged$singler_label_clean)
singler_cols   <- setNames(
  ifelse(singler_labels == "Unassigned", "#AAAAAA",
         ifelse(singler_labels %in% names(CELLTYPE_COLORS),
                CELLTYPE_COLORS[singler_labels],
                scales::hue_pal()(length(singler_labels))
                  [seq_along(singler_labels)])),
  singler_labels)
p_singler <- DimPlot(merged, group.by = "singler_label_clean", reduction = "umap",
                      cols = singler_cols,
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
  # Assign per-cell SingleR labels directly so rare types (Neutrophils, Unassigned)
  # are not absorbed into the cluster majority and are visible in the cell type UMAP
  merged$cell_type <- merged$singler_label_clean
  message("  Using per-cell SingleR labels (set CLUSTER_CELLTYPE_MAP in config.R for manual labels)")
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

# Cell type UMAP  -  Unassigned → grey, known types use CELLTYPE_COLORS, rest auto-coloured
ct_labels      <- unique(merged$cell_type)
ct_colors_used <- CELLTYPE_COLORS[names(CELLTYPE_COLORS) %in% ct_labels]
extra          <- setdiff(ct_labels, c(names(CELLTYPE_COLORS), "Unassigned", "Unknown"))
if (length(extra) > 0)
  ct_colors_used <- c(ct_colors_used, setNames(hue_pal()(length(extra)), extra))
if ("Unassigned" %in% ct_labels) ct_colors_used["Unassigned"] <- "#AAAAAA"
if ("Unknown"    %in% ct_labels) ct_colors_used["Unknown"]    <- "#AAAAAA"

p_ct <- DimPlot(merged, group.by = "cell_type", reduction = "umap",
                 cols = ct_colors_used, label = TRUE,
                 label.size = PLOT$label_size, repel = TRUE, pt.size = PLOT$pt_size) +
  labs(title = "Integrated  -  Cell Types") + theme_classic()
ggsave(file.path(DIRS$annotation, "celltype_umap.pdf"),
       p_ct, width = PLOT$width + 2, height = PLOT$height, dpi = PLOT$dpi)
report_plots[["Integrated  -  Cell Type UMAP"]] <- set_page(p_ct, 8.5, 7.5)

# =============================================================================
# PART 4: DC and monocyte sub-type presence check
# =============================================================================
message("\n--- DC & monocyte sub-type check ---")
check_markers <- function(genes, label) {
  present <- genes[genes %in% rownames(merged)]
  if (length(present) == 0) {
    message("  ", label, ": marker genes not found in dataset — may be absent or filtered")
    return(invisible(NULL))
  }
  expr <- GetAssayData(merged, assay = "RNA", layer = "data")[present, , drop = FALSE]
  pct_expressing <- round(rowMeans(expr > 0) * 100, 1)
  message("  ", label, ": ", paste0(present, " (", pct_expressing, "% cells +ve)", collapse = ", "))
}
check_markers(MARKERS$DC,          "Dendritic cells  (FCER1A, CLEC9A)")
check_markers(MARKERS$CD14_mono,   "CD14+ Monocytes  (CD14, LYZ, CST3, S100A8)")
check_markers(MARKERS$FCGR3A_mono, "FCGR3A+ Monocytes (FCGR3A, MS4A7)")

# Confirm monocyte sub-types are resolved in separate clusters
mono_genes <- c("CD14", "FCGR3A")
mono_present <- mono_genes[mono_genes %in% rownames(merged)]
if (length(mono_present) == 2) {
  expr_mat   <- GetAssayData(merged, assay = "RNA", layer = "data")[mono_genes, ]
  meta       <- merged@meta.data
  meta$cd14_hi   <- as.numeric(expr_mat["CD14",   ] > 0)
  meta$fcgr3a_hi <- as.numeric(expr_mat["FCGR3A", ] > 0)
  per_cluster <- meta %>%
    group_by(seurat_clusters) %>%
    summarise(pct_CD14   = round(mean(cd14_hi)   * 100, 1),
              pct_FCGR3A = round(mean(fcgr3a_hi) * 100, 1),
              n = n(), .groups = "drop") %>%
    filter(pct_CD14 > 10 | pct_FCGR3A > 10) %>%
    arrange(desc(pct_CD14))
  if (nrow(per_cluster) > 0) {
    message("  Monocyte marker expression by cluster:")
    print(as.data.frame(per_cluster))
    if (n_distinct(per_cluster$seurat_clusters) < 2)
      message("  NOTE: CD14+ and FCGR3A+ monocytes may share a cluster — ",
              "try CLUSTER$default_res = 0.6 or 0.8")
    else
      message("  CD14+ and FCGR3A+ monocytes appear in separate clusters — looks good.")
  }
}

# =============================================================================
# PART 5: T cell sub-clustering
# =============================================================================
if (isTRUE(SUBCLUSTER$enabled)) {
  message("\n--- T cell sub-clustering (res=", SUBCLUSTER$resolution, ") ---")

  # Identify clusters whose majority cell type matches T cell patterns
  t_clusters <- annot_table %>%
    filter(grepl(SUBCLUSTER$t_patterns, final_cell_type,
                 ignore.case = TRUE, perl = TRUE)) %>%
    pull(seurat_clusters) %>% as.character()

  n_t_cells <- sum(as.character(merged$seurat_clusters) %in% t_clusters)
  message("  T cell clusters identified: ",
          if (length(t_clusters)) paste(t_clusters, collapse = ", ") else "none",
          "  (", n_t_cells, " cells)")

  if (length(t_clusters) == 0) {
    message("  No clusters matched SUBCLUSTER$t_patterns — ",
            "set CLUSTER_CELLTYPE_MAP or adjust t_patterns in config.R")
  } else if (n_t_cells < SUBCLUSTER$min_cells) {
    message("  < ", SUBCLUSTER$min_cells, " T cells — skipping sub-clustering")
  } else {
    # Find the SNN graph name
    graph_name <- grep("_snn$", names(merged@graphs), value = TRUE)[1]
    if (is.na(graph_name)) graph_name <- "RNA_snn"
    message("  Using graph: ", graph_name)

    merged <- FindSubCluster(merged, cluster = t_clusters,
                              graph.name = graph_name,
                              resolution = SUBCLUSTER$resolution,
                              algorithm  = CLUSTER$algorithm)

    sub_labels <- merged$sub.cluster[as.character(merged$seurat_clusters) %in% t_clusters]
    n_sub <- length(unique(sub_labels))
    message("  Sub-clusters found: ", n_sub)

    # UMAP — non-T cells grey, T sub-clusters coloured
    sub_col_vals <- c(Other = "#DDDDDD",
                      setNames(scales::hue_pal()(n_sub), sort(unique(sub_labels))))
    merged$sub_plot <- ifelse(as.character(merged$seurat_clusters) %in% t_clusters,
                               merged$sub.cluster, "Other")

    p_sub <- DimPlot(merged, group.by = "sub_plot", reduction = "umap",
                     cols = sub_col_vals,
                     label = TRUE, label.size = 3, pt.size = PLOT$pt_size, repel = TRUE) +
      labs(title = paste0("T Cell Sub-clusters  (res=", SUBCLUSTER$resolution, ")")) +
      theme_classic()
    ggsave(file.path(DIRS$annotation, "tcell_subclusters_umap.pdf"),
           p_sub, width = 9, height = 7, dpi = PLOT$dpi)
    report_plots[["T Cell Sub-clusters  -  UMAP"]] <- set_page(p_sub, pw = 9, ph = 7)

    # Dot plot of T cell markers across sub-clusters
    t_markers <- unique(c(MARKERS$T_pan, MARKERS$CD4_T, MARKERS$CD8_T, MARKERS$Treg))
    t_markers <- t_markers[t_markers %in% rownames(merged)]
    t_barcodes <- colnames(merged)[as.character(merged$seurat_clusters) %in% t_clusters]
    merged_t   <- subset(merged, cells = t_barcodes)
    Idents(merged_t) <- "sub.cluster"

    if (length(t_markers) > 0 && length(unique(Idents(merged_t))) > 1) {
      p_sub_dot <- DotPlot(merged_t, features = t_markers,
                            dot.scale = 8, dot.min = 0.01) +
        scale_color_viridis_c(option = "plasma") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)) +
        labs(title = "T Cell Sub-clusters  -  CD4 / CD8 / Treg Markers",
             x = NULL, y = "Sub-cluster")
      ggsave(file.path(DIRS$annotation, "tcell_subclusters_dotplot.pdf"),
             p_sub_dot, width = 10, height = 6, dpi = PLOT$dpi)
      report_plots[["T Cell Sub-clusters  -  CD4 / CD8 / Treg Marker Dot Plot"]] <-
        set_page(p_sub_dot, pw = 10, ph = 6)
    }

    # Summary table
    sub_summary <- merged@meta.data %>%
      filter(as.character(seurat_clusters) %in% t_clusters) %>%
      group_by(sub.cluster) %>%
      summarise(n_cells          = n(),
                singler_majority = names(which.max(table(singler_label_clean))),
                .groups = "drop")
    write.csv(sub_summary,
              file.path(DIRS$annotation, "tcell_subcluster_summary.csv"),
              row.names = FALSE)
    message("  Sub-cluster summary:")
    print(as.data.frame(sub_summary))
  }
}

saveRDS(merged, file.path(DIRS$integrated, "integrated_annotated.rds"))

save_report_pdf(report_plots, file.path(DIRS$annotation, "annotation_report.pdf"))

message("\n05_annotate.R complete.")
message("Key output: ", file.path(DIRS$annotation, "canonical_markers_dotplot.pdf"))
