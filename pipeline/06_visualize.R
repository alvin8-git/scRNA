# =============================================================================
# 06_visualize.R - Publication-quality figures from the annotated object.
#   Outputs: individual PDFs + visualization_report.pdf
# =============================================================================
source("/data/alvin/scRNA/pipeline/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(pheatmap)
  library(viridis)
  library(dplyr)
  library(scales)
})

plan("multicore", workers = PARALLEL$workers)
options(future.globals.maxSize = PARALLEL$future_mem_gb * 1024^3)
message("Parallelism: ", PARALLEL$workers, " cores (future multicore)")

merged <- readRDS(file.path(DIRS$integrated, "integrated_annotated.rds"))
Idents(merged) <- "cell_type"
message("Loaded: ", ncol(merged), " cells | ",
        length(unique(merged$cell_type)), " cell types")

theme_pub <- theme_classic(base_size = 11) +
  theme(plot.title  = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 9),
        axis.text   = element_text(size = 9))

# Resolve cell type colours  -  Unassigned/Unknown forced to grey
all_ct  <- sort(unique(merged$cell_type))
ct_cols <- CELLTYPE_COLORS[all_ct[all_ct %in% names(CELLTYPE_COLORS)]]
missing_ct <- setdiff(all_ct, c(names(CELLTYPE_COLORS), "Unassigned", "Unknown"))
if (length(missing_ct) > 0)
  ct_cols <- c(ct_cols, setNames(hue_pal()(length(missing_ct)), missing_ct))
if ("Unassigned" %in% all_ct) ct_cols["Unassigned"] <- "#AAAAAA"
if ("Unknown"    %in% all_ct) ct_cols["Unknown"]    <- "#AAAAAA"

# Resolve sample colours (handles novel sample names)
sample_cols <- SAMPLE_COLORS[names(SAMPLE_COLORS) %in% SAMPLE_NAMES]
missing_s   <- setdiff(SAMPLE_NAMES, names(SAMPLE_COLORS))
if (length(missing_s) > 0)
  sample_cols <- c(sample_cols, setNames(hue_pal()(length(missing_s)), missing_s))

report_plots <- list()

# =============================================================================
# PLOT SET 1: UMAP triptych  -  cluster / sample / cell type
# =============================================================================
message("\n--- Plot Set 1: UMAP triptych ---")

p1 <- DimPlot(merged, group.by = "seurat_clusters", reduction = "umap",
               label = TRUE, label.size = PLOT$label_size, pt.size = PLOT$pt_size) +
  NoLegend() + labs(title = "Clusters") + theme_pub

p2 <- DimPlot(merged, group.by = "sample", reduction = "umap",
               cols = sample_cols, pt.size = PLOT$pt_size) +
  labs(title = "Sample") + theme_pub

p3 <- DimPlot(merged, group.by = "cell_type", reduction = "umap", cols = ct_cols,
               label = TRUE, label.size = PLOT$label_size, repel = TRUE,
               pt.size = PLOT$pt_size) +
  labs(title = "Cell Type") + theme_pub

ggsave(file.path(DIRS$integrated, "integrated_umap_cluster.pdf"),
       p1, width = PLOT$width, height = PLOT$height, dpi = PLOT$dpi)
ggsave(file.path(DIRS$integrated, "integrated_umap_sample.pdf"),
       p2, width = PLOT$width, height = PLOT$height, dpi = PLOT$dpi)
ggsave(file.path(DIRS$integrated, "integrated_umap_celltype.pdf"),
       p3, width = PLOT$width, height = PLOT$height, dpi = PLOT$dpi)

# Triptych: re-create p3 with smaller labels to avoid clutter
p3_trip <- DimPlot(merged, group.by = "cell_type", reduction = "umap", cols = ct_cols,
                    label = TRUE, label.size = 2.5, repel = TRUE,
                    pt.size = PLOT$pt_size) +
  labs(title = "Cell Type") + theme_pub
p_trip <- p1 | p2 | p3_trip
ggsave(file.path(DIRS$integrated, "umap_triptych.pdf"),
       p_trip, width = 24, height = 7, dpi = PLOT$dpi)

report_plots[["UMAP  -  Clusters"]]   <- mark_small(p1)
report_plots[["UMAP  -  Sample"]]     <- mark_small(p2)
report_plots[["UMAP  -  Cell Types"]] <- set_page(p3, 8.5, 7.5)
report_plots[["UMAP Triptych  -  Clusters / Sample / Cell Type"]] <-
  set_page(p_trip, pw = 11, ph = 4)

# =============================================================================
# PLOT SET 2: Split by sample
# =============================================================================
if (!SINGLE_SAMPLE) {
  message("\n--- Plot Set 2: Split by sample ---")
  p_split <- DimPlot(merged, group.by = "cell_type", split.by = "sample",
                      cols = ct_cols, reduction = "umap",
                      pt.size = PLOT$pt_size, ncol = 2) +
    labs(title = "Cell Types  -  split by Sample") + theme_pub
  ggsave(file.path(DIRS$integrated, "umap_split_by_sample.pdf"),
         p_split, width = 16, height = 7, dpi = PLOT$dpi)
  report_plots[["UMAP  -  Cell Types Split by Sample"]] <- set_page(p_split, pw = 11, ph = 5)
}

# =============================================================================
# PLOT SET 3: Feature plots  -  canonical markers
# =============================================================================
message("\n--- Plot Set 3: Feature plots ---")

# Paginated feature plot: 3 cols x 2 rows per page (genes_per_page = 6)
plot_features <- function(genes, title, filename, genes_per_page = 6L) {
  genes <- genes[genes %in% rownames(merged)]
  if (length(genes) == 0) { message("  Skipped ", filename); return(invisible(NULL)) }

  chunks  <- split(genes, ceiling(seq_along(genes) / genes_per_page))
  n_pages <- length(chunks)

  # Build one ggplot per page
  page_plots <- lapply(chunks, function(gchunk) {
    ncols <- min(3L, length(gchunk))
    FeaturePlot(merged, features = gchunk, reduction = "umap", ncol = ncols,
                pt.size = 0.4, order = TRUE, cols = c("lightgrey", "#E64B35")) &
      theme_pub & theme(axis.text = element_blank(), axis.ticks = element_blank(),
                        plot.title = element_text(size = 10))
  })

  # Save multi-page PDF (one page per chunk)
  pdf_path <- file.path(DIRS$integrated, filename)
  pdf(pdf_path, width = 3 * 4.5, height = 2 * 4)
  for (p in page_plots) print(p)
  dev.off()

  # Add each page to report_plots
  for (ci in seq_along(page_plots)) {
    gchunk   <- chunks[[ci]]
    ncols    <- min(3L, length(gchunk))
    nrows    <- ceiling(length(gchunk) / ncols)
    pg_title <- if (n_pages > 1) paste0(title, "  (", ci, "/", n_pages, ")") else title
    report_plots[[pg_title]] <<- set_page(page_plots[[ci]],
                                           pw = min(ncols * 4.5, 11),
                                           ph = min(nrows * 4, 9))
  }
}

plot_features(c(MARKERS$T_pan, MARKERS$CD4_T, MARKERS$CD8_T),
              "Feature Plot  -  T Cell Markers", "feature_T_cells.pdf")
plot_features(MARKERS$NK,
              "Feature Plot  -  NK Cell Markers", "feature_NK.pdf")
plot_features(MARKERS$B_cell,
              "Feature Plot  -  B Cell Markers", "feature_B_cells.pdf")
plot_features(c(MARKERS$CD14_mono, MARKERS$FCGR3A_mono),
              "Feature Plot  -  Monocyte Markers", "feature_monocytes.pdf")
plot_features(c(MARKERS$DC, MARKERS$Platelet),
              "Feature Plot  -  DC & Platelet Markers", "feature_DC_platelet.pdf")
plot_features(MARKERS$Neutrophil,
              "Feature Plot  -  Neutrophil Markers", "feature_neutrophils.pdf")
plot_features(MARKERS$RBC,
              "Feature Plot  -  RBC Markers (Contamination)", "feature_rbc.pdf")
plot_features(ALL_MARKERS,
              "Feature Plot  -  All Canonical Markers", "feature_all_canonical_markers.pdf")

# =============================================================================
# PLOT SET 4: Dot plot  -  cell types × canonical markers
# =============================================================================
message("\n--- Plot Set 4: Dot plot ---")

markers_in_data <- ALL_MARKERS[ALL_MARKERS %in% rownames(merged)]
p_dot <- DotPlot(merged, features = markers_in_data, group.by = "cell_type",
                  dot.scale = 8, dot.min = 0.01) +
  scale_color_viridis_c(option = "plasma") +
  labs(title = "Canonical PBMC Markers × Cell Type", x = NULL, y = NULL) + theme_pub +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))
ggsave(file.path(DIRS$integrated, "integrated_dotplot.pdf"),
       p_dot, width = 16, height = 7, dpi = PLOT$dpi)
report_plots[["Dot Plot  -  Canonical Markers × Cell Type"]] <-
  set_page(p_dot, pw = 11, ph = 6)

# =============================================================================
# PLOT SET 5: Heatmap  -  top 3 markers per cluster
# =============================================================================
message("\n--- Plot Set 5: Heatmap ---")

marker_file <- file.path(DIRS$integrated, "integrated_cluster_markers.csv")
if (file.exists(marker_file)) {
  all_markers_int <- read.csv(marker_file)
  top3_genes <- all_markers_int %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 3, with_ties = FALSE) %>%
    pull(gene) %>% unique()
  top3_genes <- top3_genes[top3_genes %in% rownames(merged)]
} else {
  message("  Marker file not found  -  computing...")
  Idents(merged) <- "seurat_clusters"
  all_markers_int <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.1,
                                     logfc.threshold = 0.25, verbose = FALSE)
  top3_genes <- all_markers_int %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 3, with_ties = FALSE) %>%
    pull(gene) %>% unique()
  Idents(merged) <- "cell_type"
}

cells_use <- unlist(tapply(colnames(merged), merged$cell_type,
                            function(x) sample(x, min(50L, length(x)))),
                    use.names = FALSE)

p_heatmap <- DoHeatmap(subset(merged, cells = cells_use), features = top3_genes,
                        group.by = "cell_type", group.colors = ct_cols,
                        slot = "data", size = 3) +
  scale_fill_viridis_c() +
  theme(axis.text.y = element_text(size = 7))
ggsave(file.path(DIRS$integrated, "integrated_heatmap.pdf"),
       p_heatmap, width = 14, height = 10, dpi = PLOT$dpi)
report_plots[["Heatmap  -  Top 3 Markers per Cluster"]] <-
  set_page(p_heatmap, pw = 11, ph = 8.5)
Idents(merged) <- "cell_type"

# =============================================================================
# PLOT SET 6: Cell type proportions per sample
# =============================================================================
message("\n--- Plot Set 6: Proportion bar charts ---")

prop_df <- merged@meta.data %>%
  group_by(sample, cell_type) %>% summarise(n = n(), .groups = "drop") %>%
  group_by(sample) %>% mutate(prop = n / sum(n))

p_bar <- ggplot(prop_df, aes(x = cell_type, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_wrap(~ sample, nrow = 1) +
  scale_fill_manual(values = ct_cols, guide = "none") +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Cell Type  -  % of Sample", x = NULL, y = "Proportion") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave(file.path(DIRS$integrated, "celltype_proportions_bar.pdf"),
       p_bar, width = max(7, length(unique(prop_df$sample)) * 5), height = 5,
       dpi = PLOT$dpi)

p_bar_n <- ggplot(prop_df, aes(x = cell_type, y = n, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_wrap(~ sample, nrow = 1) +
  scale_fill_manual(values = ct_cols, guide = "none") +
  labs(title = "Cell Type  -  Number of Cells", x = NULL, y = "Cell Count",
       fill = "Cell Type") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave(file.path(DIRS$integrated, "celltype_counts_bar.pdf"),
       p_bar_n, width = max(7, length(unique(prop_df$sample)) * 5), height = 5,
       dpi = PLOT$dpi)

# Combined 2-panel portrait figure: % on top, counts on bottom
p_bar_combined <- p_bar / p_bar_n +
  plot_annotation(
    title = "Cell Type Composition per Sample",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )
ggsave(file.path(DIRS$integrated, "celltype_composition_combined.pdf"),
       p_bar_combined,
       width  = max(7, length(unique(prop_df$sample)) * 5),
       height = 10, dpi = PLOT$dpi)

report_plots[["Cell Type Composition  -  % (top) and Cell Count (bottom)"]] <-
  set_page(p_bar_combined, pw = max(7, length(unique(prop_df$sample)) * 5), ph = 10)

# =============================================================================
# PLOT SET 7: Violin plots  -  key lineage markers
# =============================================================================
message("\n--- Plot Set 7: Violin plots ---")

key_markers <- c("CD3D", "CD4", "CD8A", "NKG7", "MS4A1", "CD14", "FCGR3A", "PPBP")
key_markers <- key_markers[key_markers %in% rownames(merged)]

p_vln <- VlnPlot(merged, features = key_markers, group.by = "cell_type",
                  pt.size = 0, ncol = 4, cols = ct_cols) &
  theme_pub & theme(legend.position = "none",
                    axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
vln_h <- ceiling(length(key_markers) / 4) * 4 + 1
ggsave(file.path(DIRS$integrated, "violin_key_markers.pdf"),
       p_vln, width = 18, height = vln_h, dpi = PLOT$dpi)
report_plots[["Violin  -  Key Lineage Markers per Cell Type"]] <-
  set_page(p_vln, pw = 11, ph = min(vln_h * 11 / 18, 8.5))

# =============================================================================
# SUMMARY
# =============================================================================
ct_summary <- merged@meta.data %>%
  group_by(cell_type) %>%
  summarise(total = n(), .groups = "drop") %>%
  arrange(desc(total))

message("\nFinal cell type composition:")
print(as.data.frame(ct_summary))

save_report_pdf(report_plots, file.path(DIRS$integrated, "visualization_report.pdf"))
message("\n--- All plots saved to: ", DIRS$integrated, " ---")
message("\n06_visualize.R complete.")
