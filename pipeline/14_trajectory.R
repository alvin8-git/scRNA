# =============================================================================
# 14_trajectory.R — Monocle3 pseudotime trajectory (bat_wing mode)
# Requires: integrated_annotated.rds from step 04+05
# Outputs:  pseudotime UMAPs, trajectory graph, top pseudotime DEGs,
#           trajectory_report.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(monocle3); library(Seurat); library(SeuratWrappers)
  library(dplyr); library(ggplot2); library(patchwork)
  library(pdftools); library(magick)
})
source(file.path(dirname(sys.frame(1)$ofilename %||% "."), "config.R"))

dir.create(DIRS$trajectory, showWarnings = FALSE, recursive = TRUE)

.combine_pdfs <- function(paths, out) {
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) return(invisible(NULL))
  imgs <- lapply(paths, function(p) {
    n <- tryCatch(pdftools::pdf_length(p), error = function(e) 1L)
    lapply(seq_len(n), function(i) magick::image_read_pdf(p, pages = i, density = 150))
  })
  magick::image_write(magick::image_join(unlist(imgs, recursive = FALSE)), out, format = "pdf")
}

rds_path <- file.path(DIRS$integrated, "integrated_annotated.rds")
if (!file.exists(rds_path))
  stop("Missing: ", rds_path, " — run steps 04+05 first")
merged <- readRDS(rds_path)
merged$condition <- SAMPLE_CONDITIONS[merged$sample]

# =============================================================================
# PART 1: Full dataset trajectory
# =============================================================================
message("Converting Seurat → Monocle3 CDS...")
cds <- as.cell_data_set(merged)
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)

# Root: healthy cells from the most abundant healthy cell type
healthy_cells <- colnames(merged)[merged$condition == CONDITION_LEVELS[1]]
if (length(healthy_cells) == 0) healthy_cells <- colnames(merged)
root_type  <- names(sort(table(merged$cell_type[healthy_cells]), decreasing = TRUE))[1]
root_cells <- healthy_cells[merged$cell_type[healthy_cells] == root_type]
root_cells <- head(root_cells, 50)

cds <- order_cells(cds, root_cells = root_cells)
saveRDS(cds, file.path(DIRS$trajectory, "monocle3_cds.rds"))

p_pt   <- plot_cells(cds, color_cells_by = "pseudotime",
                     show_trajectory_graph = TRUE, cell_size = 0.5,
                     label_cell_groups = FALSE) +
  scale_colour_viridis_c(option = "C") + theme_classic(base_size = 10)
p_ct   <- plot_cells(cds, color_cells_by = "cell_type",
                     show_trajectory_graph = TRUE, cell_size = 0.5,
                     label_groups_by_cluster = FALSE, group_label_size = 3) +
  theme_classic(base_size = 10)
p_cond <- plot_cells(cds, color_cells_by = "condition",
                     show_trajectory_graph = FALSE, cell_size = 0.5,
                     label_cell_groups = FALSE) +
  scale_colour_manual(values = c("steelblue", "tomato")) +
  theme_classic(base_size = 10)

full_pdf <- file.path(DIRS$trajectory, "pseudotime_full.pdf")
pdf(full_pdf, width = 14, height = 5)
print(p_pt | p_ct | p_cond)
dev.off()
message("Saved: ", full_pdf)

# =============================================================================
# PART 2: Lineage-specific trajectories
# =============================================================================
lineages <- list(
  Fibroblast = c("Fibroblast", "Myofibroblast", "Fibroblast (resting)",
                 "Fibroblast (wound)"),
  Macrophage = c("Macrophage", "Macrophage (M1/inflam)", "Macrophage (M2/repair)",
                 "Macrophage (proliferat)")
)

lineage_pdfs <- character(0)
for (lin_name in names(lineages)) {
  lin_types <- lineages[[lin_name]]
  lin_cells <- colnames(merged)[merged$cell_type %in% lin_types]
  if (length(lin_cells) < 50) {
    message("Skipping ", lin_name, " — too few cells (", length(lin_cells), ")")
    next
  }
  message("Trajectory: ", lin_name, " (", length(lin_cells), " cells)")

  sub_seu <- subset(merged, cells = lin_cells)
  sub_cds <- tryCatch({
    sc <- as.cell_data_set(sub_seu)
    sc <- cluster_cells(sc, reduction_method = "UMAP")
    sc <- learn_graph(sc, use_partition = FALSE, verbose = FALSE)

    root_lin <- intersect(lin_cells,
      colnames(merged)[merged$condition == CONDITION_LEVELS[1] &
                       merged$cell_type == lin_types[1]])
    if (length(root_lin) < 5) root_lin <- head(lin_cells, 20)
    sc <- order_cells(sc, root_cells = head(root_lin, 30))
    sc
  }, error = function(e) { message("  Lineage CDS failed: ", e$message); NULL })

  if (is.null(sub_cds)) next

  p1 <- plot_cells(sub_cds, color_cells_by = "pseudotime",
                   show_trajectory_graph = TRUE, cell_size = 1.2,
                   label_cell_groups = FALSE) +
    scale_colour_viridis_c(option = "C") +
    labs(title = paste0(lin_name, " — pseudotime")) + theme_classic()
  p2 <- plot_cells(sub_cds, color_cells_by = "cell_type",
                   show_trajectory_graph = TRUE, cell_size = 1.2,
                   group_label_size = 4) +
    labs(title = paste0(lin_name, " — cell type")) + theme_classic()
  p3 <- plot_cells(sub_cds, color_cells_by = "condition",
                   show_trajectory_graph = FALSE, cell_size = 1.2,
                   label_cell_groups = FALSE) +
    scale_colour_manual(values = c("steelblue", "tomato")) +
    labs(title = paste0(lin_name, " — condition")) + theme_classic()

  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp, width = 14, height = 5)
  print(p1 | p2 | p3)
  dev.off()
  lineage_pdfs <- c(lineage_pdfs, tmp)

  saveRDS(sub_cds, file.path(DIRS$trajectory, paste0("cds_", lin_name, ".rds")))
}

# =============================================================================
# PART 3: Pseudotime-associated DEGs (Moran's I graph autocorrelation)
# =============================================================================
message("Testing pseudotime-associated genes (graph_test)...")
pt_degs_pdfs <- character(0)
tryCatch({
  pr_test <- graph_test(cds, neighbor_graph = "principal_graph",
                        cores = min(4L, PARALLEL$workers))
  pt_sig  <- pr_test %>%
    dplyr::filter(q_value < 0.05) %>%
    dplyr::arrange(dplyr::desc(morans_I)) %>%
    head(50)
  write.csv(pt_sig,
            file.path(DIRS$trajectory, "pseudotime_DEG.csv"), row.names = FALSE)
  message("Saved pseudotime_DEG.csv (", nrow(pt_sig), " genes)")

  top_genes <- head(pt_sig$gene_short_name, 9)
  if (length(top_genes) >= 3) {
    p_gene <- plot_cells(cds, genes = top_genes,
                         show_trajectory_graph = FALSE,
                         label_cell_groups = FALSE, cell_size = 0.5) +
      theme_classic(base_size = 8)
    tmp <- tempfile(fileext = ".pdf")
    ggsave(tmp, p_gene, width = 12, height = 9)
    pt_degs_pdfs <- c(pt_degs_pdfs, tmp)
  }
}, error = function(e) message("graph_test failed: ", e$message))

all_pdfs <- c(full_pdf, lineage_pdfs, pt_degs_pdfs)
.combine_pdfs(all_pdfs, file.path(DIRS$trajectory, "trajectory_report.pdf"))

message("Step 14 complete — trajectory results in: ", DIRS$trajectory)
