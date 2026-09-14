# =============================================================================
# trajectory.R — Monocle3 pseudotime for hESC -> cardiomyocyte differentiation (cm mode)
# Requires: integrated_annotated.rds from steps 04+05 of a SCRNA_SPECIES=cm run
# Outputs (DIRS$trajectory): trajectory_cm.pdf, pseudotime_per_cell.csv,
#          pseudotime_genes.csv (Moran's I), monocle3_cds_cm.rds
# Run:  SCRNA_SPECIES=cm SCRNA_RESULTS_DIR=Results/results_H1D01_7samples_filtered \
#         Rscript pipeline/projects/cm/trajectory.R
# =============================================================================

suppressPackageStartupMessages({
  library(monocle3); library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
})
.pipeline_dir <- {
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "."
}
source(file.path(dirname(dirname(.pipeline_dir)), "config.R"))  # core config is two dirs up (pipeline/)
dir.create(DIRS$trajectory, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

# Lineage set from the first CLI arg (default "cardiac"). Excluded everywhere: Unknown (ambient
# cluster 0), Hepatic/Endoderm (other germ layer), Endothelial/Epithelial (off-lineage, too few
# cells), Proliferating (cell-cycle state spanning several lineages; short-circuits the graph).
LINEAGES <- list(
  cardiac = c("Pluripotent", "Cardiomyocyte", "Proepicardial", "Cardiac progenitor",
              "Epicardial", "Fibroblast", "Myofibroblast"),
  cm      = c("Pluripotent", "Cardiomyocyte")   # maturation axis only
)
LINEAGE <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(LINEAGE)) LINEAGE <- "cardiac"
if (!LINEAGE %in% names(LINEAGES)) stop("Unknown lineage '", LINEAGE, "'; use one of: ", paste(names(LINEAGES), collapse = ", "))
LINEAGE_TYPES <- LINEAGES[[LINEAGE]]
SUFFIX <- if (LINEAGE == "cardiac") "" else paste0("_", LINEAGE)
MAX_PER_TYPE  <- 6000L  # ponytail: downsample so Pluripotent (20k) doesn't dominate the graph; raise if memory allows
TREND_GENES   <- c("POU5F1", "MESP1", "ISL1", "NKX2-5", "GATA4", "TNNT2", "MYH6", "MYL7",
                   "MYL2", "TNNI1", "TNNI3", "NPPA", "WT1", "TCF21", "TBX18", "COL1A1", "POSTN", "ACTA2")

merged <- readRDS(file.path(DIRS$integrated, "integrated_annotated.rds"))
merged$timepoint <- sub("^.*(D[0-9]+).*$", "\\1", merged$sample)
tp_levels <- unique(merged$timepoint)[order(as.integer(sub("D", "", unique(merged$timepoint))))]
merged$timepoint <- factor(merged$timepoint, levels = tp_levels)

keep <- unlist(lapply(LINEAGE_TYPES, function(t) {
  cells <- colnames(merged)[merged$cell_type == t]
  if (length(cells) > MAX_PER_TYPE) sample(cells, MAX_PER_TYPE) else cells
}))
sub <- subset(merged, cells = keep)
rm(merged); invisible(gc())
message("Lineage subset: ", ncol(sub), " cells; ", paste(names(table(sub$cell_type)), table(sub$cell_type), sep = "=", collapse = ", "))

# Re-embed the subset on the Harmony space so off-lineage cells don't shape the layout.
sub <- RunUMAP(sub, reduction = "harmony", dims = seq_len(ncol(sub[["harmony"]])),
               reduction.name = "umap_lineage", verbose = FALSE)

# Build the CDS by hand (SeuratWrappers is not installable in this env).
counts <- LayerData(sub, assay = "RNA", layer = "counts")
cds <- new_cell_data_set(counts, cell_metadata = sub@meta.data,
                         gene_metadata = data.frame(gene_short_name = rownames(counts),
                                                    row.names = rownames(counts)))
cds <- estimate_size_factors(cds)
reducedDims(cds)[["UMAP"]] <- Embeddings(sub, "umap_lineage")
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE)  # one tree: D0 must connect to every fate

# Root = principal node most often nearest to D0 pluripotent cells (monocle3 docs recipe).
root_cells <- colnames(cds)[cds$timepoint == tp_levels[1] & cds$cell_type == "Pluripotent"]
if (length(root_cells) < 20) stop("Too few D0 Pluripotent cells to root the trajectory: ", length(root_cells))
closest <- principal_graph_aux(cds)[["UMAP"]]$pr_graph_cell_proj_closest_vertex
root_node <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
  as.integer(names(which.max(table(closest[root_cells, ]))))]
cds <- order_cells(cds, root_pr_nodes = root_node)
pt <- pseudotime(cds)
message("Root node ", root_node, "; unreachable cells (Inf pseudotime): ", sum(!is.finite(pt)))

df <- data.frame(cell = colnames(cds), sample = cds$sample, timepoint = cds$timepoint,
                 cell_type = cds$cell_type, pseudotime = pt)
write.csv(df, file.path(DIRS$trajectory, paste0("pseudotime_per_cell", SUFFIX, ".csv")), row.names = FALSE)
saveRDS(cds, file.path(DIRS$trajectory, paste0("monocle3_cds_cm", SUFFIX, ".rds")))

# --- Plots -------------------------------------------------------------------
cols <- CELLTYPE_COLORS[intersect(names(CELLTYPE_COLORS), LINEAGE_TYPES)]
p_pt <- plot_cells(cds, color_cells_by = "pseudotime", cell_size = 0.4, label_cell_groups = FALSE,
                   label_leaves = FALSE, label_branch_points = FALSE, label_roots = TRUE) +
  ggtitle("Pseudotime (root = D0 pluripotent)")
p_ct <- plot_cells(cds, color_cells_by = "cell_type", cell_size = 0.4, label_cell_groups = FALSE,
                   label_leaves = FALSE, label_branch_points = FALSE) +
  scale_colour_manual(values = cols) + ggtitle("Cell type")
p_tp <- plot_cells(cds, color_cells_by = "timepoint", cell_size = 0.4, label_cell_groups = FALSE,
                   label_leaves = FALSE, label_branch_points = FALSE, show_trajectory_graph = FALSE) +
  ggtitle("Timepoint")

fin <- df[is.finite(df$pseudotime), ]
p_box_tp <- ggplot(fin, aes(timepoint, pseudotime, fill = timepoint)) +
  geom_violin(scale = "width") + geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  facet_wrap(~ sample, nrow = 1, scales = "free_x") + theme_classic() + theme(legend.position = "none") +
  ggtitle("Pseudotime by sample — replicates should agree; pseudotime should rise with day")
p_box_ct <- ggplot(fin, aes(reorder(cell_type, pseudotime, median), pseudotime, fill = cell_type)) +
  geom_boxplot(outlier.shape = NA) + scale_fill_manual(values = cols) + coord_flip() +
  theme_classic() + theme(legend.position = "none") + labs(x = NULL, title = "Pseudotime by cell type")

genes <- intersect(TREND_GENES, rownames(sub))
expr <- as.matrix(LayerData(sub, assay = "RNA", layer = "data")[genes, df$cell])
trend <- do.call(rbind, lapply(genes, function(g)
  data.frame(gene = g, pseudotime = df$pseudotime, expr = expr[g, ], cell_type = df$cell_type)))
trend <- trend[is.finite(trend$pseudotime), ]
trend$gene <- factor(trend$gene, levels = genes)
p_trend <- ggplot(trend, aes(pseudotime, expr)) +
  geom_point(aes(colour = cell_type), size = 0.1, alpha = 0.15) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), colour = "black", se = FALSE) +
  scale_colour_manual(values = cols) + facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_classic(base_size = 9) + labs(y = "log-normalised expression", title = "Marker trends along pseudotime")

out_pdf <- file.path(DIRS$trajectory, paste0("trajectory_cm", SUFFIX, ".pdf"))
pdf(out_pdf, width = 16, height = 6)
print(p_pt | p_ct | p_tp)
print((p_box_tp | p_box_ct) + plot_layout(widths = c(2, 1)))
dev.off()
pdf(sub("\\.pdf$", "_genes.pdf", out_pdf), width = 16, height = 10)
print(p_trend)
dev.off()
message("Saved: ", out_pdf)

# --- Genes varying along the graph (Moran's I). Slow: minutes to an hour; runs last. -----
expressed <- rownames(cds)[Matrix::rowMeans(counts > 0) >= 0.01]
gt <- graph_test(cds[expressed, ], neighbor_graph = "principal_graph",
                 cores = max(1L, min(8L, PARALLEL$workers)))
gt <- gt[order(-gt$morans_I), ]
gt <- gt[!is.na(gt$q_value) & gt$q_value < 0.05, c("gene_short_name", "morans_I", "morans_test_statistic", "q_value")]
write.csv(gt, file.path(DIRS$trajectory, paste0("pseudotime_genes", SUFFIX, ".csv")), row.names = FALSE)
message("Saved pseudotime_genes.csv (", nrow(gt), " genes, q < 0.05)")
message("trajectory.R complete — ", DIRS$trajectory)
