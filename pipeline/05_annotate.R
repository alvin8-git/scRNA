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

.ref_name <- if (exists("SINGLER_REF")) SINGLER_REF else "HumanPrimaryCellAtlas"
if (.ref_name == "MonacoImmune") {
  ref        <- MonacoImmuneData()
  ref_labels <- ref$label.main
  message("  Reference: MonacoImmuneData (blood-optimised, 29 types)")
} else {
  ref        <- HumanPrimaryCellAtlasData()
  ref_labels <- ref$label.main
  message("  Reference: HumanPrimaryCellAtlasData (broad)")
}

expr_mat       <- GetAssayData(merged, assay = "RNA", layer = "data")
singler_result <- SingleR(test = expr_mat, ref = ref, labels = ref_labels,
                           fine.tune = TRUE, prune = TRUE, BPPARAM = bp_param)

merged$singler_label  <- singler_result$labels
merged$singler_pruned <- singler_result$pruned.labels
merged$singler_delta  <- apply(singler_result$scores, 1,
                                function(x) { s <- sort(x, decreasing = TRUE); s[1] - s[2] })
merged$singler_label_clean <- ifelse(is.na(merged$singler_pruned),
                                      "Unassigned", merged$singler_pruned)

# Normalise raw SingleR labels to match CELLTYPE_COLORS canonical names
SINGLER_NORM <- c(
  # --- HumanPrimaryCellAtlas labels ---
  # Platelet / megakaryocyte
  "Platelets"              = "Platelet",
  "Megakaryocyte"          = "Platelet",
  # Neutrophil / granulocyte precursors
  "Neutrophils"            = "Neutrophil",
  "Neutrophil_-G-CSF"      = "Neutrophil",
  "GMP"                    = "Neutrophil",
  "Pro-Myelocyte"          = "Neutrophil",
  "Myelocyte"              = "Neutrophil",
  # Basophil / eosinophil / mast
  "Basophils"              = "Basophil",
  "Eosinophils"            = "Eosinophil",
  "Mast_cells"             = "Mast cell",
  # RBC / erythroid
  "Erythrocyte"            = "RBC",
  "Erythroblast"           = "RBC",
  "BFU-E"                  = "RBC",
  "CFU-E"                  = "RBC",
  "MEP"                    = "RBC",
  # HSPC
  "CMP"                    = "HSPC",
  "HSC_-G-CSF"             = "HSPC",
  "HSC_CD34+"              = "HSPC",
  # Lymphoid
  "NK_cell"                = "NK",
  "T_cells"                = "CD4 T",
  "B_cell"                 = "B cell",
  "Pro-B_cell_CD34+"       = "B cell",
  "Pre-B_cell_CD34-"       = "B cell",
  # Stromal / tissue contamination
  "Endothelial_cells"      = "Endothelial",
  "Epithelial_cells"       = "Epithelial",
  "Fibroblasts"            = "Fibroblast",
  "Smooth_muscle_cells"    = "Smooth Muscle",
  # --- MonacoImmuneData labels ---
  "CD4+ T cells"           = "CD4 T",
  "CD8+ T cells"           = "CD8 T",
  "T regulatory cells"     = "Treg",
  "Vd2 gd T cells"         = "γδ T",
  "Non Vd2 gd T cells"     = "γδ T",
  "MAIT cells"             = "CD8 T",
  "NK cells"               = "NK",
  "B cells"                = "B cell",
  "Plasmablasts"           = "Plasma",
  "Classical monocytes"    = "CD14+ Mono",
  "Intermediate monocytes" = "CD14+ Mono",
  "Non-classical monocytes" = "FCGR3A+ Mono",
  "Plasmacytoid DC"        = "DC",
  "Myeloid DC"             = "DC",
  "Progenitor cells"       = "HSPC",
  "Low-density neutrophils" = "Neutrophil",
  "Low-density basophils"  = "Basophil",
  # Monaco label.main catch-alls (coarser labels returned for ambiguous cells)
  "Monocytes"              = "CD14+ Mono",
  "T cells"                = "CD4 T",
  "Dendritic cells"        = "DC",
  "B cells"                = "B cell",
  "NK"                     = "NK"
)
idx <- merged$singler_label_clean %in% names(SINGLER_NORM)
merged$singler_label_clean[idx] <- SINGLER_NORM[merged$singler_label_clean[idx]]


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
  labs(title = paste0("SingleR  -  ", .ref_name)) + theme_classic() +
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
# PART 1b: scType marker-based annotation (no cross-species reference needed)
# =============================================================================
message("\n--- Running scType (marker-based, bat overrides applied automatically) ---")

# Build positive marker gene sets from config MARKERS (bat overrides already applied)
.gs_pos <- list(
  "CD4 T"        = unique(c(MARKERS$T_pan, MARKERS$CD4_T)),
  "CD8 T"        = unique(c(MARKERS$T_pan, MARKERS$CD8_T)),
  "Treg"         = unique(c(MARKERS$T_pan, MARKERS$Treg)),
  "NK"           = MARKERS$NK,
  "B cell"       = MARKERS$B_cell,
  "Plasma"       = MARKERS$Plasma,
  "CD14+ Mono"   = MARKERS$CD14_mono,
  "FCGR3A+ Mono" = MARKERS$FCGR3A_mono,
  "Neutrophil"   = MARKERS$Neutrophil,
  "DC"           = MARKERS$DC,
  "Platelet"     = MARKERS$Platelet,
  "RBC"          = MARKERS$RBC,
  "HSPC"         = MARKERS$HSPC,
  "Eosinophil"   = MARKERS$Eosinophil,
  "Mast cell"    = MARKERS$Mast_cell
)
if (!is.null(MARKERS$gamma_delta_T))
  .gs_pos[["γδ T"]] <- MARKERS$gamma_delta_T

# Keep only genes present in the dataset
.gs_pos <- lapply(.gs_pos, function(g) intersect(g, rownames(merged)))
.gs_pos <- Filter(function(g) length(g) > 0, .gs_pos)

# Get scaled data (scType requires scaled expression)
.scale_mat <- tryCatch(
  GetAssayData(merged, assay = "RNA", layer = "scale.data"),
  error = function(e) matrix(nrow = 0, ncol = 0)
)
if (nrow(.scale_mat) == 0) {
  message("  scale.data absent — scaling HVGs for scType")
  hvg <- VariableFeatures(merged)
  if (length(hvg) == 0) hvg <- rownames(merged)  # defensive fallback only
  merged <- ScaleData(merged, features = hvg, verbose = FALSE)
  message("  scType: scaled ", length(hvg), " features (VariableFeatures)")
  .scale_mat <- GetAssayData(merged, assay = "RNA", layer = "scale.data")
}

# Marker sensitivity weight: genes shared by fewer cell types score higher
.gene_freq   <- table(unlist(.gs_pos))
.gene_weight <- setNames(1 / as.numeric(.gene_freq), names(.gene_freq))

# Per-cell scores for each cell type (weighted sum of scaled expression)
# Filter each gene set to rows actually present in scale.data (scale.data contains HVGs only)
.score_mat <- do.call(cbind, lapply(names(.gs_pos), function(ct) {
  genes <- intersect(.gs_pos[[ct]], rownames(.scale_mat))
  if (length(genes) == 0) return(rep(0, ncol(.scale_mat)))
  w <- .gene_weight[genes]
  as.numeric(colSums(sweep(.scale_mat[genes, , drop = FALSE], 1, w, "*")))
}))
colnames(.score_mat) <- names(.gs_pos)

# Assign per-cluster label: cell type with highest aggregate score wins
.cl_sctype <- do.call(rbind, lapply(levels(merged$seurat_clusters), function(cl) {
  cells     <- which(merged$seurat_clusters == cl)
  cl_scores <- colSums(.score_mat[cells, , drop = FALSE])
  data.frame(cluster = cl, sctype = names(which.max(cl_scores)),
             score = max(cl_scores), n_cells = length(cells),
             stringsAsFactors = FALSE)
}))
.cl_sctype <- .cl_sctype[order(as.integer(.cl_sctype$cluster)), ]

merged$sctype_label <- .cl_sctype$sctype[match(as.character(merged$seurat_clusters),
                                                .cl_sctype$cluster)]
message("  scType cluster assignments:")
print(.cl_sctype[, c("cluster", "sctype", "n_cells")])

# UMAP coloured by scType labels
.sctype_u    <- unique(merged$sctype_label)
.sctype_cols <- setNames(
  ifelse(.sctype_u %in% names(CELLTYPE_COLORS),
         CELLTYPE_COLORS[.sctype_u],
         scales::hue_pal()(length(.sctype_u))[seq_along(.sctype_u)]),
  .sctype_u)
p_sctype <- DimPlot(merged, group.by = "sctype_label", reduction = "umap",
                     cols = .sctype_cols, label = TRUE, label.size = 3,
                     pt.size = PLOT$pt_size, repel = TRUE) +
  labs(title = "scType  -  Marker-based Annotation") + theme_classic() +
  theme(legend.position = "right", legend.text = element_text(size = 8))
ggsave(file.path(DIRS$annotation, "sctype_labels_umap.pdf"),
       p_sctype, width = 10, height = 7, dpi = PLOT$dpi)
report_plots[["scType  -  Marker-based Labels on UMAP"]] <- set_page(p_sctype, 8.5, 7.5)

# SingleR vs scType comparison table per cluster
.comparison <- merge(
  cluster_singler[, c("seurat_clusters", "majority_singler")],
  .cl_sctype[, c("cluster", "sctype", "n_cells")],
  by.x = "seurat_clusters", by.y = "cluster"
)
.comparison$agreement <- .comparison$majority_singler == .comparison$sctype
.comparison <- .comparison[order(as.integer(as.character(.comparison$seurat_clusters))), ]
message("\n  SingleR vs scType comparison per cluster:")
print(.comparison)
write.csv(.comparison,
          file.path(DIRS$annotation, "singler_vs_sctype_comparison.csv"),
          row.names = FALSE)

# Save cluster→type mapping for cell types MonacoImmune cannot label.
# These are propagated to merged$cell_type after the SingleR-based assignment.
.monaco_blind <- c("RBC", "Platelet", "Eosinophil", "Mast cell")
.sctype_monaco_blind <- .cl_sctype[.cl_sctype$sctype %in% .monaco_blind, , drop = FALSE]
rm(.scale_mat, .score_mat, .gs_pos, .gene_freq, .gene_weight, .cl_sctype,
   .sctype_u, .sctype_cols, .comparison, .monaco_blind)

# =============================================================================
# PART 2: Canonical marker plots for manual annotation
# =============================================================================
message("\n--- Canonical marker plots ---")

markers_present <- lapply(MARKERS, function(g) g[g %in% rownames(merged)])

p_dot_manual <- suppressWarnings(
  DotPlot(merged, features = unique(unlist(markers_present)),
          group.by = "seurat_clusters", dot.scale = 8, dot.min = 0.01) +
    scale_color_viridis_c(option = "plasma") + RotatedAxis() +
    labs(title = "Canonical PBMC Markers  -  use this to fill CLUSTER_CELLTYPE_MAP",
         x = NULL, y = "Cluster")
)
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
if (!is.null(MARKERS$gamma_delta_T)) {
  plot_marker_group(MARKERS$gamma_delta_T,
                    "Feature Plot  -  γδ T Cell Markers (TRDC, TRGC1, TRGC2)",
                    "feature_gamma_delta_T.pdf")
}

# =============================================================================
# Annotation precedence (highest priority last — later steps override earlier):
#
#  [1] SingleR per-cell  ->  singler_label_clean  (baseline for all cells)
#      |
#  [2] scType per-cluster  ->  Monaco-blind types only
#      (RBC, Platelet, Eosinophil, Mast cell: MonacoImmune cannot detect these)
#      |
#  [3] CLUSTER_CELLTYPE_MAP  ->  manual cluster override (NULL = auto majority vote)
#      Unmapped clusters fall back to SingleR majority vote automatically.
#      |
#  [4] CONTAMINATION_TYPES  ->  per-cell SingleR label wins over cluster majority
#      (Neutrophil, RBC, HSPC, Platelet, Basophil... preserved as singletons)
#
#  Final label stored in: merged$cell_type -> merged$final_cell_type
#
#  To change precedence: reorder the blocks below. Each block sets merged$cell_type
#  for the cells it covers; later blocks overwrite earlier ones.
# =============================================================================
# PART 3: Assign final cell_type labels
# =============================================================================
message("\n--- Assigning cell type labels ---")

if (!is.null(CLUSTER_CELLTYPE_MAP)) {
  merged$cell_type <- unname(CLUSTER_CELLTYPE_MAP[as.character(merged$seurat_clusters)])
  # Fall back to SingleR per-cell label for any cluster not in the map
  unmapped <- is.na(merged$cell_type)
  if (any(unmapped)) {
    merged$cell_type[unmapped] <- merged$singler_label_clean[unmapped]
    n_mapped   <- sum(!unmapped)
    n_unmapped <- sum(unmapped)
    cl_unmapped <- sort(unique(merged$seurat_clusters[unmapped]))
    message("  Manual CLUSTER_CELLTYPE_MAP: ", n_mapped, " cells mapped manually; ",
            n_unmapped, " cells (clusters ", paste(cl_unmapped, collapse = ", "),
            ") fall back to SingleR labels")
  } else {
    message("  Using manual CLUSTER_CELLTYPE_MAP (all clusters mapped)")
  }
} else {
  # Auto-generate cluster-level labels from SingleR majority vote per cluster.
  # This gives coherent, cluster-level annotations without any manual config.
  # A copy-pasteable CLUSTER_CELLTYPE_MAP is printed to the log for review.
  cl_majority <- merged@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(auto_label = names(which.max(table(singler_label_clean))), .groups = "drop")
  label_map <- setNames(cl_majority$auto_label, as.character(cl_majority$seurat_clusters))

  # --- Sub-type refinement ---
  # For clusters whose majority label is a generic type (e.g. "CD4 T", "B cell",
  # "Monocyte"), score sub-type-specific markers by average expression per cluster
  # and assign the best-matching specific sub-type label.
  if (!is.null(SUBTYPE_MARKERS)) {
    generic_clusters <- names(label_map)[label_map %in% names(SUBTYPE_MARKERS)]

    if (length(generic_clusters) > 0) {
      all_st_genes <- unique(unlist(lapply(SUBTYPE_MARKERS, unlist)))
      all_st_genes <- intersect(all_st_genes, rownames(merged))

      if (length(all_st_genes) > 0) {
        avg_expr <- AverageExpression(merged, features = all_st_genes,
                                      group.by = "seurat_clusters", assay = "RNA")$RNA
        # Seurat prepends "g" to numeric cluster names — strip it for consistent indexing
        colnames(avg_expr) <- sub("^g", "", colnames(avg_expr))

        n_refined <- 0
        for (cl in generic_clusters) {
          major   <- label_map[cl]
          subtypes <- SUBTYPE_MARKERS[[major]]

          scores <- sapply(names(subtypes), function(st) {
            genes <- intersect(subtypes[[st]], rownames(avg_expr))
            if (length(genes) == 0) return(0)
            mean(avg_expr[genes, cl, drop = TRUE])
          })

          best <- names(which.max(scores))
          if (max(scores) > 0.05) {
            label_map[cl] <- best
            n_refined <- n_refined + 1
          }
        }

        if (n_refined > 0)
          message("  Sub-type refinement: ", n_refined,
                  " cluster(s) refined from generic to specific labels.")
      }
    }
  }

  merged$cell_type <- unname(label_map[as.character(merged$seurat_clusters)])

  # Override cluster-majority label for contamination/rare types — ensures scattered
  # cells (e.g. 3 Neutrophils inside a Monocyte cluster) keep their per-cell SingleR
  # identity on the UMAP instead of being silently absorbed into the majority label.
  n_overridden <- 0L
  for (ct in CONTAMINATION_TYPES) {
    ct_cells <- which(merged$singler_label_clean == ct)
    if (length(ct_cells) > 0) {
      merged$cell_type[ct_cells] <- ct
      message("  [CONTAM] Preserved ", length(ct_cells), " ", ct,
              " cell(s) using per-cell SingleR label.")
      n_overridden <- n_overridden + length(ct_cells)
    }
  }
  if (n_overridden > 0)
    message("  Total contamination/rare cells relabelled: ", n_overridden)

  message("  Auto-annotated using SingleR majority label per cluster.")
  message("  To override, copy the suggested CLUSTER_CELLTYPE_MAP below into config.R:\n")
  message("  CLUSTER_CELLTYPE_MAP <- c(")
  for (cl in sort(as.integer(names(label_map)))) {
    message("    \"", cl, "\" = \"", label_map[as.character(cl)], "\",")
  }
  message("  )")
}

# Propagate scType "RBC" label to cell_type for clusters scType identifies as RBC
# but SingleR (MonacoImmune) cannot detect (Monaco has no erythroid reference types).
# Skipped if the cluster is already set in CLUSTER_CELLTYPE_MAP.
# Propagate scType labels for cell types MonacoImmune cannot detect
# (RBC, Platelet, Eosinophil, Mast cell) to merged$cell_type.
# Only applies to clusters not already set in CLUSTER_CELLTYPE_MAP.
if (exists(".sctype_monaco_blind") && nrow(.sctype_monaco_blind) > 0) {
  .mapped <- names(CLUSTER_CELLTYPE_MAP %||% c())
  for (.i in seq_len(nrow(.sctype_monaco_blind))) {
    .ct <- .sctype_monaco_blind$sctype[.i]
    .cl <- .sctype_monaco_blind$cluster[.i]
    if (!(.cl %in% .mapped)) {
      .cells <- as.character(merged$seurat_clusters) == .cl
      merged$cell_type[.cells] <- .ct
      message("  [scType override] Cluster ", .cl, ": labelled ", sum(.cells),
              " cells as '", .ct, "' (MonacoImmune cannot detect this type)")
    }
  }
  rm(.sctype_monaco_blind, .mapped, .i, .ct, .cl, .cells)
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
# Contamination / rare cell type summary plot (per sample)
# =============================================================================
# Always include Monaco-blind types (RBC, Platelet, Eosinophil, Mast cell)
# if present in cell_type — they are labelled via scType not CONTAMINATION_TYPES,
# so they would otherwise be silently excluded from this summary.
.monaco_blind_present <- intersect(c("RBC", "Platelet", "Eosinophil", "Mast cell"),
                                    unique(merged$cell_type))
.summary_types  <- unique(c(CONTAMINATION_TYPES, .monaco_blind_present))
contam_present  <- intersect(.summary_types, unique(merged$cell_type))
rm(.monaco_blind_present, .summary_types)
if (length(contam_present) > 0) {
  contam_df <- merged@meta.data %>%
    group_by(sample) %>%
    mutate(total_cells = n()) %>%
    ungroup() %>%
    filter(cell_type %in% contam_present) %>%
    group_by(sample, cell_type) %>%
    summarise(n = n(), total_cells = first(total_cells), .groups = "drop") %>%
    mutate(pct = round(n / total_cells * 100, 2),
           label = ifelse(pct >= 1, paste0(pct, "%"), ""))

  contam_cols <- CELLTYPE_COLORS[names(CELLTYPE_COLORS) %in% contam_present]
  missing_col <- setdiff(contam_present, names(contam_cols))
  if (length(missing_col) > 0)
    contam_cols <- c(contam_cols, setNames(hue_pal()(length(missing_col)), missing_col))

  .n_per_row   <- 8L
  .n_per_page  <- 16L   # 2 rows × 8 samples per page
  .page_chunks <- split(SAMPLE_NAMES, ceiling(seq_along(SAMPLE_NAMES) / .n_per_page))
  .n_pages     <- length(.page_chunks)
  .contam_pdfs <- character(0)

  for (.pg in seq_along(.page_chunks)) {
    .samps <- .page_chunks[[.pg]]
    .pg_df <- contam_df %>%
      filter(sample %in% .samps) %>%
      mutate(
        sample    = factor(sample, levels = .samps),
        row_group = factor(
          ifelse(as.integer(factor(sample, levels = .samps)) <= .n_per_row,
                 "Row 1", "Row 2"),
          levels = c("Row 1", "Row 2")
        )
      )
    if (nrow(.pg_df) == 0) next

    .suffix <- if (.n_pages > 1) paste0(" (", .pg, "/", .n_pages, ")") else ""
    .n_rg   <- length(unique(.pg_df$row_group))
    .pg_w   <- max(5, min(.n_per_row, length(.samps)) * 1.8)
    .pg_h   <- 3.5 * .n_rg

    .p <- ggplot(.pg_df,
                 aes(x = sample, y = pct, fill = cell_type, label = label)) +
      geom_col(position = "stack", width = 0.6) +
      geom_text(position = position_stack(vjust = 0.5), size = 3, color = "white",
                fontface = "bold") +
      scale_fill_manual(values = contam_cols) +
      scale_x_discrete(drop = TRUE) +
      facet_wrap(~row_group, nrow = 2, scales = "free_x") +
      labs(title = paste0("Contamination & Rare Cell Types — % per Sample", .suffix),
           x = NULL, y = "% of cells", fill = "Cell Type") +
      theme_classic() +
      theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 10),
            strip.text       = element_blank(),
            strip.background = element_blank())

    .tmp <- tempfile(fileext = ".pdf")
    ggsave(.tmp, .p, width = .pg_w, height = .pg_h, dpi = PLOT$dpi, limitsize = FALSE)
    .contam_pdfs <- c(.contam_pdfs, .tmp)

    .key <- if (.pg == 1) "Contamination & Rare Cell Types" else
      paste0("Contamination & Rare Cell Types (", .pg, "/", .n_pages, ")")
    report_plots[[.key]] <- set_page(.p, pw = 8, ph = min(.pg_h * 8 / .pg_w, 9))
  }

  .combine_pdfs(.contam_pdfs,
                file.path(DIRS$annotation, "contamination_summary.pdf"))
  rm(.n_per_row, .n_per_page, .page_chunks, .n_pages, .contam_pdfs,
     .pg, .samps, .pg_df, .suffix, .n_rg, .pg_w, .pg_h, .p, .tmp, .key)

  message("\n  Contamination/rare cell types detected:")
  print(as.data.frame(select(contam_df, sample, cell_type, n, pct)))
} else {
  message("\n  No contamination/rare cell types detected — samples appear clean.")
}

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
      p_sub_dot <- suppressWarnings(
        DotPlot(merged_t, features = t_markers,
                dot.scale = 8, dot.min = 0.01) +
          scale_color_viridis_c(option = "plasma") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)) +
          labs(title = "T Cell Sub-clusters  -  CD4 / CD8 / Treg Markers",
               x = NULL, y = "Sub-cluster")
      )
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
