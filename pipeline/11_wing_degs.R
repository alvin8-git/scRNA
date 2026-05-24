# =============================================================================
# 11_wing_degs.R — Condition-aware DEG and wound module scores (bat_wing mode)
# Requires: integrated_annotated.rds from step 04+05
# Outputs:  DEG CSVs per cell type, volcano plots, module score UMAPs/violins,
#           top DEG heatmap, all written to DIRS$differential
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
  library(ggrepel); library(pdftools); library(magick)
})

.pipeline_dir <- {
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "."
}
source(file.path(.pipeline_dir, "config.R"))

if (length(CONDITION_LEVELS) < 2) {
  message("Only one condition found — skipping 11_wing_degs.R (set SCRNA_CONDITION)")
  quit(save = "no", status = 0)
}

dir.create(DIRS$differential, showWarnings = FALSE, recursive = TRUE)

.combine_pdfs <- function(paths, out) {
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) return(invisible(NULL))
  imgs <- lapply(paths, function(p) {
    n <- tryCatch(pdftools::pdf_length(p), error = function(e) 1L)
    lapply(seq_len(n), function(i) magick::image_read_pdf(p, pages = i, density = 150))
  })
  all_imgs <- unlist(imgs, recursive = FALSE)
  magick::image_write(magick::image_join(all_imgs), out, format = "pdf")
}

# --- Load annotated object ---------------------------------------------------
rds_path <- file.path(DIRS$integrated, "integrated_annotated.rds")
if (!file.exists(rds_path))
  stop("Missing: ", rds_path, " — run steps 04+05 first")
merged <- readRDS(rds_path)

merged$condition <- unname(SAMPLE_CONDITIONS[merged$sample])
if (anyNA(merged$condition))
  warning("Some cells have no condition label — check SCRNA_CONDITION")

cond_a <- CONDITION_LEVELS[1]
cond_b <- CONDITION_LEVELS[2]
message("Comparing: ", cond_b, " vs ", cond_a, " (positive FC = up in ", cond_b, ")")

Idents(merged) <- "cell_type"
cell_types <- sort(unique(merged$cell_type))
cell_types  <- cell_types[!cell_types %in% CONTAMINATION_TYPES]

# =============================================================================
# PART 1: Module scores — wound healing gene sets
# =============================================================================
if (length(WOUND_MODULES) > 0) {
  message("Computing wound module scores...")
  for (mod in names(WOUND_MODULES)) {
    genes <- WOUND_MODULES[[mod]]
    genes <- genes[genes %in% rownames(merged)]
    if (length(genes) < 2) {
      message("  Skipping module '", mod, "' — fewer than 2 genes present")
      next
    }
    merged <- AddModuleScore(merged, features = list(genes),
                             name = paste0(mod, "_score"), seed = 42)
    # AddModuleScore appends "1" to the name
    old_col <- paste0(mod, "_score1")
    new_col <- paste0(mod, "_score")
    if (old_col %in% colnames(merged@meta.data))
      colnames(merged@meta.data)[colnames(merged@meta.data) == old_col] <- new_col
  }
}

score_cols <- paste0(names(WOUND_MODULES), "_score")
score_cols <- score_cols[score_cols %in% colnames(merged@meta.data)]

if (length(score_cols) > 0) {
  # UMAP coloured by module score, split by condition
  mod_pdf <- file.path(DIRS$differential, "wound_module_scores_umap.pdf")
  pdf(mod_pdf, width = 10, height = 4)
  for (sc in score_cols) {
    tryCatch({
      p <- FeaturePlot(merged, features = sc, split.by = "condition",
                       pt.size = 0.3, order = TRUE) &
        scale_colour_gradientn(colours = c("grey90", "#F4A460", "#8B0000")) &
        theme_classic(base_size = 9)
      print(p)
    }, error = function(e) message("  FeaturePlot failed for ", sc, ": ", e$message))
  }
  dev.off()

  # Violin per cell type, split by condition
  vln_pdf <- file.path(DIRS$differential, "wound_module_scores_violin.pdf")
  pdf(vln_pdf, width = 12, height = 4)
  for (sc in score_cols) {
    tryCatch({
      p <- VlnPlot(merged, features = sc, group.by = "cell_type",
                   split.by = "condition", pt.size = 0,
                   cols = c("steelblue", "tomato")) +
        theme_classic(base_size = 8) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = gsub("_score", "", sc))
      print(p)
    }, error = function(e) message("  VlnPlot failed for ", sc, ": ", e$message))
  }
  dev.off()
  message("Saved module score plots.")
}

# =============================================================================
# PART 2: FindMarkers DEG per cell type (MAST, recovering vs healthy)
# =============================================================================
message("Running FindMarkers per cell type (MAST)...")
de_list <- list()

for (ct in cell_types) {
  cells_a <- colnames(merged)[merged$cell_type == ct & merged$condition == cond_a]
  cells_b <- colnames(merged)[merged$cell_type == ct & merged$condition == cond_b]
  if (length(cells_a) < 10 || length(cells_b) < 10) {
    message("  Skipping ", ct, " — too few cells (", length(cells_a),
            " vs ", length(cells_b), ")")
    next
  }
  message("  ", ct, ": n=", length(cells_b), " vs n=", length(cells_a))
  tryCatch({
    de <- FindMarkers(merged,
                      ident.1 = cells_b, ident.2 = cells_a,
                      test.use = "MAST", logfc.threshold = 0.1,
                      min.pct = 0.05, verbose = FALSE)
    de$gene      <- rownames(de)
    de$cell_type <- ct
    de$contrast  <- paste0(cond_b, "_vs_", cond_a)
    de_list[[ct]] <- de
  }, error = function(e) message("  MAST failed for ", ct, ": ", e$message))
}

# Save per-cell-type CSVs
de_dir <- file.path(DIRS$differential, "DEG_per_celltype")
dir.create(de_dir, showWarnings = FALSE)
for (ct in names(de_list)) {
  out_csv <- file.path(de_dir, paste0(gsub("[/ ]", "_", ct), "_DEG.csv"))
  write.csv(de_list[[ct]], out_csv, row.names = FALSE)
}

all_de <- if (length(de_list) > 0) dplyr::bind_rows(de_list) else data.frame()
write.csv(all_de, file.path(DIRS$differential, "all_DEG_combined.csv"), row.names = FALSE)
message("Saved DEG tables. Total sig genes: ",
        if (nrow(all_de) > 0) sum(all_de$p_val_adj < 0.05, na.rm = TRUE) else 0)

# Pseudo-bulk warning
n_per_cond <- table(
  unique(merged@meta.data[, c("sample", "condition")])$condition
)
if (any(n_per_cond < 2)) {
  message("WARNING: n=1 sample per condition detected.")
  message("  FindMarkers treats cells as independent — p-values are anti-conservative.")
  message("  Results are exploratory. Add biological replicates for pseudo-bulk DESeq2.")
}

# =============================================================================
# PART 3: Volcano plots per cell type
# =============================================================================
.volcano <- function(de, ct, cb, ca, top_n = 15) {
  de <- de %>% dplyr::mutate(
    sig = p_val_adj < 0.05 & abs(avg_log2FC) > 0.5,
    direction = dplyr::case_when(
      sig & avg_log2FC > 0  ~ paste0("Up in ", cb),
      sig & avg_log2FC < 0  ~ paste0("Up in ", ca),
      TRUE                  ~ "NS"
    ),
    label = ifelse(sig & dplyr::dense_rank(-abs(avg_log2FC)) <= top_n,
                   gene, NA_character_)
  )
  pal <- setNames(
    c("tomato", "steelblue", "grey70"),
    c(paste0("Up in ", cb), paste0("Up in ", ca), "NS")
  )
  ggplot(de, aes(avg_log2FC, -log10(p_val_adj + 1e-300), colour = direction)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                    na.rm = TRUE) +
    scale_colour_manual(values = pal) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    labs(title = ct, x = paste0("log2 FC (", cb, " / ", ca, ")"),
         y = "-log10(adj p)", colour = NULL) +
    theme_classic(base_size = 10)
}

vol_pdfs <- character(0)
for (ct in names(de_list)) {
  p   <- .volcano(de_list[[ct]], ct, cond_b, cond_a)
  tmp <- tempfile(fileext = ".pdf")
  ggsave(tmp, p, width = 6, height = 5, dpi = 150)
  vol_pdfs <- c(vol_pdfs, tmp)
}
if (length(vol_pdfs) > 0) {
  .combine_pdfs(vol_pdfs, file.path(DIRS$differential, "volcano_plots.pdf"))
  unlink(vol_pdfs)
}

# =============================================================================
# PART 4: Top DEG heatmap across cell types
# =============================================================================
if (nrow(all_de) > 0) {
  top_genes <- all_de %>%
    dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
    dplyr::group_by(cell_type) %>%
    dplyr::slice_max(abs(avg_log2FC), n = 5) %>%
    dplyr::pull(gene) %>% unique()

  if (length(top_genes) > 2) {
    top_genes <- top_genes[top_genes %in% rownames(merged)]
    avg <- AverageExpression(merged, features = top_genes,
                             group.by = c("cell_type", "condition"),
                             slot = "data")$RNA
    hm_pdf <- file.path(DIRS$differential, "top_DEG_heatmap.pdf")
    pdf(hm_pdf,
        width  = max(8, ncol(avg) * 0.5),
        height = max(6, length(top_genes) * 0.25))
    pheatmap::pheatmap(t(scale(t(as.matrix(avg)))),
                       fontsize = 7, angle_col = 45,
                       color = colorRampPalette(
                         c("navy", "white", "firebrick3"))(100),
                       main = paste0("Top DEGs: ", cond_b, " vs ", cond_a))
    dev.off()
    message("Saved top DEG heatmap.")
  }
}

message("Step 11 complete — DEG results in: ", DIRS$differential)
