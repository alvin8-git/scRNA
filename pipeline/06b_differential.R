# =============================================================================
# 06b_differential.R - Differential expression between samples per cell type.
#   Outputs: per-cell-type DE CSVs, volcano plots, inflammatory module scores,
#            differential_report.pdf
# =============================================================================
.pipeline_dir <- local({
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(f)) dirname(f)
  else {
    a <- commandArgs(trailingOnly = FALSE)
    d <- sub("--file=", "", a[grep("--file=", a)])
    if (length(d) > 0) dirname(normalizePath(d)) else "."
  }
})
source(file.path(.pipeline_dir, "config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

plan("multicore", workers = PARALLEL$merge_workers)
options(future.globals.maxSize = PARALLEL$merge_mem_gb * 1024^3)

# --- Skip if single sample ---------------------------------------------------
if (SINGLE_SAMPLE) {
  message("Single-sample run — differential expression skipped.")
  quit(status = 0)
}

merged <- readRDS(file.path(DIRS$integrated, "integrated_annotated.rds"))
Idents(merged) <- "cell_type"

sample_col <- if ("orig.ident" %in% colnames(merged@meta.data)) "orig.ident" else "sample"
samples    <- sort(unique(merged@meta.data[[sample_col]]))

if (length(samples) != 2) {
  message("DE currently supports exactly 2 samples. Found: ", paste(samples, collapse = ", "))
  message("Skipping DE.")
  quit(status = 0)
}
s1 <- samples[1]; s2 <- samples[2]
message("Loaded: ", ncol(merged), " cells | ", length(unique(merged$cell_type)),
        " cell types | Comparing: ", s1, " vs ", s2)

theme_pub <- theme_classic(base_size = 11) +
  theme(plot.title       = element_text(size = 12, face = "bold"),
        plot.subtitle    = element_text(size = 9, color = "grey40"),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 9))

# =============================================================================
# PART 1: Gene set module scores
# =============================================================================
message("\n--- Computing inflammatory module scores ---")

GENE_SETS <- list(
  Inflammatory = c("TNF", "IL1B", "IL6", "CXCL8", "NFKB1", "CCL2", "CCL3",
                   "CCL4", "S100A8", "S100A9", "IL18", "PTGS2"),
  Interferon   = c("IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "OAS1", "OAS2",
                   "IRF7", "CXCL10", "STAT1"),
  T_Activation = c("IFNG", "GZMB", "PRF1", "TNFRSF4", "CD69", "IL2RA", "HAVCR2")
)
GENE_SETS <- lapply(GENE_SETS, function(g) intersect(g, rownames(merged)))
GENE_SETS <- GENE_SETS[sapply(GENE_SETS, length) >= 3]

score_cols <- c()
for (nm in names(GENE_SETS)) {
  merged <- AddModuleScore(merged, features = list(GENE_SETS[[nm]]),
                            name = paste0("__", nm))
  real_col <- paste0(nm, "_score")
  merged@meta.data[[real_col]] <- merged@meta.data[[paste0("__", nm, "1")]]
  merged@meta.data[[paste0("__", nm, "1")]] <- NULL
  score_cols <- c(score_cols, real_col)
  message("  ", nm, ": ", length(GENE_SETS[[nm]]), " genes scored")
}

# =============================================================================
# PART 2: Differential expression per cell type
# =============================================================================
message("\n--- Differential expression per cell type ---")

MIN_CELLS  <- 20
cell_types <- sort(unique(merged$cell_type))
de_results <- list()
report_plots <- list()

for (ct in cell_types) {
  sub  <- subset(merged, cell_type == ct)
  n_s1 <- sum(sub@meta.data[[sample_col]] == s1)
  n_s2 <- sum(sub@meta.data[[sample_col]] == s2)

  if (n_s1 < MIN_CELLS || n_s2 < MIN_CELLS) {
    message("  Skipping '", ct, "': ", n_s1, " ", s1, " / ", n_s2, " ", s2,
            " (need >= ", MIN_CELLS, " per group)")
    next
  }
  message("  '", ct, "': ", n_s1, " ", s1, " vs ", n_s2, " ", s2)

  Idents(sub) <- sample_col

  tryCatch({
    markers <- FindMarkers(sub, ident.1 = s1, ident.2 = s2,
                           min.pct = 0.1, logfc.threshold = 0.1,
                           test.use = "wilcox", verbose = FALSE)
    markers$gene      <- rownames(markers)
    markers$cell_type <- ct
    markers$comparison <- paste0(s1, "_vs_", s2)

    ct_safe  <- gsub("[^A-Za-z0-9]", "_", ct)
    write.csv(markers, file.path(DIRS$differential, paste0("DE_", ct_safe, ".csv")),
              row.names = FALSE)
    de_results[[ct]] <- markers

    # --- Volcano plot ---
    markers$direction <- ifelse(
      markers$p_val_adj < 0.05 & markers$avg_log2FC >  0.5, s1,
      ifelse(markers$p_val_adj < 0.05 & markers$avg_log2FC < -0.5, s2, "NS"))

    # Top 12 genes to label (highest |FC| among significant)
    top_genes <- markers %>%
      filter(direction != "NS") %>%
      arrange(desc(abs(avg_log2FC))) %>%
      slice_head(n = 12) %>%
      pull(gene)
    markers$label <- ifelse(markers$gene %in% top_genes, markers$gene, "")

    n_up <- sum(markers$direction == s1)
    n_dn <- sum(markers$direction == s2)

    sample_pal <- setNames(
      c(unname(SAMPLE_COLORS[s1]), unname(SAMPLE_COLORS[s2]), "#CCCCCC"),
      c(s1, s2, "NS"))
    if (is.na(sample_pal[s1])) sample_pal[s1] <- "#E64B35"
    if (is.na(sample_pal[s2])) sample_pal[s2] <- "#4DBBD5"

    p_vol <- ggplot(markers,
                    aes(x = avg_log2FC,
                        y = -log10(p_val_adj + 1e-300),
                        color = direction, label = label)) +
      geom_point(alpha = 0.6, size = 1.5) +
      geom_text(size = 2.4, hjust = -0.1, check_overlap = TRUE) +
      scale_color_manual(values = sample_pal) +
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.4) +
      geom_hline(yintercept = -log10(0.05),  linetype = "dashed", alpha = 0.4) +
      labs(title    = paste0(ct, "  —  ", s1, " vs ", s2),
           subtitle = paste0("Up in ", s1, ": ", n_up,
                             "   |   Up in ", s2, ": ", n_dn,
                             "   |   NS: ", sum(markers$direction == "NS")),
           x = paste0("log2FC  (", s1, " / ", s2, ")"),
           y = "-log10(adj. p-value)", color = "Higher in") +
      theme_pub

    ggsave(file.path(DIRS$differential, paste0("volcano_", ct_safe, ".pdf")),
           p_vol, width = 7, height = 5, dpi = PLOT$dpi)
    report_plots[[paste0("Volcano — ", ct)]] <- p_vol

  }, error = function(e) {
    message("  WARNING: DE failed for '", ct, "': ", e$message)
  })
}

# --- Combined DE table -------------------------------------------------------
if (length(de_results) > 0) {
  all_de <- do.call(rbind, de_results)
  write.csv(all_de,
            file.path(DIRS$differential, "DE_all_celltypes.csv"),
            row.names = FALSE)

  summary_df <- do.call(rbind, lapply(names(de_results), function(ct) {
    df <- de_results[[ct]]
    data.frame(
      cell_type  = ct,
      n_cells_s1 = sum(merged@meta.data[[sample_col]] == s1 & merged$cell_type == ct),
      n_cells_s2 = sum(merged@meta.data[[sample_col]] == s2 & merged$cell_type == ct),
      n_DE_total = sum(df$p_val_adj < 0.05 & abs(df$avg_log2FC) > 0.5, na.rm = TRUE),
      n_up_s1    = sum(df$p_val_adj < 0.05 & df$avg_log2FC >  0.5, na.rm = TRUE),
      n_up_s2    = sum(df$p_val_adj < 0.05 & df$avg_log2FC < -0.5, na.rm = TRUE),
      top3_up_s1 = paste(head(df$gene[!is.na(df$p_val_adj) & df$p_val_adj < 0.05 & df$avg_log2FC > 0], 3), collapse = ", "),
      top3_up_s2 = paste(head(df$gene[!is.na(df$p_val_adj) & df$p_val_adj < 0.05 & df$avg_log2FC < 0], 3), collapse = ", ")
    )
  }))
  write.csv(summary_df,
            file.path(DIRS$differential, "DE_summary.csv"),
            row.names = FALSE)
  message("\nDE Summary:")
  print(summary_df, row.names = FALSE)
}

# =============================================================================
# PART 3: Module score plots per cell type
# =============================================================================
message("\n--- Module score plots ---")

meta <- merged@meta.data
meta[[sample_col]] <- factor(meta[[sample_col]], levels = c(s1, s2))

for (nm in names(GENE_SETS)) {
  score_col <- paste0(nm, "_score")
  if (!score_col %in% colnames(meta)) next

  p_mod <- ggplot(meta,
                  aes(x = .data[[sample_col]],
                      y = .data[[score_col]],
                      fill = .data[[sample_col]])) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.4) +
    facet_wrap(~ cell_type, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c(
      setNames(unname(SAMPLE_COLORS[s1]), s1),
      setNames(unname(SAMPLE_COLORS[s2]), s2))) +
    labs(title    = paste0(gsub("_", " ", nm), " Signature Score"),
         subtitle = paste0(s1, " vs ", s2, "  |  genes: ",
                           paste(GENE_SETS[[nm]], collapse = ", ")),
         x = NULL, y = "Module Score") +
    theme_pub +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(DIRS$differential, paste0("module_score_", nm, ".pdf")),
         p_mod, width = 11, height = 8, dpi = PLOT$dpi)
  report_plots[[paste0("Module Score — ", nm)]] <- p_mod
  message("  Saved: module_score_", nm, ".pdf")
}

# =============================================================================
# PART 4: Save combined report PDF
# =============================================================================
if (length(report_plots) > 0) {
  pdf(file.path(DIRS$differential, "differential_report.pdf"), width = 8, height = 6)
  for (p in report_plots) print(p)
  dev.off()
  message("\nReport saved: ", file.path(DIRS$differential, "differential_report.pdf"))
}

message("\n06b_differential.R complete. Outputs in: ", DIRS$differential)
