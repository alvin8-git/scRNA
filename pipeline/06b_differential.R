# =============================================================================
# 06b_differential.R - Differential expression between CONDITIONS per cell type,
#   on sample-level pseudobulk counts (edgeR quasi-likelihood).
#
#   Why pseudobulk: per-cell tests treat every cell as an independent replicate and
#   return p-values that are anti-conservative by orders of magnitude. The unit of
#   replication in a multi-sample design is the sample. Counts are summed per
#   sample x cell type, filtered with filterByExpr, TMM-normalised and fitted with
#   ~ condition. Each condition needs >= 2 samples; otherwise the step writes a
#   DESCRIPTIVE log2FC table (no p-values) and says so.
#
#   Conditions come from SCRNA_CONDITION (config.R -> SAMPLE_CONDITIONS). With no
#   conditions set and exactly 2 samples, each sample is its own condition and the
#   output is descriptive. Contrasts: every pair of conditions (later vs earlier in
#   SCRNA_CONDITION order), so a 2x2 design gets its within-timepoint contrasts.
#
#   Outputs (differential/): DE_<celltype>.csv, DE_all_celltypes.csv, DE_summary.csv,
#            pseudobulk_counts.csv.gz, volcano_<celltype>.pdf, module_score_*.pdf,
#            differential_report.pdf
#   Columns kept compatible with 08b_html_report.R: gene, cell_type, comparison,
#            avg_log2FC, p_val, p_val_adj (NA when descriptive).
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
  library(edgeR)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

if (SINGLE_SAMPLE) {
  message("Single-sample run — differential expression skipped.")
  quit(status = 0)
}

merged <- readRDS(file.path(DIRS$integrated, "integrated_annotated.rds"))
merged <- apply_reference_labels(merged)
Idents(merged) <- "cell_type"

sample_col <- if ("sample" %in% colnames(merged@meta.data)) "sample" else "orig.ident"
samples    <- sort(unique(merged@meta.data[[sample_col]]))

# ---- design -----------------------------------------------------------------
cond <- SAMPLE_CONDITIONS[samples]
if (length(unique(na.omit(cond))) < 2) {
  if (length(samples) == 2) {
    cond <- setNames(samples, samples)        # legacy 2-sample run: sample = condition
    message("No SCRNA_CONDITION set — treating each of the 2 samples as its own condition (descriptive only).")
  } else {
    message("DE needs >= 2 conditions. Set SCRNA_CONDITION=sample=label,... (",
            length(samples), " samples, 1 condition). Skipping DE.")
    quit(status = 0)
  }
}
cond[is.na(cond)] <- "unlabelled"
cond_levels <- unique(unname(cond))
n_per_cond  <- table(factor(cond, levels = cond_levels))
message("Loaded: ", ncol(merged), " cells | ", length(unique(merged$cell_type)), " cell types | ",
        "conditions: ", paste(names(n_per_cond), n_per_cond, sep = "=", collapse = ", "),
        if (all(n_per_cond >= 2)) " | edgeR QL, all contrasts tested"
        else paste0(" | contrasts involving ", paste(names(n_per_cond)[n_per_cond < 2], collapse = "/"),
                    " are DESCRIPTIVE (single sample)"))

theme_pub <- theme_classic(base_size = 11) +
  theme(plot.title       = element_text(size = 12, face = "bold"),
        plot.subtitle    = element_text(size = 9, color = "grey40"),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 9))

cond_pal <- setNames(c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4")[seq_along(cond_levels)],
                     cond_levels)

# =============================================================================
# PART 1: Gene set module scores (unchanged: per-cell scores, shown per condition)
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
for (nm in names(GENE_SETS)) {
  merged <- AddModuleScore(merged, features = list(GENE_SETS[[nm]]), name = paste0("__", nm))
  merged@meta.data[[paste0(nm, "_score")]] <- merged@meta.data[[paste0("__", nm, "1")]]
  merged@meta.data[[paste0("__", nm, "1")]] <- NULL
  message("  ", nm, ": ", length(GENE_SETS[[nm]]), " genes scored")
}
merged$condition <- unname(cond[merged@meta.data[[sample_col]]])

# =============================================================================
# PART 2: Pseudobulk per sample x cell type
# =============================================================================
message("\n--- Pseudobulk (sample x cell type) ---")
MIN_CELLS <- 20     # a pseudobulk column needs this many cells to count as a replicate
# Summed directly from the count layer: AggregateExpression rewrites group names
# ("_" -> "-", a leading "g" when a value starts with a digit), which makes mapping
# columns back to sample/cell type fragile. A sparse indicator product is exact.
.cnt <- SeuratObject::LayerData(merged, assay = "RNA", layer = "counts")
grp  <- factor(paste(merged@meta.data[[sample_col]], merged$cell_type, sep = "|"))
pb   <- as.matrix(.cnt %*% Matrix::t(Matrix::fac2sparse(grp)))   # genes x (sample|cell_type)
colnames(pb) <- levels(grp)
key  <- data.frame(col = levels(grp), stringsAsFactors = FALSE)
key[[sample_col]] <- sub("\\|.*$", "", key$col)
key$cell_type     <- sub("^[^|]*\\|", "", key$col)
key$n_cells       <- as.integer(table(grp)[key$col])
key$condition     <- unname(cond[key[[sample_col]]])
rm(.cnt)

write.csv(data.frame(gene = rownames(pb), pb, check.names = FALSE),
          gzfile(file.path(DIRS$differential, "pseudobulk_counts.csv.gz")), row.names = FALSE)
write.csv(key[, c("col", sample_col, "cell_type", "condition", "n_cells")],
          file.path(DIRS$differential, "pseudobulk_samples.csv"), row.names = FALSE)

# =============================================================================
# PART 3: edgeR QL per cell type, each condition vs the first
# =============================================================================
message("\n--- Differential expression per cell type ---")
de_results <- list(); report_plots <- list()

for (ct in sort(unique(key$cell_type))) {
  k <- key[key$cell_type == ct & key$n_cells >= MIN_CELLS, ]
  have <- table(factor(k$condition, levels = cond_levels))
  if (sum(have >= 1) < 2 || have[1] < 1) {
    message("  Skipping '", ct, "': ", paste(names(have), have, sep = "=", collapse = ", "),
            " samples with >= ", MIN_CELLS, " cells")
    next
  }
  # Dispersion can be estimated whenever the design leaves residual df (samples >
  # conditions present). Each CONTRAST then gets p-values only if both of its arms
  # have >= 2 samples; a 1-vs-n contrast stays descriptive while the rest are tested.
  ct_inferential <- (nrow(k) - sum(have > 0)) > 0
  y <- DGEList(pb[, k$col, drop = FALSE], group = factor(k$condition, levels = cond_levels))
  keep <- filterByExpr(y, group = y$samples$group)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~ 0 + group, data = y$samples)
  colnames(design) <- sub("^group", "", colnames(design))
  present <- names(have)[have > 0]
  pairs   <- combn(present, 2, simplify = FALSE)            # every pair: later vs earlier level
  contr   <- lapply(pairs, function(p) { v <- setNames(numeric(ncol(design)), colnames(design)); v[p[2]] <- 1; v[p[1]] <- -1; v })
  message("  '", ct, "': ", paste(names(have), have, sep = "=", collapse = ", "),
          " | ", sum(keep), " genes", if (!ct_inferential) " | descriptive (no replicates)")

  res_ct <- tryCatch({
    if (ct_inferential) {
      y   <- estimateDisp(y, design, robust = TRUE)
      fit <- glmQLFit(y, design, robust = TRUE)
    } else {
      # No replicates in some arm: fixed BCV 0.4 gives a sensible log2FC shrinkage but
      # NO valid p-value, so p columns are set to NA below.
      fit <- glmFit(y, design, dispersion = 0.4^2)
    }
    bind_rows(lapply(seq_along(pairs), function(i) {
      ref <- pairs[[i]][1]; alt <- pairs[[i]][2]
      inf_i <- ct_inferential && have[ref] >= 2 && have[alt] >= 2
      tt <- if (ct_inferential) topTags(glmQLFTest(fit, contrast = contr[[i]]), n = Inf)$table
            else topTags(glmLRT(fit, contrast = contr[[i]]), n = Inf)$table   # rows sorted best-first
      data.frame(gene = rownames(tt), cell_type = ct,
                 comparison = paste0(alt, "_vs_", ref),
                 avg_log2FC = tt$logFC, logCPM = tt$logCPM,
                 p_val = if (inf_i) tt$PValue else NA_real_,
                 p_val_adj = if (inf_i) tt$FDR else NA_real_,
                 n_samples_ref = unname(have[ref]), n_samples_alt = unname(have[alt]),
                 inferential = inf_i, stringsAsFactors = FALSE)
    }))
  }, error = function(e) { message("  WARNING: DE failed for '", ct, "': ", e$message); NULL })
  if (is.null(res_ct) || !nrow(res_ct)) next

  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  write.csv(res_ct, file.path(DIRS$differential, paste0("DE_", ct_safe, ".csv")), row.names = FALSE)
  de_results[[ct]] <- res_ct

  # ---- volcano (or MA when descriptive) per contrast ----
  for (cmp in unique(res_ct$comparison)) {
    d <- res_ct[res_ct$comparison == cmp, ]
    alt <- sub("_vs_.*", "", cmp); ref <- sub(".*_vs_", "", cmp)
    if (d$inferential[1]) {
      d$direction <- ifelse(d$p_val_adj < 0.05 & d$avg_log2FC >  1, alt,
                     ifelse(d$p_val_adj < 0.05 & d$avg_log2FC < -1, ref, "NS"))
      top <- d %>% filter(direction != "NS") %>% arrange(p_val_adj) %>% slice_head(n = 12) %>% pull(gene)
      d$label <- ifelse(d$gene %in% top, d$gene, "")
      p <- ggplot(d, aes(avg_log2FC, -log10(pmax(p_val_adj, 1e-300)), colour = direction, label = label)) +
        geom_point(alpha = 0.6, size = 1.5) +
        geom_text(size = 2.4, hjust = -0.1, check_overlap = TRUE) +
        scale_colour_manual(values = c(setNames(unname(cond_pal[c(alt, ref)]), c(alt, ref)), NS = "#CCCCCC")) +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.4) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.4) +
        labs(title = paste0(ct, "  —  ", alt, " vs ", ref, "  (pseudobulk edgeR QL)"),
             subtitle = sprintf("Up in %s: %d   |   Up in %s: %d   |   samples %d vs %d",
                                alt, sum(d$direction == alt), ref, sum(d$direction == ref),
                                d$n_samples_alt[1], d$n_samples_ref[1]),
             x = paste0("log2FC (", alt, " / ", ref, ")"), y = "-log10(FDR)", colour = "Higher in") +
        theme_pub
    } else {
      p <- ggplot(d, aes(logCPM, avg_log2FC)) +
        geom_point(alpha = 0.4, size = 1.2, colour = "grey40") +
        geom_hline(yintercept = c(-1, 0, 1), linetype = c("dashed", "solid", "dashed"), alpha = 0.4) +
        labs(title = paste0(ct, "  —  ", alt, " vs ", ref, "  (DESCRIPTIVE: no replicates)"),
             subtitle = "No p-values: a condition has a single sample. log2FC only, BCV fixed at 0.4.",
             x = "Average log2 CPM", y = paste0("log2FC (", alt, " / ", ref, ")")) +
        theme_pub
    }
    ggsave(file.path(DIRS$differential, paste0("volcano_", ct_safe, "_", gsub("[^A-Za-z0-9]", "_", cmp), ".pdf")),
           p, width = 7, height = 5, dpi = PLOT$dpi)
    report_plots[[paste0("DE — ", ct, " — ", cmp)]] <- p
  }
}

# ---- combined tables ----
if (length(de_results)) {
  all_de <- bind_rows(de_results)
  write.csv(all_de, file.path(DIRS$differential, "DE_all_celltypes.csv"), row.names = FALSE)
  summary_df <- all_de %>% group_by(cell_type, comparison) %>%
    summarise(inferential = first(inferential), n_samples_ref = first(n_samples_ref),
              n_samples_alt = first(n_samples_alt), n_genes_tested = n(),
              n_DE_fdr05 = sum(p_val_adj < 0.05 & abs(avg_log2FC) > 1, na.rm = TRUE),
              n_up_alt = sum(p_val_adj < 0.05 & avg_log2FC >  1, na.rm = TRUE),
              n_up_ref = sum(p_val_adj < 0.05 & avg_log2FC < -1, na.rm = TRUE),
              top3_up_alt = paste(head(gene[avg_log2FC > 0], 3), collapse = ", "),   # rows are in topTags order
              top3_up_ref = paste(head(gene[avg_log2FC < 0], 3), collapse = ", "),
              .groups = "drop")
  write.csv(summary_df, file.path(DIRS$differential, "DE_summary.csv"), row.names = FALSE)
  message("\nDE Summary:"); print(as.data.frame(summary_df), row.names = FALSE)
} else message("\nNo cell type had enough samples for DE.")

# =============================================================================
# PART 4: Module score plots per condition
# =============================================================================
message("\n--- Module score plots ---")
meta <- merged@meta.data
meta$condition <- factor(meta$condition, levels = cond_levels)
for (nm in names(GENE_SETS)) {
  sc <- paste0(nm, "_score")
  p_mod <- ggplot(meta, aes(x = condition, y = .data[[sc]], fill = condition)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.4) +
    facet_wrap(~ cell_type, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = cond_pal) +
    labs(title = paste0(gsub("_", " ", nm), " Signature Score"),
         subtitle = paste0("Per cell, by condition  |  genes: ", paste(GENE_SETS[[nm]], collapse = ", ")),
         x = NULL, y = "Module Score") +
    theme_pub + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(DIRS$differential, paste0("module_score_", nm, ".pdf")), p_mod, width = 11, height = 8, dpi = PLOT$dpi)
  report_plots[[paste0("Module Score — ", nm)]] <- p_mod
}

if (length(report_plots)) {
  pdf(file.path(DIRS$differential, "differential_report.pdf"), width = 8, height = 6)
  for (p in report_plots) print(p)
  dev.off()
  message("\nReport saved: ", file.path(DIRS$differential, "differential_report.pdf"))
}
message("\n06b_differential.R complete. Outputs in: ", DIRS$differential)
