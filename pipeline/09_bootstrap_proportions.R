# =============================================================================
# 09_bootstrap_proportions.R
# Normalised cell type proportion comparison across samples using:
#   1. Bootstrap downsampling  — resample each sample to the smallest n, 1000x
#   2. Multinomial 95% CI      — analytical CI for each sample's observed proportions
#   3. Pairwise chi-squared     — descriptive only (cells are not replicates)
#   4. Composition test         — replicate-aware, sample-level, propeller-style
#                                 (needs SCRNA_CONDITION and >= 2 samples per condition)
# Output: bootstrap_proportions_report.pdf + bootstrap_summary.csv + composition_test.csv
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
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

seu <- readRDS(file.path(DIRS$integrated, "integrated_annotated.rds"))
seu <- apply_reference_labels(seu)   # use frozen-reference labels when 05r has run (additive)
.lblsrc <- getOption("scrna.label_source", "de-novo")
message("Loaded: ", ncol(seu), " cells across samples: ",
        paste(unique(seu$sample), collapse = ", "), "  | labels: ", .lblsrc)

meta       <- seu@meta.data
n_boot     <- 1000
set.seed(42)

# ── Observed proportions per sample ─────────────────────────────────────────
obs <- meta %>%
  group_by(sample, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(total  = sum(n),
         prop   = n / total,
         # Multinomial 95% CI (Wilson-score per cell type)
         ci_lo  = pmax(0, prop - 1.96 * sqrt(prop * (1 - prop) / total)),
         ci_hi  = pmin(1, prop + 1.96 * sqrt(prop * (1 - prop) / total))) %>%
  ungroup()

sample_sizes <- obs %>% group_by(sample) %>% summarise(total = first(total))
n_min        <- min(sample_sizes$total)
message("Smallest sample: ", n_min, " cells — bootstrapping all samples to this depth")

# ── Bootstrap downsampling ───────────────────────────────────────────────────
all_cell_types <- sort(unique(meta$cell_type))

boot_results <- lapply(unique(meta$sample), function(s) {
  cells <- meta$cell_type[meta$sample == s]
  draws <- replicate(n_boot, {
    samp  <- sample(cells, size = n_min, replace = FALSE)
    table(factor(samp, levels = all_cell_types)) / n_min
  })
  # draws is n_celltypes × n_boot
  data.frame(
    sample    = s,
    cell_type = all_cell_types,
    boot_mean = rowMeans(draws),
    boot_lo   = apply(draws, 1, quantile, 0.025),
    boot_hi   = apply(draws, 1, quantile, 0.975),
    stringsAsFactors = FALSE
  )
})
boot_df <- bind_rows(boot_results)

# ── Summary CSV ──────────────────────────────────────────────────────────────
summary_df <- obs %>%
  select(sample, cell_type, total, n, prop, ci_lo, ci_hi) %>%
  left_join(boot_df %>% select(sample, cell_type, boot_mean, boot_lo, boot_hi),
            by = c("sample", "cell_type")) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

write.csv(summary_df,
          file.path(DIRS$reports, "bootstrap_summary.csv"),
          row.names = FALSE)
message("Summary CSV written.")

# ── Pairwise chi-squared test (observed counts) ──────────────────────────────
samples     <- unique(meta$sample)
ct_counts   <- obs %>% select(sample, cell_type, n) %>%
               pivot_wider(names_from = cell_type, values_from = n, values_fill = 0)
pair_tests  <- combn(samples, 2, simplify = FALSE)

chisq_results <- lapply(pair_tests, function(pair) {
  tmp <- ct_counts %>% filter(sample %in% pair) %>% as.data.frame()
  rownames(tmp) <- tmp$sample
  tmp$sample    <- NULL
  mat           <- as.matrix(tmp)
  # Drop cell types absent in both
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  test <- tryCatch(chisq.test(mat), error = function(e) NULL)
  data.frame(
    sample_A  = pair[1],
    sample_B  = pair[2],
    chi_sq    = if (!is.null(test)) round(test$statistic, 2) else NA,
    df        = if (!is.null(test)) test$parameter else NA,
    p_value   = if (!is.null(test)) signif(test$p.value, 3) else NA,
    stringsAsFactors = FALSE
  )
})
chisq_df <- bind_rows(chisq_results)
message("\nPairwise chi-squared tests (observed cell counts — descriptive only, cells are not replicates):")
print(chisq_df)

# ── Replicate-aware composition test (propeller-style) ───────────────────────
# Unit of replication is the SAMPLE, not the cell. Per-sample proportions are
# arcsin-sqrt transformed (variance-stabilising for proportions) and fitted with
# limma (~ condition) with empirical-Bayes moderation — the propeller method of
# Phipson et al. 2022 without the speckle dependency. Needs >= 2 conditions and
# >= 2 samples in every condition; otherwise the run is descriptive only.
# ponytail: t-test on asin(sqrt(p)) per cell type; add logit variant if ever needed.
cond <- SAMPLE_CONDITIONS[samples]
if (any(is.na(cond))) cond[is.na(cond)] <- "condition1"
n_per_cond <- table(cond)
comp_test  <- NULL
if (length(n_per_cond) >= 2 && all(n_per_cond >= 2)) {
  prop_mat <- obs %>% select(sample, cell_type, prop) %>%
    pivot_wider(names_from = sample, values_from = prop, values_fill = 0) %>%
    as.data.frame()
  rownames(prop_mat) <- prop_mat$cell_type; prop_mat$cell_type <- NULL
  prop_mat <- as.matrix(prop_mat[, samples, drop = FALSE])
  lv       <- unique(cond)
  design   <- model.matrix(~ 0 + factor(cond, levels = lv)); colnames(design) <- lv
  pairs    <- combn(lv, 2, simplify = FALSE)                       # every pair: later vs earlier
  cm       <- sapply(pairs, function(p) { v <- setNames(numeric(length(lv)), lv); v[p[2]] <- 1; v[p[1]] <- -1; v })
  colnames(cm) <- sapply(pairs, function(p) paste0(p[2], "_vs_", p[1]))
  fit <- limma::lmFit(asin(sqrt(prop_mat)), design)
  fit <- limma::eBayes(limma::contrasts.fit(fit, cm), robust = TRUE, trend = FALSE)
  comp_test <- bind_rows(lapply(colnames(cm), function(cf) {
    tt  <- limma::topTable(fit, coef = cf, number = Inf, sort.by = "none")
    alt <- sub("_vs_.*", "", cf); ref <- sub(".*_vs_", "", cf)
    data.frame(cell_type = rownames(tt), contrast = cf,
               mean_prop_ref = rowMeans(prop_mat[rownames(tt), cond == ref, drop = FALSE]),
               mean_prop_alt = rowMeans(prop_mat[rownames(tt), cond == alt, drop = FALSE]),
               t = tt$t, p_value = tt$P.Value, fdr = tt$adj.P.Val, stringsAsFactors = FALSE)
  })) %>% mutate(across(where(is.numeric), ~ signif(.x, 4))) %>% arrange(fdr)
  write.csv(comp_test, file.path(DIRS$reports, "composition_test.csv"), row.names = FALSE)
  message("\nReplicate-aware composition test (", paste(names(n_per_cond), n_per_cond, sep = "=", collapse = ", "), "):")
  print(comp_test, row.names = FALSE)
} else {
  message("\nComposition test skipped — needs >= 2 conditions with >= 2 samples each ",
          "(have: ", paste(names(n_per_cond), n_per_cond, sep = "=", collapse = ", "),
          "). Set SCRNA_CONDITION=sample=label,... Proportions are DESCRIPTIVE only.")
}

# ── Plots ────────────────────────────────────────────────────────────────────
ct_colors_use <- CELLTYPE_COLORS[names(CELLTYPE_COLORS) %in% all_cell_types]
extra <- setdiff(all_cell_types, c(names(CELLTYPE_COLORS), "Unassigned", "Unknown"))
if (length(extra) > 0)
  ct_colors_use <- c(ct_colors_use, setNames(hue_pal()(length(extra)), extra))

report_plots <- list()

# 1. Observed stacked bar with multinomial CI (error bars on each segment top)
p_obs <- ggplot(obs, aes(x = sample, y = prop, fill = cell_type)) +
  geom_col(width = 0.6, position = "stack") +
  geom_text(aes(label = ifelse(prop >= 0.02, paste0(round(prop * 100, 1), "%"), "")),
            position = position_stack(vjust = 0.5),
            size = 2.8, color = "white", fontface = "bold") +
  scale_fill_manual(values = ct_colors_use, name = "Cell Type") +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Observed Cell Type Proportions per Sample",
       subtitle = paste0("Labels: ", .lblsrc, "  |  Error bars = multinomial 95% CI  |  n cells: ",
                         paste(paste0(sample_sizes$sample, "=", sample_sizes$total),
                               collapse = ", ")),
       x = NULL, y = "Proportion") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
        legend.text = element_text(size = 8))
report_plots[["Observed Proportions (all cells)"]] <- p_obs

# 2. Bootstrap-downsampled proportions (all samples at n_min depth)
p_boot <- ggplot(boot_df, aes(x = sample, y = boot_mean, fill = cell_type)) +
  geom_col(width = 0.6, position = "stack") +
  geom_text(aes(label = ifelse(boot_mean >= 0.02,
                               paste0(round(boot_mean * 100, 1), "%"), "")),
            position = position_stack(vjust = 0.5),
            size = 2.8, color = "white", fontface = "bold") +
  scale_fill_manual(values = ct_colors_use, name = "Cell Type") +
  scale_y_continuous(labels = percent_format()) +
  labs(title = paste0("Bootstrap-Normalised Proportions (", n_min, " cells × ", n_boot, " draws)"),
       subtitle = paste0("Labels: ", .lblsrc, " — all samples downsampled to smallest n (removes capture-depth bias)"),
       x = NULL, y = "Proportion") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
        legend.text = element_text(size = 8))
report_plots[["Bootstrap-Normalised Proportions"]] <- p_boot

# 3. Per-cell-type dot plot with bootstrap CI across samples
p_dot <- ggplot(boot_df %>% filter(boot_mean > 0),
                aes(x = sample, y = boot_mean, color = cell_type)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = boot_lo, ymax = boot_hi), width = 0.2) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 4) +
  scale_color_manual(values = ct_colors_use, guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(title = "Bootstrap 95% CI per Cell Type × Sample",
       subtitle = paste0("n = ", n_min, " cells per draw, ", n_boot, " iterations"),
       x = NULL, y = "Proportion") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7),
        strip.text = element_text(size = 7, face = "bold"))
report_plots[["Bootstrap CI per Cell Type"]] <- p_dot

# 4. Chi-squared summary table as plot
if (nrow(chisq_df) > 0) {
  chisq_df$label <- ifelse(!is.na(chisq_df$p_value),
                            paste0("χ²=", chisq_df$chi_sq,
                                   "  p=", chisq_df$p_value),
                            "test failed")
  p_chisq <- ggplot(chisq_df, aes(x = sample_A, y = sample_B, fill = p_value)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = label), size = 3) +
    scale_fill_gradient2(low = "#E64B35", mid = "#FFFBCC", high = "#4DBBD5",
                         midpoint = 0.05, name = "p-value",
                         limits = c(0, 1), na.value = "grey80") +
    labs(title = "Pairwise Chi-squared Test — Cell Type Composition (descriptive)",
         subtitle = "Treats cells as replicates, so p-values are anti-conservative. Use the composition test below for inference.",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text = element_text(size = 9))
  report_plots[["Pairwise Chi-squared Test"]] <- p_chisq
}

# 5. Replicate-aware composition test
if (!is.null(comp_test) && nrow(comp_test) > 0) {
  ct_df <- comp_test %>%
    mutate(cell_type = factor(cell_type, levels = rev(unique(cell_type[order(fdr)]))),
           sig = ifelse(fdr < 0.05, "FDR < 0.05", "n.s."),
           diff_pp = 100 * (mean_prop_alt - mean_prop_ref))
  p_comp <- ggplot(ct_df, aes(x = diff_pp, y = cell_type, colour = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_point(size = 3) +
    geom_text(aes(label = paste0("FDR=", signif(fdr, 2))), size = 2.6, vjust = -1) +
    facet_wrap(~ contrast, scales = "free_x") +
    scale_colour_manual(values = c("FDR < 0.05" = "#E64B35", "n.s." = "grey50"), name = NULL) +
    labs(title = "Composition Test — sample-level, replicate-aware (propeller-style)",
         subtitle = paste0("asin(sqrt(prop)) per sample, limma moderated t; replicates per condition: ",
                           paste(names(n_per_cond), n_per_cond, sep = "=", collapse = ", ")),
         x = "Difference in mean proportion (percentage points)", y = NULL) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")
  report_plots[["Composition Test (replicate-aware)"]] <- p_comp
}

save_report_pdf(report_plots,
                file.path(DIRS$reports, "bootstrap_proportions_report.pdf"))

message("\n09_bootstrap_proportions.R complete.")
message("Outputs:")
message("  ", file.path(DIRS$reports, "bootstrap_proportions_report.pdf"))
message("  ", file.path(DIRS$reports, "bootstrap_summary.csv"))
if (!is.null(comp_test)) message("  ", file.path(DIRS$reports, "composition_test.csv"))
