# =============================================================================
# 09_bootstrap_proportions.R
# Normalised cell type proportion comparison across samples using:
#   1. Bootstrap downsampling  — resample each sample to the smallest n, 1000x
#   2. Multinomial 95% CI      — analytical CI for each sample's observed proportions
#   3. Pairwise proportion test — chi-squared goodness-of-fit per cell type pair
# Output: bootstrap_proportions_report.pdf + bootstrap_summary.csv
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
message("\nPairwise chi-squared tests (observed cell counts):")
print(chisq_df)

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
    labs(title = "Pairwise Chi-squared Test — Cell Type Composition",
         subtitle = "p < 0.05 = significantly different composition",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text = element_text(size = 9))
  report_plots[["Pairwise Chi-squared Test"]] <- p_chisq
}

save_report_pdf(report_plots,
                file.path(DIRS$reports, "bootstrap_proportions_report.pdf"))

message("\n09_bootstrap_proportions.R complete.")
message("Outputs:")
message("  ", file.path(DIRS$reports, "bootstrap_proportions_report.pdf"))
message("  ", file.path(DIRS$reports, "bootstrap_summary.csv"))
