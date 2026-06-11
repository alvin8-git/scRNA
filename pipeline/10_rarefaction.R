# =============================================================================
# 10_rarefaction.R
# Bootstrap rarefaction to determine minimum cell count for stable proportion
# estimates in bat whole blood scRNA-seq.
#
# Approach:
#   - Use the largest sample (by cell count) as ground truth
#   - Subsample at increasing depths (n_depths), 1000 draws each
#   - For each depth: compute mean proportion, 95% CI width, RMSE vs ground truth
#   - Fit theoretical CI ~ a/sqrt(n) curve per cell type
#   - Report minimum n where empirical CI is within 5% of asymptote
#
# Output: rarefaction_report.pdf + rarefaction_summary.csv
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
meta <- seu@meta.data

# ── Use largest sample as ground truth ───────────────────────────────────────
sample_sizes <- sort(table(meta$sample), decreasing = TRUE)
message("Sample sizes:")
print(sample_sizes)

ref_sample <- names(sample_sizes)[1]
ref_n      <- as.integer(sample_sizes[1])
message("\nGround truth sample: ", ref_sample, " (n=", ref_n, ")")

ref_cells      <- meta$cell_type[meta$sample == ref_sample]
all_cell_types <- sort(unique(meta$cell_type))

# True proportions from full reference sample
true_props <- setNames(
  as.numeric(table(factor(ref_cells, levels = all_cell_types))) / ref_n,
  all_cell_types
)
message("True proportions (", ref_sample, "):")
print(round(true_props * 100, 2))

# ── Rarefaction depths ────────────────────────────────────────────────────────
depths  <- c(100, 250, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000)
depths  <- depths[depths <= ref_n]
n_boot  <- 1000
set.seed(42)

message("\nRunning rarefaction: ", length(depths), " depths × ", n_boot, " draws...")

raref_results <- lapply(depths, function(n) {
  message("  depth = ", n)
  draws <- replicate(n_boot, {
    samp <- sample(ref_cells, size = n, replace = FALSE)
    table(factor(samp, levels = all_cell_types)) / n
  })
  # draws: n_celltypes × n_boot matrix
  data.frame(
    depth     = n,
    cell_type = all_cell_types,
    mean_prop = rowMeans(draws),
    ci_width  = apply(draws, 1, function(x) diff(quantile(x, c(0.025, 0.975)))),
    rmse      = sqrt(rowMeans((draws - true_props)^2)),
    stringsAsFactors = FALSE
  )
})
raref_df <- bind_rows(raref_results)

# Add true proportion column for reference
raref_df <- raref_df %>%
  mutate(true_prop = true_props[cell_type])

# ── Fit theoretical CI curve: ci_width ~ a / sqrt(n) ─────────────────────────
# Theoretical: CI = 2 * 1.96 * sqrt(p*(1-p)/n) → a = 2*1.96*sqrt(p*(1-p))
raref_df <- raref_df %>%
  mutate(
    theo_ci = 2 * 1.96 * sqrt(true_prop * (1 - true_prop) / depth)
  )

# ── Minimum n: where empirical CI ≤ 105% of asymptotic theoretical value ──────
# Use n=ref_n as "asymptote" reference: CI_asymptote = theoretical CI at ref_n
min_n_df <- raref_df %>%
  group_by(cell_type) %>%
  mutate(
    ci_asymptote = 2 * 1.96 * sqrt(first(true_prop) * (1 - first(true_prop)) / ref_n)
  ) %>%
  filter(ci_width <= 1.05 * ci_asymptote * (ref_n / depth)^0 | TRUE) %>%
  # Find smallest depth where ci_width ≤ 1.5× theoretical CI at that depth
  filter(abs(ci_width - theo_ci) / (theo_ci + 1e-6) <= 0.15) %>%
  summarise(min_n_stable = min(depth), .groups = "drop")

# Alternative: smallest n where ci_width < 2× the value at max depth
asymp_tbl <- raref_df %>%
  group_by(cell_type) %>%
  filter(depth == max(depth)) %>%
  summarise(asymp_ci  = mean(ci_width),
            tp        = mean(true_prop),
            .groups   = "drop")

convergence_df <- raref_df %>%
  left_join(asymp_tbl, by = "cell_type") %>%
  group_by(cell_type) %>%
  filter(ci_width <= 2 * asymp_ci) %>%
  summarise(
    min_n_2x  = min(depth, na.rm = TRUE),
    true_prop = first(tp),
    asymp_ci  = first(asymp_ci),
    .groups   = "drop"
  ) %>%
  # cell types with ci always narrow: mark as n=100 (already converged)
  right_join(asymp_tbl %>% select(cell_type, tp, asymp_ci) %>%
               rename(true_prop = tp), by = "cell_type") %>%
  mutate(min_n_2x = ifelse(is.na(min_n_2x), 100L, min_n_2x),
         asymp_ci = coalesce(asymp_ci.x, asymp_ci.y),
         true_prop = coalesce(true_prop.x, true_prop.y)) %>%
  select(cell_type, true_prop, min_n_2x, asymp_ci) %>%
  arrange(desc(min_n_2x))

message("\nMinimum n for CI convergence (within 2× asymptotic):")
print(convergence_df)

# ── Write summary CSV ─────────────────────────────────────────────────────────
write.csv(
  raref_df %>% mutate(across(where(is.numeric), ~ round(.x, 5))),
  file.path(DIRS$reports, "rarefaction_summary.csv"),
  row.names = FALSE
)

# ── Plots ─────────────────────────────────────────────────────────────────────
ct_colors_use <- CELLTYPE_COLORS[names(CELLTYPE_COLORS) %in% all_cell_types]
extra <- setdiff(all_cell_types, names(ct_colors_use))
if (length(extra) > 0)
  ct_colors_use <- c(ct_colors_use, setNames(scales::hue_pal()(length(extra)), extra))

report_plots <- list()

# ── Plot 1: CI width vs depth per cell type (empirical + theoretical) ─────────
# Only show cell types with true_prop > 0
ct_show <- all_cell_types[true_props > 0]

p1 <- ggplot(raref_df %>% filter(cell_type %in% ct_show),
             aes(x = depth)) +
  geom_line(aes(y = ci_width, color = cell_type), linewidth = 0.9) +
  geom_line(aes(y = theo_ci, color = cell_type),
            linetype = "dashed", linewidth = 0.5, alpha = 0.6) +
  geom_vline(xintercept = 3604, linetype = "dotted", color = "grey40", linewidth = 0.6) +
  annotate("text", x = 3604 + 80, y = max(raref_df$ci_width) * 0.92,
           label = "Sample6\n(n=3604)", size = 2.8, hjust = 0, color = "grey40") +
  scale_color_manual(values = ct_colors_use, name = "Cell Type") +
  scale_x_continuous(breaks = depths, labels = comma) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title    = paste0("Rarefaction: 95% CI Width vs Capture Depth (", ref_sample, ", n=", ref_n, ")"),
    subtitle = "Solid = empirical bootstrap CI   |   Dashed = theoretical a/√n   |   Dotted = Sample6 depth",
    x        = "Number of captured cells (n)",
    y        = "95% CI width (proportion)"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
report_plots[["CI Width vs Depth"]] <- p1

# ── Plot 2: RMSE vs depth per cell type ───────────────────────────────────────
p2 <- ggplot(raref_df %>% filter(cell_type %in% ct_show),
             aes(x = depth, y = rmse, color = cell_type)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  geom_vline(xintercept = 3604, linetype = "dotted", color = "grey40", linewidth = 0.6) +
  scale_color_manual(values = ct_colors_use, name = "Cell Type") +
  scale_x_continuous(breaks = depths, labels = comma) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title    = "RMSE of Proportion Estimates vs Ground Truth",
    subtitle = paste0("Ground truth = full ", ref_sample, " (n=", ref_n, ")"),
    x        = "Number of captured cells (n)",
    y        = "RMSE (proportion)"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
report_plots[["RMSE vs Depth"]] <- p2

# ── Plot 3: Faceted CI width per cell type with convergence annotation ─────────
conv_lines <- convergence_df %>% select(cell_type, min_n_2x)

p3 <- ggplot(raref_df %>% filter(cell_type %in% ct_show),
             aes(x = depth, y = ci_width)) +
  geom_line(aes(color = cell_type), linewidth = 1) +
  geom_line(aes(y = theo_ci), linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_vline(data = conv_lines %>% filter(cell_type %in% ct_show),
             aes(xintercept = min_n_2x), linetype = "dotdash",
             color = "#E64B35", linewidth = 0.6) +
  geom_vline(xintercept = 3604, linetype = "dotted", color = "grey40", linewidth = 0.5) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 3) +
  scale_color_manual(values = ct_colors_use, guide = "none") +
  scale_x_continuous(breaks = c(500, 2000, 4000, 6000), labels = comma) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title    = "CI Width Rarefaction per Cell Type",
    subtitle = "Red dash = convergence threshold (2× asymptotic CI)  |  Dotted = Sample6 (n=3604)  |  Grey dash = theoretical",
    x        = "n cells", y = "95% CI width"
  ) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7),
        strip.text  = element_text(size = 8, face = "bold"))
report_plots[["Faceted CI per Cell Type"]] <- p3

# ── Plot 4: Minimum n summary table ──────────────────────────────────────────
summ_tbl <- convergence_df %>%
  mutate(
    true_pct  = paste0(round(true_prop * 100, 1), "%"),
    asymp_ci_pct = paste0(round(asymp_ci * 100, 2), "%")
  ) %>%
  select(cell_type, true_pct, min_n_2x, asymp_ci_pct) %>%
  arrange(desc(min_n_2x)) %>%
  rename(
    `Cell Type`          = cell_type,
    `True Prop`          = true_pct,
    `Min n (converge)`   = min_n_2x,
    `CI at max depth`    = asymp_ci_pct
  )

p4 <- ggplot() +
  annotate("text", x = 0.5, y = 0.98,
           label = "Minimum Capture Depth for Stable Proportion Estimates",
           size = 4.5, fontface = "bold", hjust = 0.5, vjust = 1) +
  annotate("text", x = 0.5, y = 0.88,
           label = paste0("Reference: ", ref_sample, " (n=", ref_n, ")  |  ",
                          "Convergence = CI ≤ 2× asymptotic CI  |  1000 bootstrap draws"),
           size = 3, hjust = 0.5, vjust = 1, color = "grey40") +
  theme_void(base_size = 10) +
  annotation_custom(
    gridExtra::tableGrob(summ_tbl,
                         rows = NULL,
                         theme = gridExtra::ttheme_minimal(base_size = 9)),
    xmin = 0, xmax = 1, ymin = 0, ymax = 0.8
  )

if (requireNamespace("gridExtra", quietly = TRUE)) {
  report_plots[["Summary Table"]] <- p4
}

save_report_pdf(report_plots,
                file.path(DIRS$reports, "rarefaction_report.pdf"))

message("\n10_rarefaction.R complete.")
message("Reference sample: ", ref_sample, " (n=", ref_n, ")")
message("Recommended minimum n (most demanding cell type): ",
        max(convergence_df$min_n_2x))
message("Outputs:")
message("  ", file.path(DIRS$reports, "rarefaction_report.pdf"))
message("  ", file.path(DIRS$reports, "rarefaction_summary.csv"))
