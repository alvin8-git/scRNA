# =============================================================================
# 01_load_qc.R - Load expression matrices, compute QC metrics,
#                generate diagnostic plots, filter cells.
#                Outputs: per-sample *_filtered.rds + results/qc/qc_report.pdf
# =============================================================================
source("/data/alvin/scRNA/pipeline/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

# =============================================================================
# FUNCTIONS
# =============================================================================

load_sample <- function(sample_name, data_path) {
  message("Loading: ", sample_name, " from ", data_path)
  counts <- Read10X(data.dir = data_path, gene.column = 2,
                    cell.column = 1, unique.features = TRUE)
  if (is.list(counts)) counts <- counts[["Gene Expression"]]
  seu <- CreateSeuratObject(counts = counts, project = sample_name,
                            min.cells = 3, min.features = 200)
  seu$sample     <- sample_name
  seu$orig.ident <- sample_name
  seu <- RenameCells(seu, add.cell.id = sample_name)
  message("  Loaded: ", ncol(seu), " cells x ", nrow(seu), " genes")
  seu
}

add_qc_metrics <- function(seu) {
  seu[["percent.mt"]]   <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
  seu
}

# Returns a named list of ggplot objects (also saves individual files)
plot_qc <- function(seu, sample_name) {
  outdir <- DIRS$qc
  plots  <- list()

  p_vln <- VlnPlot(seu,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol     = 3, pt.size = 0.2, group.by = "orig.ident"
  ) & theme(legend.position = "none", plot.title = element_text(size = 11))
  ggsave(file.path(outdir, paste0(sample_name, "_violin_qc.pdf")),
         p_vln, width = 12, height = 5)
  plots[[paste0(sample_name, "  -  QC Violin (nFeature / nCount / %MT)")]] <-
    set_page(p_vln, pw = 11, ph = 5)

  p_sc1 <- FeatureScatter(seu, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident") +
    geom_hline(yintercept = c(QC$min_features, QC$max_features),
               linetype = "dashed", color = "red", linewidth = 0.5) +
    theme_classic() + theme(legend.position = "none") +
    labs(title = paste0(sample_name, ": UMI vs Genes"))

  p_sc2 <- FeatureScatter(seu, "nCount_RNA", "percent.mt", group.by = "orig.ident") +
    geom_hline(yintercept = QC$max_percent_mt,
               linetype = "dashed", color = "red", linewidth = 0.5) +
    theme_classic() + theme(legend.position = "none") +
    labs(title = paste0(sample_name, ": UMI vs %MT"))

  p_scatter <- p_sc1 + p_sc2
  ggsave(file.path(outdir, paste0(sample_name, "_scatter_qc.pdf")),
         p_scatter, width = 12, height = 5)
  plots[[paste0(sample_name, "  -  QC Scatter (UMI vs Genes & UMI vs %MT)")]] <-
    set_page(p_scatter, pw = 11, ph = 5)

  message("  Pre-filter: ", ncol(seu), " cells | ",
          "nFeature median: ", round(median(seu$nFeature_RNA)),
          " | %MT median: ", round(median(seu$percent.mt), 2), "%")
  plots
}

filter_cells <- function(seu, sample_name) {
  n_before <- ncol(seu)
  seu <- subset(seu,
    subset =
      nFeature_RNA >= QC$min_features &
      nFeature_RNA <= QC$max_features &
      nCount_RNA   >= QC$min_counts   &
      nCount_RNA   <= QC$max_counts   &
      percent.mt   <= QC$max_percent_mt
  )
  n_after   <- ncol(seu)
  n_removed <- n_before - n_after
  message("  Filter: ", n_before, " -> ", n_after,
          " cells (removed ", n_removed,
          ", ", round(n_removed / n_before * 100, 1), "%)")
  seu
}

# =============================================================================
# MAIN
# =============================================================================

seu_list <- mapply(load_sample, names(SAMPLE_PATHS), SAMPLE_PATHS,
                   SIMPLIFY = FALSE)
seu_list <- lapply(seu_list, add_qc_metrics)

message("\n--- Generating pre-filter QC plots ---")
report_plots <- list()
for (nm in SAMPLE_NAMES) {
  plots <- plot_qc(seu_list[[nm]], nm)
  report_plots <- c(report_plots, plots)
}

message("\n--- Applying QC filters ---")
message("Thresholds: nFeature [", QC$min_features, ",", QC$max_features, "] | ",
        "nCount [", QC$min_counts, ",", QC$max_counts, "] | ",
        "MT <= ", QC$max_percent_mt, "%")
seu_list <- mapply(filter_cells, seu_list, SAMPLE_NAMES, SIMPLIFY = FALSE)

for (nm in SAMPLE_NAMES) {
  outdir <- file.path(DIRS$individual, nm)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(seu_list[[nm]], file.path(outdir, paste0(nm, "_filtered.rds")))
}

qc_summary <- data.frame(
  Sample          = SAMPLE_NAMES,
  Cells           = sapply(seu_list, ncol),
  Genes           = sapply(seu_list, nrow),
  Median_nFeature = sapply(seu_list, function(x) round(median(x$nFeature_RNA))),
  Median_nCount   = sapply(seu_list, function(x) round(median(x$nCount_RNA))),
  Median_pct_mt   = sapply(seu_list, function(x) round(median(x$percent.mt), 2))
)
write.csv(qc_summary, file.path(DIRS$qc, "qc_summary_table.csv"), row.names = FALSE)

message("\n--- QC Summary ---")
print(qc_summary)

save_report_pdf(report_plots, file.path(DIRS$qc, "qc_report.pdf"))
message("\n01_load_qc.R complete. Outputs in: ", DIRS$qc)
