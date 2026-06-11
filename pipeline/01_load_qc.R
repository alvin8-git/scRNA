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
  if (is.list(counts)) {
    if (!"Gene Expression" %in% names(counts))
      stop("'Gene Expression' modality not found in: ", paste(names(counts), collapse = ", "))
    counts <- counts[["Gene Expression"]]
  }

  n_gem <- ncol(counts)   # raw GEM barcodes before any Seurat filtering
  seu <- CreateSeuratObject(counts = counts, project = sample_name,
                            min.cells = 3, min.features = 200)
  seu$sample     <- sample_name
  seu$orig.ident <- sample_name
  seu$n_gem      <- n_gem   # store per-cell so it survives QC subsetting
  seu <- RenameCells(seu, add.cell.id = sample_name)

  n_loaded      <- ncol(seu)
  n_low_gene    <- n_gem - n_loaded
  message("  GEM barcodes (CellRanger): ", n_gem)
  message("  After CreateSeuratObject  (min.features >= 200): ", n_loaded,
          "  [-", n_low_gene, " cells with < 200 genes]")
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
  ) & theme(legend.position = "none", plot.title = element_text(size = 11),
            axis.text.x = element_text(angle = 0, hjust = 0.5))
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

seu_list <- setNames(vector("list", length(SAMPLE_NAMES)), SAMPLE_NAMES)
for (nm in SAMPLE_NAMES) {
  cache_path <- file.path(SAMPLE_CACHE_DIR, nm, paste0(nm, "_filtered.rds"))
  hash_path  <- sub("\\.rds$", ".hash", cache_path)
  if (file.exists(cache_path)) {
    stored <- if (file.exists(hash_path)) readLines(hash_path) else ""
    if (stored == cache_hash(nm, "01")) {
      message("[CACHE HIT] ", nm, ": loading _filtered.rds from sample_cache/ (skipping load + QC filter)")
      seu_list[[nm]] <- readRDS(cache_path)
    } else {
      message("[CACHE STALE] ", nm, ": config changed — recomputing")
      seu_list[[nm]] <- add_qc_metrics(load_sample(nm, SAMPLE_PATHS[[nm]]))
    }
  } else {
    seu_list[[nm]] <- add_qc_metrics(load_sample(nm, SAMPLE_PATHS[[nm]]))
  }
}

# Capture counts BEFORE QC (after CreateSeuratObject load-filter)
n_gem_vec    <- sapply(seu_list, function(s) s$n_gem[1])
n_loaded_vec <- sapply(seu_list, ncol)

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

# Capture counts AFTER QC
n_qc_vec <- sapply(seu_list, ncol)

for (nm in SAMPLE_NAMES) {
  outdir <- file.path(DIRS$individual, nm)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  rds_path <- file.path(outdir, paste0(nm, "_filtered.rds"))
  saveRDS(seu_list[[nm]], rds_path)
  cache_path <- file.path(SAMPLE_CACHE_DIR, nm, paste0(nm, "_filtered.rds"))
  hash_path  <- sub("\\.rds$", ".hash", cache_path)
  dir.create(file.path(SAMPLE_CACHE_DIR, nm), recursive = TRUE, showWarnings = FALSE)
  file.copy(rds_path, cache_path, overwrite = TRUE)
  writeLines(cache_hash(nm, "01"), hash_path)
  message("  [CACHE] Saved ", nm, "_filtered.rds to sample_cache/")
}
gc()

qc_summary <- data.frame(
  Sample          = SAMPLE_NAMES,
  Cells           = sapply(seu_list, ncol),
  Genes           = sapply(seu_list, nrow),
  Median_nFeature = sapply(seu_list, function(x) round(median(x$nFeature_RNA))),
  Median_nCount   = sapply(seu_list, function(x) round(median(x$nCount_RNA))),
  Median_pct_mt   = sapply(seu_list, function(x) round(median(x$percent.mt), 2))
)
write.csv(qc_summary, file.path(DIRS$qc, "qc_summary_table.csv"), row.names = FALSE)

# --- Cell fate tracking: accounts for every barcode from GEM input ---
cell_fate <- data.frame(
  Sample                   = SAMPLE_NAMES,
  GEM_barcodes             = n_gem_vec,
  After_load               = n_loaded_vec,
  Removed_load_low_genes   = n_gem_vec - n_loaded_vec,   # < 200 genes (CreateSeuratObject)
  After_QC                 = n_qc_vec,
  Removed_QC               = n_loaded_vec - n_qc_vec,    # nFeature/nCount/MT thresholds
  # Doublet column filled in by 02_doublets.R
  After_doublet_removal    = NA_integer_,
  Removed_doublets         = NA_integer_,
  row.names                = NULL,
  stringsAsFactors         = FALSE
)
cell_fate$Total_removed  <- cell_fate$GEM_barcodes - cell_fate$After_QC
cell_fate$Pct_retained   <- round(cell_fate$After_QC / cell_fate$GEM_barcodes * 100, 1)
write.csv(cell_fate, file.path(DIRS$qc, "cell_fate.csv"), row.names = FALSE)

message("\n--- QC Summary ---")
print(qc_summary)
message("\n--- Cell Fate (load + QC) ---")
message("  Removal reasons:")
message("    load  : cells with < 200 genes (CreateSeuratObject min.features=200)")
message("    QC    : nFeature [", QC$min_features, ",", QC$max_features, "]  |  ",
        "nCount [", QC$min_counts, ",", QC$max_counts, "]  |  ",
        "MT <= ", QC$max_percent_mt, "%")
print(cell_fate[, c("Sample", "GEM_barcodes", "Removed_load_low_genes",
                    "After_load", "Removed_QC", "After_QC",
                    "Total_removed", "Pct_retained")])

save_report_pdf(report_plots, file.path(DIRS$qc, "qc_report.pdf"))
message("\n01_load_qc.R complete. Outputs in: ", DIRS$qc)
