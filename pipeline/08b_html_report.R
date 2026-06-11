#!/usr/bin/env Rscript
# 08b_html_report.R — assemble an interactive per-run HTML report from existing outputs.
#
# Reads a finished run directory (Results/results_*/) and renders one self-contained
# <run>_report.html via R Markdown + plotly. Additive: it never recomputes the
# analysis, only re-presents integrated_annotated.rds + the QC/annotation CSVs.
#
# Usage:
#   Rscript 08b_html_report.R <run_dir> [output.html] [--samples=A,B,...] [--max-cells=4000]
#
# <run_dir> is a Results/results_*/ directory (or its integrated/ subdir; both work).
# Default output: <run_dir>/reports/<run>_report.html

suppressWarnings(suppressMessages({
  library(Seurat)
}))

# ---- locate this script's directory (for the .Rmd template) ----
.pipeline_dir <- local({
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(f)) return(dirname(normalizePath(f)))
  a <- commandArgs(trailingOnly = FALSE)
  d <- sub("--file=", "", a[grep("--file=", a)])
  if (length(d) > 0) dirname(normalizePath(d[1])) else "."
})
TEMPLATE <- file.path(.pipeline_dir, "report_template.Rmd")

# ---- args ----
args <- commandArgs(trailingOnly = TRUE)
flags <- grep("^--", args, value = TRUE)
pos   <- args[!grepl("^--", args)]
getflag <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), flags, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
if (length(pos) < 1) stop("usage: Rscript 08b_html_report.R <run_dir> [output.html] [--samples=A,B] [--max-cells=N]")

run_dir <- normalizePath(pos[1], mustWork = TRUE)
# Accept either the run root or its integrated/ subdir.
if (basename(run_dir) == "integrated") run_dir <- dirname(run_dir)
integ_dir <- file.path(run_dir, "integrated")
run_name  <- sub("^results_", "", basename(run_dir))
run_name  <- sub("_filtered$", "", run_name)

max_cells   <- as.integer(getflag("max-cells", "4000"))
want_samples <- getflag("samples", NULL)
if (!is.null(want_samples)) want_samples <- trimws(strsplit(want_samples, ",")[[1]])

out_html <- if (length(pos) >= 2) pos[2] else
  file.path(run_dir, "reports", paste0(run_name, "_report.html"))
dir.create(dirname(out_html), showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(sprintf("[08b %s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

# ---- load the integrated object ----
rds <- file.path(integ_dir, "integrated_annotated.rds")
if (!file.exists(rds)) stop("not found: ", rds, " (run steps 04-05 first)")
msg("loading %s", rds)
obj <- readRDS(rds)
md  <- obj@meta.data

col <- function(cands) { hit <- intersect(cands, colnames(md)); if (length(hit)) hit[1] else NA_character_ }
sample_col <- col(c("sample", "orig.ident"))
ct_col     <- col(c("cell_type", "celltype", "consensus_label", "singler_label_clean", "singler_label"))
if (is.na(sample_col)) stop("no sample column in metadata")
if (is.na(ct_col))     stop("no cell-type column in metadata")

md$.sample    <- as.character(md[[sample_col]])
md$.cell_type <- as.character(md[[ct_col]])
md$.cell_type[is.na(md$.cell_type) | md$.cell_type == ""] <- "Unassigned"

samples_all <- sort(unique(md$.sample))
if (!is.null(want_samples)) {
  miss <- setdiff(want_samples, samples_all)
  if (length(miss)) warning("samples not in run, ignored: ", paste(miss, collapse = ", "))
  samples_all <- intersect(samples_all, want_samples)
  if (!length(samples_all)) stop("no requested samples present in run")
  md <- md[md$.sample %in% samples_all, , drop = FALSE]
}
msg("%d cells across %d samples: %s", nrow(md), length(samples_all), paste(samples_all, collapse = ", "))

# ---- UMAP embedding (integrated / shared coordinate space) ----
reds <- Reductions(obj)
umap_name <- {
  pref <- intersect(c("umap", "umap_harmony", "umap_integrated"), reds)
  if (length(pref)) pref[1] else grep("umap", reds, ignore.case = TRUE, value = TRUE)[1]
}
if (is.na(umap_name)) stop("no UMAP reduction found in object")
emb <- as.data.frame(Embeddings(obj, umap_name))[, 1:2]
colnames(emb) <- c("UMAP_1", "UMAP_2")
msg("using reduction '%s' for the atlas", umap_name)

# ---- per-cell frame for plotting, downsampled per sample ----
num <- function(cands) { c <- col(cands); if (is.na(c)) rep(NA_real_, nrow(md)) else suppressWarnings(as.numeric(md[[c]])) }
cells <- data.frame(
  sample        = md$.sample,
  cell_type     = md$.cell_type,
  UMAP_1        = emb[rownames(md), "UMAP_1"],
  UMAP_2        = emb[rownames(md), "UMAP_2"],
  nFeature_RNA  = num(c("nFeature_RNA")),
  nCount_RNA    = num(c("nCount_RNA")),
  percent.mt    = num(c("percent.mt", "percent_mt")),
  doublet_score = num(c("doublet_score", "scDblFinder.score")),
  stringsAsFactors = FALSE
)
set.seed(1)
if (is.finite(max_cells) && max_cells > 0) {
  keep <- unlist(lapply(split(seq_len(nrow(cells)), cells$sample), function(idx)
    if (length(idx) > max_cells) sample(idx, max_cells) else idx), use.names = FALSE)
  cells_plot <- cells[sort(keep), , drop = FALSE]
} else cells_plot <- cells
msg("plotting frame: %d cells (cap %d/sample)", nrow(cells_plot), max_cells)

# ---- proportions from the FULL data (not downsampled) ----
tab <- table(cells$cell_type, cells$sample)
prop_wide <- sweep(tab, 2, colSums(tab), "/") * 100   # cell_type x sample, % within sample
prop_long <- as.data.frame(as.table(prop_wide), stringsAsFactors = FALSE)
colnames(prop_long) <- c("cell_type", "sample", "pct")

delta <- NULL
if (length(samples_all) == 2) {
  s1 <- samples_all[1]; s2 <- samples_all[2]
  cts <- rownames(prop_wide)
  delta <- data.frame(
    cell_type = cts,
    a_pct = round(prop_wide[, s1], 1),
    b_pct = round(prop_wide[, s2], 1),
    delta_pts = round(prop_wide[, s2] - prop_wide[, s1], 1),
    row.names = NULL, stringsAsFactors = FALSE
  )
  names(delta)[2:3] <- c(paste0(s1, "_pct"), paste0(s2, "_pct"))
  delta <- delta[order(-abs(delta$delta_pts)), ]
  attr(delta, "samples") <- c(s1, s2)
}

# ---- summary table (prefer the CSVs the pipeline already wrote) ----
read_csv_safe <- function(p) if (file.exists(p)) tryCatch(read.csv(p, check.names = FALSE), error = function(e) NULL) else NULL
qc_sum   <- read_csv_safe(file.path(run_dir, "qc", "qc_summary_table.csv"))
fate     <- read_csv_safe(file.path(run_dir, "qc", "cell_fate.csv"))
clann    <- read_csv_safe(file.path(integ_dir, "integrated_cluster_markers.csv"))  # optional

# per-sample n cell types + top type computed from data (authoritative for this object)
by_s <- split(cells$cell_type, cells$sample)
n_types  <- sapply(by_s, function(x) length(unique(x)))
top_type <- sapply(by_s, function(x) { t <- sort(table(x), decreasing = TRUE); names(t)[1] })

summary <- data.frame(Sample = samples_all, stringsAsFactors = FALSE)
g <- function(df, key, val) if (!is.null(df) && all(c(key, val) %in% names(df))) df[[val]][match(summary$Sample, df[[key]])] else NA
summary$Cells           <- as.integer(table(cells$sample)[summary$Sample])
summary$Median_genes    <- g(qc_sum, "Sample", "Median_nFeature")
summary$Median_UMI      <- g(qc_sum, "Sample", "Median_nCount")
summary$Median_pct_MT   <- g(qc_sum, "Sample", "Median_pct_mt")
summary$Pct_retained    <- g(fate, "Sample", "Pct_retained")
if (!is.null(fate) && all(c("Removed_doublets", "After_QC") %in% names(fate))) {
  dr <- round(fate$Removed_doublets / pmax(fate$After_QC, 1) * 100, 1)
  summary$Doublet_pct <- dr[match(summary$Sample, fate$Sample)]
}
summary$N_cell_types <- n_types[summary$Sample]
summary$Top_type     <- top_type[summary$Sample]

# ---- bundle + render ----
bundle <- list(
  run_name = run_name, generated = as.character(Sys.time()),
  samples = samples_all, n_cells_total = nrow(cells),
  cells = cells_plot, prop_long = prop_long, prop_wide = prop_wide,
  delta = delta, summary = summary, umap_name = umap_name
)
bpath <- file.path(tempdir(), paste0("report_bundle_", Sys.getpid(), ".rds"))
saveRDS(bundle, bpath)

# point rmarkdown at the conda env's pandoc (not on PATH under a bare Rscript call)
if (!rmarkdown::pandoc_available()) {
  env_bin <- file.path(dirname(dirname(Sys.getenv("R_HOME"))), "bin")
  if (file.exists(file.path(env_bin, "pandoc"))) {
    Sys.setenv(RSTUDIO_PANDOC = env_bin); rmarkdown::find_pandoc(dir = env_bin)
  }
}
msg("pandoc: %s", tryCatch(as.character(rmarkdown::pandoc_version()), error = function(e) "NOT FOUND"))
msg("rendering -> %s", out_html)

rmarkdown::render(
  TEMPLATE,
  output_file = basename(out_html),
  output_dir  = dirname(out_html),
  params      = list(bundle = bpath, run_name = run_name),
  quiet       = TRUE,
  envir       = new.env()
)
unlink(bpath)
msg("DONE: %s", out_html)
