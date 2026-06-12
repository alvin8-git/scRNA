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
  library(ggplot2)
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

max_cells   <- as.integer(getflag("max-cells", "6000"))
want_samples <- getflag("samples", NULL)
if (!is.null(want_samples)) want_samples <- trimws(strsplit(want_samples, ",")[[1]])
lite <- "--lite" %in% flags

# Long run names (every sample dash-joined) make a >300-char path that Windows
# cannot open over a share (MAX_PATH = 260). Shorten the default report filename to
# <firstSample>_<N>samples when the run name is long; an explicit output path wins.
short_report_name <- function(name) {
  if (nchar(name) <= 64) return(name)
  toks <- strsplit(name, "-", fixed = TRUE)[[1]]
  paste0(gsub("[^A-Za-z0-9]", "", toks[1]), "_", length(toks), "samples")
}
out_html <- if (length(pos) >= 2) pos[2] else
  file.path(run_dir, "reports", paste0(short_report_name(run_name), "_report.html"))
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
# round coords/QC so the interactive plotly widgets are not bloated by 17-digit floats
cells_plot$UMAP_1 <- round(cells_plot$UMAP_1, 3)
cells_plot$UMAP_2 <- round(cells_plot$UMAP_2, 3)
for (cc in c("nFeature_RNA","nCount_RNA","percent.mt","doublet_score"))
  if (cc %in% names(cells_plot)) cells_plot[[cc]] <- round(cells_plot[[cc]], 3)
msg("plotting frame: %d cells (cap %d/sample)", nrow(cells_plot), max_cells)

# shared palette (mirrors the template) for any PNGs we pre-render here
PAL_base <- c("#4e79a7","#f28e2b","#e15759","#76b7b2","#59a14f","#edc948","#b07aa1",
              "#ff9da7","#9c755f","#bab0ac","#86bcb6","#d37295","#fabfd2","#8cd17d","#499894")
ctlev <- sort(unique(cells_plot$cell_type))
PAL   <- setNames(rep(PAL_base, length.out = length(ctlev)), ctlev)
to_uri <- function(g, w, h, dpi = 72) {
  tf <- tempfile(fileext = ".png")
  ggplot2::ggsave(tf, plot = g, width = w, height = h, dpi = dpi, bg = "white", limitsize = FALSE)
  on.exit(unlink(tf), add = TRUE); knitr::image_uri(tf)
}

# Past ~8 samples, 14+ live WebGL UMAP widgets blow up the file (16 MB of JSON) and
# exceed the browser's WebGL-context limit, so the report won't open. For large runs
# render the per-sample UMAP panels as light static PNGs instead.
umap_static <- length(samples_all) > 8
umap_panels <- NULL
if (umap_static) {
  umap_panels <- lapply(samples_all, function(s) {
    dd <- cells_plot[cells_plot$sample == s, , drop = FALSE]
    g <- ggplot(dd, aes(UMAP_1, UMAP_2, color = cell_type)) +
      geom_point(size = 0.25, alpha = 0.7, stroke = 0) +
      scale_color_manual(values = PAL, guide = "none") +
      theme_void(base_size = 9)
    to_uri(g, w = 3.1, h = 2.9)
  })
  names(umap_panels) <- samples_all
  msg("static UMAP panels: %d (n_samples > 8)", length(umap_panels))
}

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

# median %ribo straight from metadata (qc CSV may not carry it)
if ("percent.ribo" %in% colnames(md)) {
  rib <- tapply(suppressWarnings(as.numeric(md$percent.ribo)), md$.sample, median, na.rm = TRUE)
  summary$Median_pct_ribo <- round(rib[summary$Sample], 2)
}
# highly-variable genes per sample (from the per-sample objects, if present)
hvg_n <- vapply(samples_all, function(s) {
  f <- list.files(file.path(run_dir, "individual"),
                  pattern = paste0("^", s, "_(filtered|seurat|singlets)\\.rds$"),
                  recursive = TRUE, full.names = TRUE)
  if (!length(f)) return(NA_integer_)
  tryCatch(length(Seurat::VariableFeatures(readRDS(f[1]))), error = function(e) NA_integer_)
}, integer(1))
summary$HVG <- hvg_n[summary$Sample]

# ---- marker dot plot data (canonical human-PBMC panel, computed off the object) ----
# Mirrors the families the pipeline draws (T / NK / B / mono / DC / platelet / etc.).
marker_panel <- c("CD3D","CD3E","TRAC","IL7R","CCR7","CD4","CD8A","GZMK","CCL5",
                  "NKG7","GNLY","KLRD1","NCAM1","KLRB1","MS4A1","CD79A","CD19",
                  "CD14","LYZ","S100A8","S100A9","FCGR3A","MS4A7","FCER1A","CLEC9A","CST3",
                  "PPBP","PF4","CD34","CSF3R","GATA2","MS4A2")
marker_panel <- marker_panel[marker_panel %in% rownames(obj)]
markers_dot <- tryCatch({
  d <- Seurat::DotPlot(obj, features = marker_panel, group.by = ct_col)$data
  data.frame(gene = as.character(d$features.plot), cell_type = as.character(d$id),
             avg_scaled = d$avg.exp.scaled, pct_exp = d$pct.exp,
             stringsAsFactors = FALSE)
}, error = function(e) { msg("dotplot skipped: %s", conditionMessage(e)); NULL })
if (!is.null(markers_dot)) {
  attr(markers_dot, "gene_order") <- marker_panel
  msg("dot plot: %d markers x %d cell types", length(marker_panel), length(unique(markers_dot$cell_type)))
}

# ---- differential expression (step 06b), optional ----
de <- NULL
de_path <- file.path(run_dir, "differential", "DE_all_celltypes.csv")
if (file.exists(de_path)) {
  de <- tryCatch(read.csv(de_path, check.names = FALSE), error = function(e) NULL)
  need <- c("avg_log2FC", "p_val_adj", "gene", "cell_type")
  if (!is.null(de) && all(need %in% names(de))) {
    keep_cols <- c(need, intersect(c("pct.1", "pct.2", "comparison"), names(de)))
    de <- de[is.finite(de$avg_log2FC) & is.finite(de$p_val_adj), keep_cols, drop = FALSE]
    de$neglog10p <- -log10(pmax(de$p_val_adj, .Machine$double.xmin))
    de$dir <- ifelse(de$p_val_adj < 0.05 & abs(de$avg_log2FC) >= 1,
                     ifelse(de$avg_log2FC > 0, "up", "down"), "ns")
    # cap per cell type (top 300 by |log2FC|) so the linked widget stays light
    de <- do.call(rbind, lapply(split(de, de$cell_type),
                  function(d) utils::head(d[order(-abs(d$avg_log2FC)), ], 300)))
    de$uid <- paste(de$cell_type, de$gene, seq_len(nrow(de)))  # unique crosstalk key
    rownames(de) <- NULL
    msg("DE: %d rows across %d cell types", nrow(de), length(unique(de$cell_type)))
  } else de <- NULL
}

# ---- static plot galleries (regenerated natively in R; skipped under --lite) ----
galleries <- NULL
if (!lite) {
  msg("building static galleries (pass --lite to skip) ...")
  mk_uri <- function(g, w = 7, h = 4.4, dpi = 96) {
    tf <- tempfile(fileext = ".png")
    ggplot2::ggsave(tf, plot = g, width = w, height = h, dpi = dpi, bg = "white", limitsize = FALSE)
    on.exit(unlink(tf), add = TRUE)
    knitr::image_uri(tf)
  }
  drop_null <- function(x) x[!vapply(x, is.null, logical(1))]

  ## per-sample diagnostics (reached by clicking a sample card)
  by_sample <- list()
  for (s in samples_all) {
    cs <- cells[cells$sample == s, , drop = FALSE]
    imgs <- list()
    imgs[["QC: UMIs vs genes"]] <- tryCatch(mk_uri(
      ggplot(cs, aes(nCount_RNA, nFeature_RNA, color = percent.mt)) +
        geom_point(size = 0.35, alpha = 0.6) + scale_x_log10() + scale_y_log10() +
        scale_color_viridis_c(name = "% MT") + theme_minimal(base_size = 11) +
        labs(x = "UMIs / cell (log)", y = "genes / cell (log)", title = paste(s, "- QC")),
      w = 6.4, h = 4.4), error = function(e) NULL)
    if (!all(is.na(cs$doublet_score)))
      imgs[["Doublet score (UMAP)"]] <- tryCatch(mk_uri(
        ggplot(cs, aes(UMAP_1, UMAP_2, color = doublet_score)) + geom_point(size = 0.4) +
          scale_color_viridis_c(option = "magma", name = "score") + theme_void(base_size = 11) +
          labs(title = paste(s, "- doublet score")),
        w = 5.8, h = 5), error = function(e) NULL)
    hf <- list.files(file.path(run_dir, "individual"),
                     pattern = paste0("^", s, "_(seurat|filtered|singlets)\\.rds$"),
                     recursive = TRUE, full.names = TRUE)
    hf <- hf[order(!grepl("_seurat", hf))]   # prefer the object that has HVG computed
    for (cand in hf) {
      img <- tryCatch(mk_uri(Seurat::VariableFeaturePlot(readRDS(cand)), w = 6.4, h = 4.4),
                      error = function(e) NULL)
      if (!is.null(img)) { imgs[["Highly variable genes"]] <- img; break }
    }
    imgs <- drop_null(imgs)
    if (length(imgs)) by_sample[[s]] <- imgs
  }

  ## run-level galleries (one family per left-nav entry)
  stage <- list()
  fmark <- intersect(c("CD3D","MS4A1","NKG7","CD14","LYZ","FCER1A","PPBP","GNLY"), rownames(obj))
  if (length(fmark)) stage[["Marker feature plots"]] <- tryCatch(
    mk_uri(Seurat::FeaturePlot(obj, features = fmark, reduction = umap_name, order = FALSE,
                               raster = TRUE, raster.dpi = c(200, 200)),
           w = 9, h = 2.7 * ceiling(length(fmark) / 3), dpi = 72), error = function(e) NULL)
  if (!is.null(de)) stage[["Top DE marker heatmap"]] <- tryCatch({
    top <- do.call(rbind, lapply(split(de, de$cell_type),
              function(d) utils::head(d[order(-d$avg_log2FC), ], 5)))
    genes <- unique(top$gene); genes <- genes[genes %in% rownames(obj)]
    ave <- as.matrix(Seurat::AverageExpression(obj, features = genes, group.by = ct_col, assays = "RNA")[[1]])
    z <- t(scale(t(log1p(ave))))
    dfh <- as.data.frame(as.table(z)); colnames(dfh) <- c("gene","cell_type","z")
    dfh$gene <- factor(dfh$gene, levels = rev(genes))
    mk_uri(ggplot(dfh, aes(cell_type, gene, fill = z)) + geom_tile() +
      scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", name = "z") +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank()),
      w = 7.5, h = max(4, 0.17 * length(genes))) }, error = function(e) NULL)
  if (any(c("percent.mt","percent.ribo") %in% colnames(md))) stage[["Contamination (MT / ribo)"]] <- tryCatch({
    parts <- list()
    if ("percent.mt"   %in% colnames(md)) parts[["% MT"]]   <- suppressWarnings(as.numeric(md$percent.mt))
    if ("percent.ribo" %in% colnames(md)) parts[["% ribo"]] <- suppressWarnings(as.numeric(md$percent.ribo))
    ccl <- do.call(rbind, lapply(names(parts), function(nm)
             data.frame(cell_type = md$.cell_type, metric = nm, value = parts[[nm]])))
    mk_uri(ggplot(ccl, aes(cell_type, value, fill = cell_type)) +
      geom_violin(scale = "width", linewidth = 0.2) +
      facet_wrap(~ metric, scales = "free_y", ncol = 1) + theme_minimal(base_size = 10) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank()),
      w = 8, h = 6) }, error = function(e) NULL)
  if ("singler_delta" %in% colnames(md)) stage[["SingleR delta (annotation confidence)"]] <- tryCatch({
    cd <- data.frame(cell_type = md$.cell_type, delta = suppressWarnings(as.numeric(md$singler_delta)))
    cd <- cd[is.finite(cd$delta), , drop = FALSE]
    mk_uri(ggplot(cd, aes(cell_type, delta, fill = cell_type)) +
      geom_violin(scale = "width", linewidth = 0.2) + theme_minimal(base_size = 10) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank()) +
      labs(y = "SingleR delta (higher = more confident)"),
      w = 8, h = 4) }, error = function(e) NULL)
  if ("singler_label_clean" %in% colnames(md)) stage[["SingleR vs final label"]] <- tryCatch({
    tb <- as.data.frame(table(SingleR = as.character(md$singler_label_clean), Final = md$.cell_type))
    tb <- tb[tb$Freq > 0, , drop = FALSE]
    mk_uri(ggplot(tb, aes(Final, SingleR, fill = Freq)) + geom_tile() +
      scale_fill_viridis_c(trans = "log10", name = "cells") + theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)),
      w = 7.6, h = 5.6) }, error = function(e) NULL)
  stage <- drop_null(stage)

  galleries <- list(by_sample = by_sample, stage = stage)
  msg("galleries: %d per-sample sets, %d run-level panels", length(by_sample), length(stage))
}

# ---- bundle + render ----
bundle <- list(
  run_name = run_name, generated = as.character(Sys.time()),
  samples = samples_all, n_cells_total = nrow(cells),
  cells = cells_plot, prop_long = prop_long, prop_wide = prop_wide,
  delta = delta, summary = summary, umap_name = umap_name, de = de,
  markers_dot = markers_dot,
  galleries = galleries, lite = lite,
  umap_static = umap_static, umap_panels = umap_panels
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
