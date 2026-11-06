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

# CELLTYPE_COLORS and QC come from config.R, never local copies: the copies drifted (new cell
# types fell back to arbitrary colours; species QC overlays such as cm never reached the QC
# scatters). Pin config to this run dir so it doesn't derive a default H1/H2 one. Species comes
# from SCRNA_SPECIES (run_pipeline.sh exports it; set it by hand for build_report.sh on a non-human run).
Sys.setenv(SCRNA_RESULTS_DIR = run_dir)
source(file.path(.pipeline_dir, "config.R"))

max_cells   <- as.integer(getflag("max-cells", "6000"))
want_samples <- getflag("samples", NULL)
if (!is.null(want_samples)) want_samples <- trimws(strsplit(want_samples, ",")[[1]])
lite <- "--lite" %in% flags
# --customer: the customer-facing deliverable. Implies --lite (no screenshot
# galleries) and adds the value sections computed below. Default off, so the
# internal report is unaffected.
customer <- "--customer" %in% flags
if (customer) lite <- TRUE

# Long run names (every sample dash-joined) make a >300-char path that Windows
# cannot open over a share (MAX_PATH = 260). Shorten the default report filename to
# <firstSample>_<N>samples when the run name is long; an explicit output path wins.
short_report_name <- function(name) {
  toks <- strsplit(name, "-", fixed = TRUE)[[1]]
  if (length(toks) <= 4) return(name)        # abbreviate once a run has >4 samples
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
  doublet_class = { dc <- col(c("doublet_class","scDblFinder.class"))
                    if (is.na(dc)) NA_character_ else as.character(md[[dc]]) },
  stringsAsFactors = FALSE
)

# ---- frozen-reference per-cell labels: make them the PRIMARY view when present ----
# Additive: keep the de-novo call as cell_type_denovo; the report toggles between them.
# cell_type_ref is run-INDEPENDENT (label transfer onto a frozen reference) and is the
# correct view for cross-species runs where de-novo SingleR mislabels (e.g. bat neutrophils
# collapsing into CD14+ Mono). See 05r_reference_transfer.R / docs/frozen_reference_scope.md.
cells$cell_type_denovo <- cells$cell_type
label_source <- "denovo"
ref_model <- NA_character_
ref_path <- file.path(run_dir, "annotation", "reference_transfer_cells.csv.gz")
if (file.exists(ref_path)) {
  rt <- tryCatch(read.csv(gzfile(ref_path), stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(rt) && all(c("barcode", "cell_type_ref") %in% names(rt))) {
    refmap <- setNames(rt$cell_type_ref, rt$barcode)
    rv <- unname(refmap[rownames(md)])                       # aligned to cells rows
    hit <- sum(!is.na(rv))
    if (hit > 0.5 * nrow(md)) {
      rv[is.na(rv) | rv == ""] <- "Unassigned"
      cells$cell_type_ref <- rv
      cells$cell_type     <- rv                              # PRIMARY = frozen reference
      label_source        <- "ref"
      # mirror onto the Seurat object so the marker dot plot / DE heatmap group by ref too
      ov <- refmap[colnames(obj)]; ov[is.na(ov) | ov == ""] <- "Unassigned"
      obj$cell_type_ref <- unname(ov); ct_col <- "cell_type_ref"
      msg("frozen-ref labels joined (%d/%d cells); main views now use cell_type_ref", hit, nrow(md))
      # runs labelled before 05r wrote its sidecar fall back to the configured model
      mf <- file.path(run_dir, "annotation", "reference_transfer_model.txt")
      ref_model <- if (file.exists(mf)) readLines(mf, n = 1L, warn = FALSE)
                   else if (exists("REFERENCE_MODEL") && nzchar(REFERENCE_MODEL)) REFERENCE_MODEL else NA
    } else msg("frozen-ref barcode overlap too low (%d/%d) — keeping de-novo labels", hit, nrow(md))
  }
}

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

# canonical cell-type colours = config.R CELLTYPE_COLORS (sourced above), so the report matches
# the PDFs; types missing from the palette fall back to a tableau ramp.
PAL_fallback <- c("#4e79a7","#f28e2b","#e15759","#76b7b2","#59a14f","#edc948","#b07aa1",
                  "#ff9da7","#9c755f","#bab0ac","#86bcb6","#d37295","#fabfd2","#8cd17d","#499894")
ctlev <- sort(unique(c(cells_plot$cell_type, cells_plot$cell_type_denovo)))
PAL <- CELLTYPE_COLORS[ctlev]
miss <- is.na(PAL)
if (any(miss)) PAL[miss] <- rep(PAL_fallback, length.out = sum(miss))
names(PAL) <- ctlev
PAL["Unassigned"] <- "#B0B0B0"

# QC acceptance thresholds for the per-sample QC scatters — config.R QC, species overlay applied
QC_THRESH <- QC

# Overview "HVG" count = genes whose vst standardized variance exceeds this (more meaningful
# than reporting the fixed top-2000 selection). The HVG plot draws a dashed line here.
# 1.3 ≈ "clearly above the expected variance" and lands near the real HVG magnitude
# (e.g. ~1460 genes); >5 is far too stringent (only ~17), >1 too inclusive (~6300).
HVG_VAR_THRESHOLD <- 1.3
hvg_count <- function(so, thr = HVG_VAR_THRESHOLD) {
  hvi <- tryCatch(Seurat::HVFInfo(so), error = function(e) NULL)
  if (is.null(hvi)) return(NA_integer_)
  col <- grep("variance\\.standardized$", colnames(hvi), value = TRUE)
  if (!length(col)) return(NA_integer_)
  as.integer(sum(hvi[[col[1]]] > thr, na.rm = TRUE))
}

if (length(samples_all) > 8)
  msg("WARNING: %d samples (>8). UMAP comparison panels render as static images; interactive zoom is available only for runs of 8 or fewer.", length(samples_all))
to_uri <- function(g, w, h, dpi = 150) {
  tf <- tempfile(fileext = ".png")
  # ragg's agg_png antialiases small points better and compresses tighter than the
  # base png device, so we get crisp individual cell dots at modest file size.
  dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggplot2::ggsave(tf, plot = g, width = w, height = h, dpi = dpi,
                  bg = "white", limitsize = FALSE, device = dev)
  on.exit(unlink(tf), add = TRUE); knitr::image_uri(tf)
}

# Past ~8 samples, 14+ live WebGL UMAP widgets blow up the file (16 MB of JSON) and
# exceed the browser's WebGL-context limit, so the report won't open. For large runs
# render the per-sample UMAP panels as light static PNGs instead.
umap_static <- length(samples_all) > 8
umap_panels <- NULL
if (umap_static) {
  n_full <- table(cells$sample)   # true per-sample cell count (not the plotting cap)
  have_repel <- requireNamespace("ggrepel", quietly = TRUE)
  umap_panels <- lapply(samples_all, function(s) {
    dd <- cells_plot[cells_plot$sample == s, , drop = FALSE]
    # cluster labels at each cell type's centroid (median UMAP coords), as in DimPlot(label=TRUE)
    ctr <- aggregate(cbind(UMAP_1, UMAP_2) ~ cell_type, data = dd, FUN = median)
    g <- ggplot(dd, aes(UMAP_1, UMAP_2, color = cell_type)) +
      geom_point(size = 0.35, alpha = 0.85, stroke = 0) +
      scale_color_manual(values = PAL, guide = "none") +
      labs(x = "UMAP 1", y = "UMAP 2",
           title = sprintf("%s  (%s cells)", s, format(as.integer(n_full[s]), big.mark = ","))) +
      theme_classic(base_size = 11) +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            plot.title = element_text(face = "bold", size = 12))
    g <- if (have_repel)
      g + ggrepel::geom_text_repel(data = ctr, aes(UMAP_1, UMAP_2, label = cell_type),
            inherit.aes = FALSE, size = 2.7, fontface = "bold", colour = "black",
            bg.color = "white", bg.r = 0.12, max.overlaps = Inf, seed = 1,
            min.segment.length = 0, segment.size = 0.2, segment.colour = "grey55")
    else
      g + geom_text(data = ctr, aes(UMAP_1, UMAP_2, label = cell_type),
            inherit.aes = FALSE, size = 2.7, fontface = "bold", colour = "black")
    to_uri(g, w = 4.8, h = 4.4, dpi = 150)   # ~720x660 px: enough to resolve single cells
  })
  names(umap_panels) <- samples_all
  msg("static UMAP panels: %d (n_samples > 8)", length(umap_panels))
}

# ---- proportions from the FULL data (not downsampled) ----
tab <- table(cells$cell_type, cells$sample)
prop_wide <- sweep(tab, 2, colSums(tab), "/") * 100   # cell_type x sample, % within sample
prop_long <- as.data.frame(as.table(prop_wide), stringsAsFactors = FALSE)
colnames(prop_long) <- c("cell_type", "sample", "pct")
# same shape but ABSOLUTE counts (for the Cell-numbers stacked bar, cf. Overall_Report.pdf p43)
counts_long <- as.data.frame(as.table(tab), stringsAsFactors = FALSE)
colnames(counts_long) <- c("cell_type", "sample", "n")

# de-novo equivalents (for the label-source toggle) — only differ from the above when ref is primary
tab_dn         <- table(cells$cell_type_denovo, cells$sample)
prop_long_dn   <- as.data.frame(as.table(sweep(tab_dn, 2, colSums(tab_dn), "/") * 100), stringsAsFactors = FALSE)
colnames(prop_long_dn)   <- c("cell_type", "sample", "pct")
counts_long_dn <- as.data.frame(as.table(tab_dn), stringsAsFactors = FALSE)
colnames(counts_long_dn) <- c("cell_type", "sample", "n")

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
qc_thr   <- read_csv_safe(file.path(run_dir, "qc", "qc_thresholds.csv"))   # realised per-sample cut-offs (01)
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
# %MT exactly 0 in every cell means the reference has no mitochondrial genes (the Parse
# split-pipe hg38 reference ships none), so it was not measured: show "—", not a clean 0%.
mt_absent <- vapply(summary$Sample, function(s) {
  v <- cells$percent.mt[cells$sample == s]; length(v) > 0 && all(is.na(v) | v == 0) }, logical(1))
summary$Median_pct_MT[mt_absent] <- NA
summary$Pct_retained    <- g(fate, "Sample", "Pct_retained")
# Retention among barcodes that passed the load gate (>= QC$min_features genes).
# Retention is NOT a quality criterion. A permissive caller's empty droplets fail the
# gene floor and drag retention down while the surviving cells are fine (raw H1D20_2:
# 17.8k called, 8.7k under the 1,000-gene cm floor, 45.7% retained, yet 2,863 median
# genes and rho 0.09 after QC). The low-QC tag uses median genes and the SoupX ambient
# fraction instead; rho is a per-sample quantity read from qc/soupx_rho.csv when present.
.rho_tab <- read_csv_safe(file.path(run_dir, "qc", "soupx_rho.csv"))
if (!is.null(.rho_tab) && all(c("sample", "rho") %in% names(.rho_tab))) {
  .key <- sub("_sx$", "", .rho_tab$sample)
  summary$Ambient_rho <- as.numeric(.rho_tab$rho)[match(sub("_sx$", "", summary$Sample), .key)]
}
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
# highly-variable genes per sample. FindVariableFeatures runs in step 03 on the
# _seurat/_singlets object; the _filtered object pre-dates it (0 HVGs), so prefer the
# computed objects and take the first non-zero count (matches the HVG plot's selection).
hvg_n <- vapply(samples_all, function(s) {
  f <- list.files(file.path(run_dir, "individual"),
                  pattern = paste0("^", s, "_(seurat|singlets|filtered)\\.rds$"),
                  recursive = TRUE, full.names = TRUE)
  if (!length(f)) return(NA_integer_)
  f <- f[order(!grepl("_(seurat|singlets)\\.rds$", f))]   # HVG-bearing objects first
  for (cand in f) {
    n <- tryCatch(hvg_count(readRDS(cand)), error = function(e) NA_integer_)
    if (!is.na(n)) return(n)
  }
  NA_integer_
}, integer(1))
summary$HVG <- hvg_n[summary$Sample]

# ---- marker dot plot data (the species overlay's MARKERS panels, computed off the object) ----
# Panel genes come from config (MARKERS after the species overlay), never from this file,
# so a cardiac or mouse run plots its own lineages rather than a human PBMC set. Up to 3
# genes per panel that are present in the data, in panel order.
marker_panel <- unique(unlist(lapply(Filter(is.character, MARKERS), function(g)
  head(g[g %in% rownames(obj)], 3L)), use.names = FALSE))
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
    if (is.null(de) || !nrow(de)) {
      # every contrast descriptive (NA p-values): nothing to plot, section omitted
      msg("DE: table has no tested contrasts — section omitted"); de <- NULL
    } else {
      de$uid <- paste(de$cell_type, de$gene, seq_len(nrow(de)))  # unique crosstalk key
      rownames(de) <- NULL
      msg("DE: %d rows across %d cell types", nrow(de), length(unique(de$cell_type)))
    }
  } else de <- NULL
}

# ---- frozen-reference benchmark (steps 05r / 08c), optional ----
# Additive: present only when the run was label-transferred onto a frozen reference.
benchmark <- local({
  conc_p <- file.path(run_dir, "benchmark", "concordance.csv")
  sig_p  <- file.path(run_dir, "benchmark", "wholeblood_signature.csv")
  comp_p <- file.path(run_dir, "annotation", "reference_transfer_composition.csv")
  rep_p  <- file.path(run_dir, "benchmark", "benchmark_report.md")
  if (!file.exists(conc_p) && !file.exists(comp_p)) return(NULL)
  rd <- function(p) if (file.exists(p)) tryCatch(read.csv(p, check.names = FALSE), error = function(e) NULL) else NULL
  verdict <- model <- NA_character_
  if (file.exists(rep_p)) {
    ln <- readLines(rep_p, warn = FALSE)
    v <- grep("Verdict", ln, value = TRUE)
    if (length(v)) verdict <- trimws(gsub("\\*\\*|^#+\\s*|Verdict:\\s*", "", v[1]))
    m <- grep("^Model:", ln, value = TRUE)
    if (length(m)) model <- trimws(gsub("Model:\\s*|`", "", m[1]))
  }
  list(concordance = rd(conc_p), signature = rd(sig_p),
       composition = rd(comp_p), verdict = verdict, model = model)
})
if (!is.null(benchmark))
  msg("benchmark: verdict='%s', %d concordance rows", benchmark$verdict,
      if (is.null(benchmark$concordance)) 0L else nrow(benchmark$concordance))

# ---- static plot galleries (regenerated natively in R; skipped under --lite) ----
galleries <- NULL
if (!lite) {
  msg("building static galleries (pass --lite to skip) ...")
  # uniform render target for every per-sample figure. The report CSS caps gallery images at
  # 680 px wide; render a touch larger (~840 px) so they stay crisp without bloating the file.
  GAL_W   <- 5.6    # inches  (5.6 x 150 = 840 px)
  GAL_DPI <- 150
  mk_uri <- function(g, w = 7, h = 4.4, dpi = 96) {
    tf <- tempfile(fileext = ".png")
    ggplot2::ggsave(tf, plot = g, width = w, height = h, dpi = dpi, bg = "white", limitsize = FALSE)
    on.exit(unlink(tf), add = TRUE)
    knitr::image_uri(tf)
  }
  drop_null <- function(x) x[!vapply(x, is.null, logical(1))]
  # rasterise an existing pipeline PDF page (for plots not reconstructable from the rds)
  rasterize_pdf <- function(path, page = 1, dpi = 120) {
    if (!file.exists(path)) return(NULL)
    tf <- tempfile(fileext = ".png")
    ok <- tryCatch({ im <- magick::image_read_pdf(path, pages = page, density = dpi)
                     magick::image_write(im, tf, "png"); TRUE }, error = function(e) FALSE)
    if (!ok) ok <- tryCatch({ bm <- pdftools::pdf_render_page(path, page = page, dpi = dpi)
                     magick::image_write(magick::image_read(bm), tf, "png"); TRUE }, error = function(e) FALSE)
    if (ok) { on.exit(unlink(tf), add = TRUE); knitr::image_uri(tf) } else NULL
  }
  # locate a per-sample doublet PDF. Step 02 computes it ONCE and caches it under
  # sample_cache/<s>/, restoring it into this run's doublets/ on a cache hit. Prefer the
  # run's doublets/, then the sample_cache copy, then (legacy) the newest sibling run.
  results_root <- dirname(run_dir)
  sample_cache <- file.path(dirname(results_root), "sample_cache")
  find_dbl_pdf <- function(s, kind) {
    f <- sprintf("%s_doublet_%s.pdf", s, kind)
    local <- file.path(run_dir, "doublets", f);        if (file.exists(local))  return(local)
    cached <- file.path(sample_cache, s, f);           if (file.exists(cached)) return(cached)
    cand <- Sys.glob(file.path(results_root, "*", "doublets", f))
    if (!length(cand)) return(NA_character_)
    cand[order(file.info(cand)$mtime, decreasing = TRUE)][1]
  }

  ## per-sample diagnostics (reached by clicking a sample card) — match Overall_Report.pdf
  by_sample <- list()
  for (s in samples_all) {
    cs <- cells_plot[cells_plot$sample == s, , drop = FALSE]
    imgs <- list()
    # QC: UMIs vs genes — linear axes from 0, dashed acceptance box (config QC thresholds)
    imgs[["QC: UMIs vs genes"]] <- tryCatch(mk_uri(
      ggplot(cs, aes(nCount_RNA, nFeature_RNA)) +
        annotate("rect", xmin = QC_THRESH$min_counts, xmax = QC_THRESH$max_counts,
                 ymin = QC_THRESH$min_features, ymax = QC_THRESH$max_features,
                 fill = NA, colour = "#E64B35", linetype = "dashed", linewidth = 0.5) +
        geom_point(aes(colour = percent.mt), size = 0.35, alpha = 0.6) +
        scale_colour_viridis_c(name = "% MT") +
        scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .03))) +
        scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .03))) +
        theme_classic(base_size = 11) +
        labs(x = "UMIs / cell", y = "genes / cell", title = paste(s, "- UMIs vs genes")),
      w = GAL_W, h = GAL_W * 0.71, dpi = GAL_DPI), error = function(e) NULL)
    # QC: UMIs vs % MT — axes from 0, dashed %MT threshold
    imgs[["QC: UMIs vs % MT"]] <- tryCatch(mk_uri(
      ggplot(cs, aes(nCount_RNA, percent.mt)) +
        geom_point(size = 0.35, alpha = 0.5, colour = "#4DBBD5") +
        geom_hline(yintercept = QC_THRESH$max_percent_mt, linetype = "dashed", colour = "#E64B35") +
        scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .03))) +
        scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .03))) +
        theme_classic(base_size = 11) +
        labs(x = "UMIs / cell", y = "% MT", title = paste(s, "- UMIs vs % MT")),
      w = GAL_W, h = GAL_W * 0.71, dpi = GAL_DPI), error = function(e) NULL)
    # Doublets — rasterise the step-02 PDFs (the integrated object has doublets ALREADY
    # removed, so recomputing here gives 0/sample; 02_doublets.R drew these pre-removal
    # with the true counts). Samples served from sample_cache/ were not re-run here, so
    # their PDFs live in the run that first processed them — fall back across Results/*.
    # 02_doublets.R saves these at 7x6 in (umap) and 6x4 in (hist); pick a dpi that yields
    # the same target pixel width as the ggplot panels so every per-sample figure matches.
    dbl_umap <- find_dbl_pdf(s, "umap")
    if (!is.na(dbl_umap)) imgs[["Doublets (UMAP)"]] <- rasterize_pdf(dbl_umap, dpi = round(GAL_W * GAL_DPI / 7))
    dbl_hist <- find_dbl_pdf(s, "score_hist")
    if (!is.na(dbl_hist)) imgs[["Doublet score distribution"]] <- rasterize_pdf(dbl_hist, dpi = round(GAL_W * GAL_DPI / 6))
    # HVG — VariableFeaturePlot + the top-10 labelled dots (03_individual.R)
    hf <- list.files(file.path(run_dir, "individual"),
                     pattern = paste0("^", s, "_(seurat|filtered|singlets)\\.rds$"),
                     recursive = TRUE, full.names = TRUE)
    hf <- hf[order(!grepl("_seurat", hf))]   # prefer the object that has HVG computed
    for (cand in hf) {
      img <- tryCatch({
        so <- readRDS(cand)
        hvf <- Seurat::VariableFeatures(so)
        if (!length(hvf)) NULL else {       # skip objects without HVG (e.g. _filtered)
          ph <- Seurat::VariableFeaturePlot(so)
          ph <- Seurat::LabelPoints(plot = ph, points = utils::head(hvf, 10),
                                    repel = TRUE, xnudge = 0, ynudge = 0, size = 3)
          # the default legend sits on the right and eats ~half the panel; move it under the
          # plot and shrink it so the scatter owns the space (aspect matched to the doublet UMAP)
          ph <- ph +
            ggplot2::geom_hline(yintercept = HVG_VAR_THRESHOLD, linetype = "dashed",
                                colour = "#888888", linewidth = 0.4) +
            ggplot2::theme(
              legend.position = "bottom", legend.direction = "horizontal",
              legend.title = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = 8),
              legend.key.size = grid::unit(9, "pt"), legend.margin = ggplot2::margin(t = -2)) +
            ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2)))
          mk_uri(ph + labs(title = paste(s, "- Variable genes (top 10 labelled)")),
                 w = GAL_W, h = GAL_W * 0.70, dpi = GAL_DPI)
        }
      }, error = function(e) NULL)
      if (!is.null(img)) { imgs[["Highly variable genes"]] <- img; break }
    }
    imgs <- drop_null(imgs)
    # 2-up gallery rows: QC pair, then the two doublet panels together, then HVG. Figures use a
    # one-line caption above the image (description below), so image tops are flush within a row.
    ord <- c("QC: UMIs vs genes", "QC: UMIs vs % MT",
             "Doublets (UMAP)", "Doublet score distribution", "Highly variable genes")
    imgs <- imgs[c(intersect(ord, names(imgs)), setdiff(names(imgs), ord))]
    if (length(imgs)) by_sample[[s]] <- imgs
  }

  ## run-level galleries (one family per left-nav entry)
  stage <- list()
  # Top DE marker heatmap (kept) — top genes per type x scaled mean expression
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
  # Key lineage marker expression (replaces SingleR delta) — 05_annotate.R FeaturePlot style
  # lead marker of each MARKERS panel present in the data, so the species overlay sets it
  lin <- head(unique(unlist(lapply(Filter(is.character, MARKERS), function(g)
    head(g[g %in% rownames(obj)], 1L)), use.names = FALSE)), 16L)
  if (length(lin)) stage[["Key lineage marker expression"]] <- tryCatch(
    mk_uri(Seurat::FeaturePlot(obj, features = lin, reduction = umap_name, ncol = 4, order = TRUE,
                               cols = c("lightgrey", "#E64B35"), raster = TRUE, raster.dpi = c(200, 200)) &
             theme(plot.title = element_text(size = 10), axis.text = element_blank(), axis.ticks = element_blank()),
           w = 10, h = 2.5 * ceiling(length(lin) / 4), dpi = 72), error = function(e) NULL)
  # SingleR reference score heatmap (replaces SingleR-vs-final) — the raw SingleR result is
  # not stored in the rds, so rasterise the pipeline's annotation/singler_scores_heatmap.pdf.
  stage[["SingleR reference heatmap"]] <-
    rasterize_pdf(file.path(run_dir, "annotation", "singler_scores_heatmap.pdf"), dpi = 130)
  stage <- drop_null(stage)

  galleries <- list(by_sample = by_sample, stage = stage)
  msg("galleries: %d per-sample sets, %d run-level panels", length(by_sample), length(stage))
}

# ---- scRNA analysis tool versions (only the single-cell stack; report-gen tools excluded) ----
# Clustering algorithm as the run actually used it (Seurat logs FindClusters' arguments),
# not the current config: 1 = Louvain, 4 = Leiden.
fc_cmd   <- grep("^FindClusters", names(obj@commands), value = TRUE)
clus_alg <- if (length(fc_cmd)) obj@commands[[tail(fc_cmd, 1)]]$algorithm else NULL
clus_lbl <- if (identical(as.numeric(clus_alg), 4)) "Leiden" else if (identical(as.numeric(clus_alg), 1)) "Louvain" else NULL
tool_roles <- c(
  Seurat               = sprintf("Core single-cell toolkit (QC, normalization, %sclustering, UMAP)",
                                 if (is.null(clus_lbl)) "" else paste(clus_lbl, "")),
  SeuratObject         = "Single-cell data structures for Seurat",
  harmony              = "Batch integration across samples",
  scDblFinder          = "Doublet detection",
  SingleR              = "Reference-based cell-type annotation",
  celldex              = "Reference expression datasets for SingleR",
  presto               = "Fast Wilcoxon tests for cluster marker genes",
  uwot                 = "UMAP embedding",
  SingleCellExperiment = "Bioconductor single-cell data container",
  # listed because the customer report's Ambient RNA section reports its rho;
  # the row self-omits when the package is not installed
  SoupX                = "Ambient (free-floating) RNA estimation and correction",
  scattermore          = "Fast rasterised scatter rendering for large cell counts")
if (identical(clus_lbl, "Leiden"))
  tool_roles <- c(tool_roles, leidenbase = "Leiden community detection for clustering",
                  igraph = "Graph backend for Leiden clustering")
# step 06b (edgeR) only runs with 2+ conditions; list it only when it left results
if (length(list.files(file.path(run_dir, "differential"), pattern = "\\.csv$", recursive = TRUE)))
  tool_roles <- c(tool_roles, edgeR = "Pseudobulk differential expression between conditions")
tool_urls <- c(
  R                    = "https://www.r-project.org",
  Seurat               = "https://satijalab.org/seurat",
  SeuratObject         = "https://github.com/satijalab/seurat-object",
  harmony              = "https://github.com/immunogenomics/harmony",
  scDblFinder          = "https://bioconductor.org/packages/scDblFinder",
  SingleR              = "https://bioconductor.org/packages/SingleR",
  celldex              = "https://bioconductor.org/packages/celldex",
  SingleCellExperiment = "https://bioconductor.org/packages/SingleCellExperiment",
  SoupX                = "https://github.com/constantAmateur/SoupX",
  scattermore          = "https://github.com/exaexa/scattermore",
  presto               = "https://github.com/immunogenomics/presto",
  uwot                 = "https://github.com/jlmelville/uwot",
  edgeR                = "https://bioconductor.org/packages/edgeR",
  leidenbase           = "https://github.com/cole-trapnell-lab/leidenbase",
  igraph               = "https://r.igraph.org")
tv_rows <- lapply(names(tool_roles), function(p) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  if (is.na(v)) NULL else data.frame(Tool = p, Version = v, Role = unname(tool_roles[p]),
                                     Source = unname(tool_urls[p]), stringsAsFactors = FALSE)
})
tool_versions <- do.call(rbind, c(
  list(data.frame(Tool = "R", Version = paste(R.version$major, R.version$minor, sep = "."),
                  Role = "Statistical computing environment", Source = unname(tool_urls["R"]),
                  stringsAsFactors = FALSE)),
  tv_rows))
msg("tool versions: %s", paste(tool_versions$Tool, tool_versions$Version, collapse = ", "))

# ---- customer-report extras (--customer; implies --lite) -------------------
# Additions that carry real value for a customer with no bioinformatics staff,
# all rendered from data — never a rasterised PDF page. Each is NULL unless
# --customer is set, and every template chunk is guarded on non-NULL, so the
# internal report is byte-identical when this flag is absent.
cust <- NULL
if (customer) {
  msg("building customer extras ...")
  md <- obj@meta.data
  ct <- as.character(md[[if ("cell_type_ref" %in% colnames(md) &&
                            !all(is.na(md$cell_type_ref))) "cell_type_ref" else "cell_type"]])
  smp <- as.character(md$sample)
  cmat <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")

  # (1) per-sample QC detail: ONE raster per sample, so the report can drive them
  # with the same sample-chip controls the UMAP section already uses. Rasterised
  # (scattermore) so every cell is drawn for ~40 KB each rather than ~350 KB of
  # widget JSON per sample.
  qc_df <- data.frame(sample = smp, nCount = md$nCount_RNA, nFeature = md$nFeature_RNA,
                      mt = if ("percent.mt" %in% colnames(md)) md$percent.mt else NA_real_)
  mt_cap <- if (all(is.na(qc_df$mt))) NA_real_ else
    max(1, as.numeric(quantile(qc_df$mt, 0.99, na.rm = TRUE)))
  qc_png <- setNames(vector("list", length(samples_all)), samples_all)
  for (s in samples_all) {
    d <- qc_df[qc_df$sample == s, , drop = FALSE]
    if (!nrow(d)) next
    g <- ggplot2::ggplot(d, ggplot2::aes(nCount, nFeature, colour = mt)) +
      scattermore::geom_scattermore(pointsize = 2.4) +
      ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
      # shared limits across samples so panels are directly comparable
      ggplot2::coord_cartesian(xlim = range(qc_df$nCount), ylim = range(qc_df$nFeature)) +
      ggplot2::labs(x = "UMIs per cell", y = "Genes per cell") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(legend.position = "none",
                     panel.grid.minor = ggplot2::element_blank())
    g <- if (is.na(mt_cap)) g + ggplot2::scale_colour_continuous(guide = "none") else
      g + ggplot2::scale_colour_viridis_c(option = "magma", direction = -1,
                                          limits = c(0, mt_cap), oob = scales::squish)
    tf <- file.path(tempdir(), paste0("cust_qc_", gsub("[^A-Za-z0-9]", "_", s), ".png"))
    # sized for a 3-4 column small-multiples grid in the report
    ok <- tryCatch({ ggplot2::ggsave(tf, g, width = 3.4, height = 2.8, dpi = 150,
                                     bg = "white"); TRUE }, error = function(e) FALSE)
    if (ok) qc_png[[s]] <- knitr::image_uri(tf)
  }
  qc_png <- qc_png[!vapply(qc_png, is.null, logical(1))]
  if (!length(qc_png)) qc_png <- NULL
  # the colour scale is described to the template as a range plus the viridis
  # stops, so it renders as a CSS gradient — no extra image, no new dependency
  qc_scale <- if (is.na(mt_cap)) NULL else
    list(max = round(mt_cap, 1),
         stops = scales::viridis_pal(option = "magma", direction = -1)(5))

  # (2) pseudobulk sample-sample correlation. Catches a swapped or mislabelled
  # library, and shows whether replicates actually agree.
  pb <- sapply(split(seq_len(ncol(cmat)), smp), function(i)
          Matrix::rowSums(cmat[, i, drop = FALSE]))
  cpm  <- t(t(pb) / colSums(pb)) * 1e6
  keep <- rowMeans(cpm) > 1
  cor_pb <- round(cor(log1p(cpm[keep, , drop = FALSE]), method = "spearman"), 4)

  # (3) per-CLUSTER curation evidence. Deliberately NOT a cell-type x marker
  # detection matrix: the existing dot plot already encodes % expressing as dot
  # size, so that would be duplicate. What is missing from the report is which
  # clusters were merged into each label and the markers that decided it.
  panel <- unique(unlist(if (exists("ALL_MARKERS") && length(ALL_MARKERS)) ALL_MARKERS
                         else if (exists("MARKERS")) MARKERS else list()))
  panel <- intersect(panel, rownames(cmat))
  evidence <- NULL
  if (length(panel) && "seurat_clusters" %in% colnames(md)) {
    clus <- as.character(md$seurat_clusters)
    sub  <- cmat[panel, , drop = FALSE]
    ord  <- order(suppressWarnings(as.integer(unique(clus))), unique(clus))
    evidence <- do.call(rbind, lapply(unique(clus)[ord], function(k) {
      i   <- which(clus == k)
      pc  <- round(100 * Matrix::rowSums(sub[, i, drop = FALSE] > 0) / length(i), 0)
      top <- head(sort(pc, decreasing = TRUE), 5L)
      lab <- names(sort(table(ct[i]), decreasing = TRUE))[1]
      sl  <- if ("singler_label_clean" %in% colnames(md))
               names(sort(table(as.character(md$singler_label_clean[i])), decreasing = TRUE))[1]
             else NA_character_
      data.frame(Cluster = k, Cells = length(i), `Cell type` = lab,
                 `Deciding markers (% of cluster's cells)` =
                   paste(sprintf("%s %d%%", names(top), as.integer(top)), collapse = ", "),
                 `Reference call` = if (is.null(sl)) NA_character_ else sl,
                 check.names = FALSE, stringsAsFactors = FALSE)
    }))
  }

  # (4) ambient RNA per sample. Read ONLY from the run directory — never from a
  # scratch path, so a report can never silently carry another cohort's numbers.
  # The SoupX step is expected to write qc/soupx_rho.csv with sample + rho columns.
  rho <- NULL
  rho_csv <- file.path(run_dir, "qc", "soupx_rho.csv")
  if (file.exists(rho_csv)) {
    r <- tryCatch(read.csv(rho_csv, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(r) && all(c("sample", "rho") %in% names(r))) {
      # keep only rows that belong to THIS cohort, matching with and without _sx
      key <- sub("_sx$", "", samples_all)
      r <- r[sub("_sx$", "", r$sample) %in% key, , drop = FALSE]
      if (nrow(r)) rho <- r
    }
    if (is.null(rho)) msg("  ambient: %s present but no rows match this cohort — section omitted",
                          basename(rho_csv))
  }

  # (5) target-population panel: the species overlay's TARGET_POPULATION names the
  # population, its marker + UMI floor, subtypes and expression ratios. No spec, no panel.
  eff <- eff_spec <- NULL; tp <-if (exists("TARGET_POPULATION")) TARGET_POPULATION else NULL
  if (!is.null(tp) && tp$marker %in% rownames(cmat)) {
    is_t <- grepl(tp$label, ct)
    ratios <- Filter(function(g) all(g %in% rownames(cmat)), tp$ratios)
    mcol <- sprintf("%s+ cells (%%)", tp$marker)
    eff <- do.call(rbind, lapply(samples_all, function(s) {
      i <- which(smp == s); j <- i[is_t[i]]
      row <- list(Sample = s, Cells = length(i),
                  round(100 * mean(cmat[tp$marker, i] >= tp$min_umi), 1),
                  round(100 * length(j) / length(i), 1))
      names(row)[3:4] <- c(mcol, sprintf("Labelled %s (%%)", tp$name))
      for (k in names(tp$subtypes))  # a curation without subtype labels is dropped below
        row[[sprintf("%s (%% of %s)", k, tp$short)]] <-
          if (length(j)) round(100 * mean(grepl(tp$subtypes[[k]], ct[j])), 1) else NA_real_
      for (k in names(ratios)) { g <- ratios[[k]]
        row[[k]] <- if (length(j)) round(sum(cmat[g[1], j]) / max(1, sum(cmat[g[2], j])), 2) else NA_real_ }
      as.data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
    }))
    for (k in names(tp$subtypes)) if (!any(grepl(tp$subtypes[[k]], ct)))
      eff[[sprintf("%s (%% of %s)", k, tp$short)]] <- NULL
    eff_spec <- list(title = tp$title, desc = tp$desc, marker_col = mcol)
  }

  # (5b) unexpected types: labels outside the overlay's EXPECTED_TYPES (prefix match) at
  # >= UNEXPECTED_FLAG_PCT % of a sample. A 0-row table is reported as "none found".
  unexp <- NULL; thr <- if (exists("UNEXPECTED_FLAG_PCT")) UNEXPECTED_FLAG_PCT else 1
  if (exists("EXPECTED_TYPES") && length(EXPECTED_TYPES)) {
    u <- unique(ct)
    bad <- u[!vapply(u, function(l) any(startsWith(l, c(EXPECTED_TYPES, "Unknown", "Unassigned"))), TRUE)]
    unexp <- do.call(rbind, c(list(data.frame(Sample = character(0), `Cell type` = character(0),
                                              Cells = integer(0), `% of sample` = numeric(0),
                                              check.names = FALSE)),
      lapply(samples_all, function(s) {
        v <- ct[smp == s]; n <- table(factor(v[v %in% bad], levels = bad)); p <- 100 * n / length(v)
        k <- p >= thr
        if (!any(k)) return(NULL)
        data.frame(Sample = s, `Cell type` = names(n)[k], Cells = as.integer(n[k]),
                   `% of sample` = round(as.numeric(p[k]), 1), check.names = FALSE)
      })))
  }

  # (6) composition statistics (step 09): replicate-aware proportion test + per-sample CIs
  csv_or_null <- function(p) if (file.exists(p)) tryCatch(read.csv(p, stringsAsFactors = FALSE),
                                                           error = function(e) NULL) else NULL
  # Optional <run>/report_contrasts.txt, one contrast per line as named in the step 06b/09
  # outputs (e.g. 24hMed_vs_24hCell). Steps 06b/09 test every pair of groups; on a design
  # where some pairs differ in two factors at once (arm AND time, or a flowcell boundary)
  # those pairs are confounded, so the customer sees only the listed ones. Absent = all.
  ctr_f <- file.path(run_dir, "report_contrasts.txt")
  contrasts <- if (file.exists(ctr_f)) { v <- trimws(readLines(ctr_f, warn = FALSE))
                                         v[nzchar(v) & !startsWith(v, "#")] } else NULL
  comp <- csv_or_null(file.path(run_dir, "composition_test.csv"))
  if (!is.null(comp) && !is.null(contrasts)) comp <- comp[comp$contrast %in% contrasts, , drop = FALSE]
  if (!is.null(comp) && !nrow(comp)) comp <- NULL
  boot <- csv_or_null(file.path(run_dir, "bootstrap_summary.csv"))

  # (7) timepoint DE volcano (step 06b, pseudobulk edgeR QL). Contrasts without replicates on
  # both sides have no p-values; they are named, never plotted. Per panel: every significant
  # gene up to the 300 strongest, plus 300 random non-significant genes as the backdrop.
  de_vol <- de_top <- NULL; de_desc <- character(0)
  dd <- if (file.exists(de_path)) tryCatch(read.csv(de_path, stringsAsFactors = FALSE),
                                           error = function(e) NULL) else NULL
  if (!is.null(dd) && !all(c("gene", "cell_type", "comparison", "avg_log2FC", "p_val_adj") %in% names(dd))) dd <- NULL
  if (!is.null(dd)) {
    if (!is.null(contrasts)) dd <- dd[dd$comparison %in% contrasts, , drop = FALSE]
    inf <- is.finite(dd$p_val_adj) & (if ("inferential" %in% names(dd)) as.logical(dd$inferential) else TRUE)
    de_desc <- setdiff(unique(dd$comparison), unique(dd$comparison[inf]))
    dd <- dd[inf & is.finite(dd$avg_log2FC), ]
  }
  if (!is.null(dd) && nrow(dd)) {   # nothing left once filtered = section omitted
    dd$nl10 <- round(pmin(-log10(pmax(dd$p_val_adj, 1e-300)), 300), 2)
    dd$sig  <- dd$p_val_adj < 0.05 & abs(dd$avg_log2FC) >= 1
    dd$avg_log2FC <- round(dd$avg_log2FC, 2)
    set.seed(1)
    parts <- split(dd, paste(dd$cell_type, dd$comparison, sep = "|"))
    de_vol <- do.call(rbind, lapply(parts, function(x) {
      s <- x[x$sig, ]; s <- head(s[order(s$p_val_adj, -abs(s$avg_log2FC)), ], 300)
      b <- x[!x$sig, ]; b <- b[sample.int(nrow(b), min(nrow(b), 300)), ]
      rbind(s, b)[, c("cell_type", "comparison", "gene", "avg_log2FC", "nl10", "sig")]
    }))
    de_n <- do.call(rbind, lapply(parts, function(x)
      data.frame(cell_type = x$cell_type[1], comparison = x$comparison[1],
                 up = sum(x$sig & x$avg_log2FC > 0), down = sum(x$sig & x$avg_log2FC < 0))))
    attr(de_vol, "n") <- de_n
    de_top <- do.call(rbind, lapply(parts, function(x) {
      s <- x[x$sig, ]; head(s[order(s$p_val_adj), ], 50) }))
    if (!is.null(de_top) && nrow(de_top))
      de_top <- data.frame(Gene = de_top$gene, `Cell type` = de_top$cell_type,
                           Contrast = de_top$comparison, log2FC = de_top$avg_log2FC,
                           FDR = signif(de_top$p_val_adj, 3), check.names = FALSE)
    rownames(de_vol) <- NULL; rm(dd, parts)
  }

  # (8) the run's Excel workbook (08f + build_export_xlsx.py), embedded as a download
  xl <- Sys.glob(file.path(run_dir, "export", "*_analysis.xlsx"))
  xlsx <- if (length(xl)) list(name = basename(xl[1]), size = file.size(xl[1]),
                               uri = xfun::base64_uri(xl[1])) else NULL

  # (9) per-cell table (barcode, sample, cluster, label, UMAP) as a gzip CSV download.
  # UMAP is rounded to 3 decimals: plenty for plotting, ~25% smaller on 100k+ cells.
  percell <- tryCatch({
    pm <- obj@meta.data
    keep <- as.character(if (!is.na(sample_col)) pm[[sample_col]] else pm$sample) %in% samples_all
    pm <- pm[keep, , drop = FALSE]
    pc <- data.frame(barcode = rownames(pm),
                     sample  = as.character(pm[[sample_col]]),
                     cluster = if ("seurat_clusters" %in% colnames(pm)) as.character(pm$seurat_clusters)
                               else NA_character_,
                     stringsAsFactors = FALSE)
    pc$cell_type <- as.character(pm[[if (label_source == "ref" && "cell_type_ref" %in% colnames(pm)) "cell_type_ref"
                                     else ct_col]])
    if ("cell_type_ref" %in% colnames(pm) && !all(is.na(pm$cell_type_ref)))
      pc$cell_type_ref <- as.character(pm$cell_type_ref)
    if ("cell_type_denovo" %in% colnames(pm) && !all(is.na(pm$cell_type_denovo)))
      pc$cell_type_denovo <- as.character(pm$cell_type_denovo)
    pc$UMAP_1 <- round(emb[rownames(pm), "UMAP_1"], 3)
    pc$UMAP_2 <- round(emb[rownames(pm), "UMAP_2"], 3)
    if (all(is.na(pc$cluster))) pc$cluster <- NULL
    tf <- tempfile(fileext = ".csv.gz")
    gz <- gzfile(tf, "w", compression = 9); utils::write.csv(pc, gz, row.names = FALSE); close(gz)
    list(name = paste0(run_name, "_per_cell.csv.gz"), size = file.size(tf), n = nrow(pc),
         cols = colnames(pc), uri = xfun::base64_uri(tf))
  }, error = function(e) { msg("per-cell table skipped: %s", conditionMessage(e)); NULL })
  if (!is.null(percell)) msg("per-cell table: %d rows, %.1f MB gz", percell$n, percell$size / 1e6)

  cust <- list(qc_png = qc_png, qc_scale = qc_scale, cor_pb = cor_pb,
               evidence = evidence, rho = rho, n_panel = length(panel),
               eff = eff, eff_spec = eff_spec, unexp = unexp, unexp_thr = thr,
               expected = if (exists("EXPECTED_TYPES")) EXPECTED_TYPES, comp = comp, boot = boot, de_vol = de_vol, de_top = de_top,
               de_desc = de_desc, xlsx = xlsx, percell = percell, contrasts = contrasts)
  msg("customer extras: qc panels=%d cor=%dx%d evidence=%s rho=%s",
      if (is.null(qc_png)) 0L else length(qc_png), nrow(cor_pb), ncol(cor_pb),
      if (is.null(evidence)) "none" else paste0(nrow(evidence), " clusters"),
      if (is.null(rho)) "none" else paste0(nrow(rho), " samples"))
  rm(cmat); gc(verbose = FALSE)
}

# ---- bundle + render ----
bundle <- list(
  customer = customer, cust = cust,
  run_name = run_name, generated = as.character(Sys.time()),
  samples = samples_all, n_cells_total = nrow(cells),
  cells = cells_plot, prop_long = prop_long, prop_wide = prop_wide, counts_long = counts_long,
  label_source = label_source, ref_model = ref_model,
  singler_ref = if (exists("SINGLER_REF")) SINGLER_REF else NA,
  species = Sys.getenv("SCRNA_SPECIES", "human"), prop_long_denovo = prop_long_dn, counts_long_denovo = counts_long_dn,
  delta = delta, summary = summary, qc_thresholds = qc_thr, umap_name = umap_name, de = de,
  markers_dot = markers_dot,
  galleries = galleries, lite = lite, benchmark = benchmark,
  umap_static = umap_static, umap_panels = umap_panels,
  cell_colors = PAL, tool_versions = tool_versions,
  hvg_var_threshold = HVG_VAR_THRESHOLD
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
