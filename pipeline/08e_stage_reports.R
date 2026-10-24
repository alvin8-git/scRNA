# =============================================================================
# 08e_stage_reports.R — independent self-contained HTML per pipeline stage
#
# Companions to (not replacements for) the stage PDFs that 07_finalize_reports.R
# writes. 07 concatenates existing PDFs with .combine_pdfs(), so a stage PDF is a
# *collection* of composed pages; these HTML reports instead render every figure
# from the run's own data and CSVs, so nothing is a rasterised PDF page.
#
# Writes into <run>/reports/:
#   01-QC.html  02-Doublet.html  03-Individual.html  04-Annotation.html
#   05-Integrated.html  index.html
# index.html is the Overall equivalent: Overall_report.pdf is every stage
# concatenated, and the interactive run report already covers that in HTML, so
# duplicating it would add a third copy of the same content.
#
# Usage: Rscript pipeline/08e_stage_reports.R <run_dir> [stage ...]
#        stages default to all; e.g. ... <run_dir> 01 04
# =============================================================================
.pipeline_dir <- local({
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(f)) dirname(f) else {
    a <- commandArgs(trailingOnly = FALSE)
    d <- sub("--file=", "", a[grep("--file=", a)])
    if (length(d) > 0) dirname(normalizePath(d)) else "."
  }
})
args    <- commandArgs(trailingOnly = TRUE)
pos     <- args[!grepl("^--", args)]
run_dir <- if (length(pos)) pos[1] else stop("usage: 08e_stage_reports.R <run_dir> [stage ...]")
want    <- if (length(pos) > 1) pos[-1] else c("01", "02", "03", "04", "05")

Sys.setenv(SCRNA_RESULTS_DIR = run_dir)
source(file.path(.pipeline_dir, "config.R"))
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(ggplot2) })

msg <- function(fmt, ...) message(sprintf(paste0("[08e] ", fmt), ...))
rp  <- function(...) file.path(run_dir, ...)
out_dir <- rp("reports"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
TEMPLATE <- file.path(.pipeline_dir, "stage_report.Rmd")
stopifnot(file.exists(TEMPLATE))

# pandoc is not on PATH under a bare Rscript call (same shim as 08b:576)
if (!rmarkdown::pandoc_available()) {
  env_bin <- file.path(dirname(dirname(Sys.getenv("R_HOME"))), "bin")
  if (file.exists(file.path(env_bin, "pandoc"))) {
    Sys.setenv(RSTUDIO_PANDOC = env_bin); rmarkdown::find_pandoc(dir = env_bin)
  }
}
msg("pandoc: %s", tryCatch(as.character(rmarkdown::pandoc_version()),
                           error = function(e) "NOT FOUND"))

csv <- function(p) if (file.exists(p)) tryCatch(read.csv(p, stringsAsFactors = FALSE),
                                                error = function(e) NULL) else NULL
uri <- function(g, w = 7, h = 4.4, dpi = 140) {
  tf <- tempfile(fileext = ".png")
  ok <- tryCatch({ ggsave(tf, g, width = w, height = h, dpi = dpi, bg = "white"); TRUE },
                 error = function(e) FALSE)
  if (ok) knitr::image_uri(tf) else NULL
}

# ---- shared: the integrated object drives 01, 02, 04, 05 -------------------
need_obj <- any(want %in% c("01", "02", "04", "05"))
obj <- NULL; md <- NULL; PAL <- NULL
if (need_obj) {
  f <- rp("integrated", "integrated_annotated.rds")
  if (!file.exists(f)) stop("missing ", f, " — run step 05 first")
  msg("reading %s", f)
  obj <- readRDS(f)
  md  <- obj@meta.data
  lab <- if ("cell_type_ref" %in% colnames(md) && !all(is.na(md$cell_type_ref)))
           "cell_type_ref" else "cell_type"
  md$.ct <- as.character(md[[lab]]); md$.sample <- as.character(md$sample)
  sl <- if (exists("SAMPLE_NAMES")) c(intersect(SAMPLE_NAMES, unique(md$.sample)),
                                      sort(setdiff(unique(md$.sample), SAMPLE_NAMES)))
        else sort(unique(md$.sample))
  md$.sample <- factor(md$.sample, levels = sl)
  u <- sort(unique(md$.ct))
  known <- if (exists("CELLTYPE_COLORS")) CELLTYPE_COLORS[u] else rep(NA, length(u))
  PAL <- setNames(ifelse(is.na(known), scales::hue_pal()(length(u)), known), u)
  msg("%s cells | %d samples | %d types | labels from %s",
      format(nrow(md), big.mark = ","), length(sl), length(u), lab)
}

mt_col <- if (!is.null(md)) intersect(c("percent.mt", "percent_mt"), colnames(md)) else character(0)

render_stage <- function(id, title, blurb, bundle, outfile) {
  bp <- file.path(tempdir(), paste0("stage_", id, ".rds"))
  saveRDS(c(bundle, list(.id = id, .title = title, .blurb = blurb,
                         .run = basename(run_dir),
                         .generated = format(Sys.time(), "%Y-%m-%d %H:%M"))), bp)
  rmarkdown::render(TEMPLATE, output_file = basename(outfile), output_dir = dirname(outfile),
                    params = list(bundle = bp, title = title), quiet = TRUE,
                    envir = new.env())
  msg("%-18s %5.2f MB", basename(outfile), file.info(outfile)$size / 1e6)
  unlink(bp)
}

# ============================ 01 — QC =======================================
if ("01" %in% want) {
  long <- do.call(rbind, lapply(intersect(c("nFeature_RNA", "nCount_RNA", mt_col),
                                          colnames(md)), function(m)
    data.frame(sample = md$.sample, metric = m, value = md[[m]])))
  long$metric <- factor(long$metric,
    levels = c("nFeature_RNA", "nCount_RNA", mt_col),
    labels = c("Genes / cell", "UMIs / cell", "% mitochondrial")[seq_along(c("nFeature_RNA","nCount_RNA",mt_col))])
  gv <- ggplot(long, aes(sample, value, fill = sample)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
    facet_wrap(~ metric, scales = "free_y") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "none",
          axis.title = element_blank(), strip.text = element_text(face = "bold"))
  # per-sample UMI-vs-gene rasters, every cell drawn
  qcp <- list()
  for (s in levels(md$.sample)) {
    d <- data.frame(x = md$nCount_RNA[md$.sample == s], y = md$nFeature_RNA[md$.sample == s],
                    mt = if (length(mt_col)) md[[mt_col[1]]][md$.sample == s] else NA_real_)
    if (!nrow(d)) next
    g <- ggplot(d, aes(x, y, colour = mt)) + scattermore::geom_scattermore(pointsize = 2.4) +
      scale_x_log10() + scale_y_log10() +
      coord_cartesian(xlim = range(md$nCount_RNA), ylim = range(md$nFeature_RNA)) +
      labs(x = "UMIs", y = "Genes") + theme_minimal(base_size = 9) +
      theme(legend.position = "none", panel.grid.minor = element_blank())
    if (length(mt_col)) g <- g + scale_colour_viridis_c(option = "magma", direction = -1)
    qcp[[s]] <- uri(g, 4.4, 3.4)
  }
  render_stage("01", "Quality control",
    paste("Cell-level quality metrics for every sample, and how many cells each filter removed.",
          "Rendered from the run's own cell metadata — every cell is plotted."),
    list(violin = uri(gv, 9.6, 3.8), qc_panels = qcp,
         fate = csv(rp("qc", "cell_fate.csv")),
         summary = csv(rp("qc", "qc_summary_table.csv")),
         rho = csv(rp("qc", "soupx_rho.csv")),
         thresholds = list(min_features = QC$min_features, max_features = QC$max_features,
                           max_mt = if (!is.null(QC$max_mt)) QC$max_mt else NA)),
    rp("reports", "01-QC.html"))
}

# ============================ 02 — Doublets =================================
if ("02" %in% want) {
  fate <- csv(rp("qc", "cell_fate.csv"))
  gsc <- NULL
  if ("doublet_score" %in% colnames(md)) {
    gsc <- ggplot(data.frame(s = md$doublet_score, sample = md$.sample),
                  aes(s, fill = sample)) +
      geom_histogram(bins = 60, alpha = 0.85) + facet_wrap(~ sample, scales = "free_y") +
      labs(x = "scDblFinder score (retained singlets only)", y = "cells") +
      theme_minimal(base_size = 10) +
      theme(legend.position = "none", strip.text = element_text(face = "bold"))
  }
  render_stage("02", "Doublet detection",
    paste("Droplets that captured two cells are identified with scDblFinder and removed",
          "before integration. The counts below are what each sample lost."),
    list(fate = fate, score_hist = if (!is.null(gsc)) uri(gsc, 9.6, 4.6) else NULL,
         # stated rather than quietly omitted: the pre-removal per-cell calls are
         # not persisted by step 02, so the classic doublet UMAP cannot be rebuilt
         # from data here. See docs/report_product_design.md.
         limitation = paste("Per-cell doublet calls are not saved by step 02, so the",
                            "pre-removal doublet UMAP cannot be rendered from data in this",
                            "report. The score distribution shown covers retained singlets",
                            "only. The per-sample doublet UMAP remains in",
                            "02-Doublet_report.pdf.")),
    rp("reports", "02-Doublet.html"))
}

# ============================ 03 — Individual ===============================
if ("03" %in% want) {
  samps <- list.dirs(rp("individual"), full.names = FALSE, recursive = FALSE)
  samps <- samps[nzchar(samps)]
  per <- list()
  for (s in samps) {
    f <- rp("individual", s, paste0(s, "_seurat.rds"))
    mk <- csv(rp("individual", s, paste0(s, "_cluster_markers.csv")))
    pan <- NULL; nclust <- NA_integer_; ncell <- NA_integer_
    if (file.exists(f)) {
      so <- tryCatch(readRDS(f), error = function(e) NULL)
      if (!is.null(so) && "umap" %in% names(so@reductions)) {
        em <- Embeddings(so, "umap")
        cl <- as.character(so$seurat_clusters)
        d  <- data.frame(u1 = em[, 1], u2 = em[, 2], cl = cl)
        nclust <- length(unique(cl)); ncell <- nrow(d)
        g <- ggplot(d, aes(u1, u2, colour = cl)) +
          scattermore::geom_scattermore(pointsize = 2.6) +
          scale_colour_discrete(guide = "none") + theme_void(base_size = 9) +
          labs(title = sprintf("%s — %d clusters, %s cells", s, nclust,
                               format(ncell, big.mark = ","))) +
          theme(plot.title = element_text(size = 10, face = "bold"))
        pan <- uri(g, 4.4, 3.6)
      }
      rm(so); gc(verbose = FALSE)
    }
    top <- NULL
    if (!is.null(mk)) {
      gn <- intersect(c("gene", "feature"), names(mk)); cl <- intersect(c("cluster"), names(mk))
      sc <- intersect(c("avg_log2FC", "avg_logFC"), names(mk))
      if (length(gn) && length(cl) && length(sc)) {
        mk <- mk[order(mk[[cl]], -mk[[sc]]), ]
        top <- do.call(rbind, lapply(split(mk, mk[[cl]]), function(x)
          data.frame(Cluster = x[[cl]][1],
                     `Top markers` = paste(head(x[[gn]], 8), collapse = ", "),
                     check.names = FALSE)))
        rownames(top) <- NULL
      }
    }
    per[[s]] <- list(panel = pan, n_clusters = nclust, n_cells = ncell, top = top)
    msg("  03: %s (%s clusters)", s, nclust)
  }
  render_stage("03", "Per-sample analysis",
    paste("Each sample clustered on its own, before integration. Cluster counts differing",
          "widely between samples usually reflects differing cell numbers or complexity",
          "rather than biology."),
    list(per_sample = per), rp("reports", "03-Individual.html"))
}

# ============================ 04 — Annotation ===============================
if ("04" %in% want) {
  render_stage("04", "Cell-type annotation",
    paste("How each cluster was assigned an identity: the reference-based call, the",
          "marker-based call, and where they agreed."),
    list(clusters = csv(rp("annotation", "cluster_annotation_table.csv")),
         consensus = csv(rp("annotation", "consensus_annotation.csv")),
         vs = csv(rp("annotation", "singler_vs_sctype_comparison.csv")),
         label_src = if ("cell_type_ref" %in% colnames(md) &&
                         !all(is.na(md$cell_type_ref))) "frozen reference" else "de-novo"),
    rp("reports", "04-Annotation.html"))
}

# ============================ 05 — Integrated ===============================
if ("05" %in% want) {
  emb <- Embeddings(obj, "umap")
  d   <- data.frame(u1 = emb[, 1], u2 = emb[, 2], ct = md$.ct, sample = md$.sample)
  rng <- list(x = range(d$u1), y = range(d$u2))
  ras <- function(colour_by, pal) {
    g <- ggplot(d, aes(u1, u2, colour = .data[[colour_by]])) +
      scattermore::geom_scattermore(pointsize = 2) +
      scale_colour_manual(values = pal, guide = "none") +
      scale_x_continuous(limits = rng$x, expand = c(0, 0)) +
      scale_y_continuous(limits = rng$y, expand = c(0, 0)) +
      theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
    uri(g, 6.4, 5.4)
  }
  spal <- setNames(scales::hue_pal()(nlevels(d$sample)), levels(d$sample))
  cnt <- as.data.frame(table(sample = d$sample, cell_type = d$ct), stringsAsFactors = FALSE)
  names(cnt)[3] <- "n"; cnt$pct <- 100 * cnt$n / ave(cnt$n, cnt$sample, FUN = sum)
  cent <- do.call(rbind, lapply(split(d, d$ct), function(x)
    data.frame(ct = x$ct[1], u1 = median(x$u1), u2 = median(x$u2), n = nrow(x),
               stringsAsFactors = FALSE)))
  render_stage("05", "Integrated analysis",
    paste("All samples combined and corrected for sample-to-sample technical variation,",
          "then clustered and annotated together. This is the view the headline results",
          "come from."),
    list(umap_ct = ras("ct", PAL), umap_sample = ras("sample", spal),
         range = rng, centroids = cent, pal = PAL, spal = spal, counts = cnt,
         markers = csv(rp("integrated", "integrated_cluster_markers.csv"))),
    rp("reports", "05-Integrated.html"))
}

# ============================ index =========================================
idx <- rp("reports", "index.html")
items <- list(
  c("01-QC.html", "Quality control", "Per-cell metrics and what each filter removed"),
  c("02-Doublet.html", "Doublet detection", "Two-cell droplets identified and removed"),
  c("03-Individual.html", "Per-sample analysis", "Each sample clustered before integration"),
  c("04-Annotation.html", "Cell-type annotation", "How each cluster got its identity"),
  c("05-Integrated.html", "Integrated analysis", "All samples combined — the headline results"))
have <- Filter(function(x) file.exists(rp("reports", x[1])), items)
extra <- Filter(function(x) file.exists(rp("reports", x[1])), list(
  c("Customer_report.html", "Customer report", "The client-facing deliverable")))
# the internal report is named <run_name>_report.html, which is NOT
# basename(run_dir) — find it by pattern rather than guessing the stem
.int <- setdiff(list.files(rp("reports"), pattern = "_report\\.html$"), "Customer_report.html")
if (length(.int))
  extra <- c(extra, list(c(.int[1], "Interactive run report",
                           "Full internal report — every stage, with live plots")))
writeLines(c(
  "<!doctype html><meta charset='utf-8'><title>Reports</title>",
  "<style>body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,Helvetica,Arial,sans-serif;",
  "max-width:760px;margin:40px auto;padding:0 20px;color:#222;line-height:1.5}",
  "h1{font-size:1.5em;border-bottom:2px solid #eee;padding-bottom:.3em}",
  "a.card{display:block;text-decoration:none;color:inherit;border:1px solid #e4e4e4;",
  "border-radius:6px;padding:12px 16px;margin:10px 0}a.card:hover{border-color:#1a4f7a;",
  "background:#f7fafc}.t{font-weight:600;color:#1a4f7a}.d{font-size:.9em;color:#666}",
  ".meta{color:#666;font-size:.9em}</style>",
  sprintf("<h1>%s</h1><p class='meta'>Stage reports · generated %s</p>",
          htmltools::htmlEscape(basename(run_dir)), format(Sys.time(), "%Y-%m-%d %H:%M")),
  vapply(c(have, extra), function(x)
    sprintf("<a class='card' href='%s'><span class='t'>%s</span><br><span class='d'>%s</span></a>",
            x[1], htmltools::htmlEscape(x[2]), htmltools::htmlEscape(x[3])), character(1)),
  "<p class='meta'>The matching PDFs remain in this run's root for printing.</p>"), idx)
msg("index.html        %5.2f MB  (%d stage links)", file.info(idx)$size / 1e6, length(have))
cat("STAGE_REPORTS_DONE\n")
