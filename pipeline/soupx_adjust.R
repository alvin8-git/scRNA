# =============================================================================
# soupx_adjust.R — estimate and remove ambient RNA for every sample in a run
#
# Reads a FINISHED run (needs integrated_annotated.rds for its cluster
# assignments, which SoupX uses to build the soup profile), estimates the
# contamination fraction rho per sample with autEstCont, and writes corrected
# matrices to Samples_soupx/<sample>_sx/filter_matrix for re-running the
# pipeline on.
#
# Two outputs, and the first one matters for the reports:
#   <run_dir>/qc/soupx_rho.csv   sample + rho + diagnostics — this is what the
#                                customer report's Ambient RNA section reads
#   Samples_soupx/<s>_sx/        corrected counts, as their own sample dirs
#
# This replaces an ad-hoc script that wrote one HARDCODED filename
# (/data/alvin/tmp/soupx_adjust_pbmc_report.csv) for every cohort. Running it
# for a second cohort silently overwrote the first — which destroyed the PBMC
# rho values on 2026-10-24. Writing into the run directory makes that
# impossible: each run owns its own file.
#
# Originals are never modified; the _sx names give corrected data its own run
# directory and its own sample_cache entries.
#
# Usage: Rscript pipeline/soupx_adjust.R <run_dir> [sample_base_dir]
#   sample_base_dir defaults to Samples/ — pass e.g. Samples/Cardio when the
#   run's inputs live in a subdirectory.
# =============================================================================
suppressPackageStartupMessages({
  library(Matrix); library(SoupX); library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: soupx_adjust.R <run_dir> [sample_base_dir]")
RUN  <- args[1]
base <- if (length(args) >= 2) args[2] else "Samples"
dest <- "Samples_soupx"
dir.create(dest, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(RUN, "qc"), showWarnings = FALSE, recursive = TRUE)

rds <- file.path(RUN, "integrated", "integrated_annotated.rds")
if (!file.exists(rds)) stop("no integrated_annotated.rds in ", RUN,
                            " — SoupX needs this run's cluster assignments")
msg <- function(fmt, ...) message(sprintf(paste0("[soupx] ", fmt), ...))

s    <- readRDS(rds)
SAMP <- sort(unique(as.character(s$sample)))
msg("run %s | %d samples | base %s", basename(RUN), length(SAMP), base)

filtered_dir <- function(sm) {
  d <- file.path(base, sm)
  for (sub in c("filter_matrix", "filtered_feature_bc_matrix"))
    if (dir.exists(file.path(d, sub))) return(file.path(d, sub))
  if (file.exists(file.path(d, "matrix.mtx.gz"))) return(d)
  NA_character_
}
read_mtx <- function(d) {
  m  <- readMM(gzfile(file.path(d, "matrix.mtx.gz")))
  ft <- read.delim(gzfile(file.path(d, "features.tsv.gz")), header = FALSE,
                   stringsAsFactors = FALSE)
  bc <- readLines(gzfile(file.path(d, "barcodes.tsv.gz")))
  rownames(m) <- make.unique(if (ncol(ft) >= 2) ft[[2]] else ft[[1]])
  colnames(m) <- bc
  list(m = as(m, "CsparseMatrix"), ft = ft)
}
write_mtx <- function(m, ft, bc, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  writeMM(as(m, "CsparseMatrix"), file.path(outdir, "matrix.mtx"))
  system2("gzip", c("-f", shQuote(file.path(outdir, "matrix.mtx"))))
  write.table(ft, gzfile(file.path(outdir, "features.tsv.gz")), sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  writeLines(bc, gzfile(file.path(outdir, "barcodes.tsv.gz")))
}

rep <- list()
for (sm in SAMP) {
  msg("=== %s ===", sm)
  rawd <- file.path(base, sm, "raw_matrix"); fd <- filtered_dir(sm)
  if (!dir.exists(rawd) || is.na(fd)) {
    msg("  SKIP — needs both raw_matrix/ and a filtered matrix under %s",
        file.path(base, sm))
    next
  }
  R <- read_mtx(rawd); F <- read_mtx(fd)
  incell <- colnames(s)[s$sample == sm]
  bc     <- sub(paste0("^", sm, "_"), "", incell)
  keep   <- bc %in% colnames(F$m)
  toc    <- F$m[, bc[keep], drop = FALSE]
  clus   <- as.character(s$seurat_clusters[incell[keep]]); names(clus) <- colnames(toc)
  g      <- intersect(rownames(R$m), rownames(toc))
  sc <- SoupChannel(R$m[g, , drop = FALSE], toc[g, , drop = FALSE])
  sc <- setClusters(sc, clus)
  sc <- autoEstCont(sc, doPlot = FALSE, forceAccept = TRUE)
  rho <- sc$fit$rhoEst
  adj <- adjustCounts(sc, roundToInt = TRUE)
  before <- sum(toc[g, ]); after <- sum(adj)
  msg("  rho=%.3f  counts %.1fM -> %.1fM (removed %.2f%%)  cells=%d",
      rho, before/1e6, after/1e6, 100*(before-after)/before, ncol(adj))
  outdir <- file.path(dest, paste0(sm, "_sx"), "filter_matrix")
  ftk <- F$ft[match(rownames(adj),
                    make.unique(if (ncol(F$ft) >= 2) F$ft[[2]] else F$ft[[1]])), , drop = FALSE]
  write_mtx(adj, ftk, colnames(adj), outdir)
  rep[[sm]] <- data.frame(sample = sm, rho = round(rho, 4), cells = ncol(adj),
                          counts_before = before, counts_after = after,
                          pct_counts_removed = round(100*(before-after)/before, 2),
                          stringsAsFactors = FALSE)
  msg("  wrote %s", outdir)
}

if (!length(rep)) stop("no sample could be corrected — check raw_matrix/ paths under ", base)
res <- do.call(rbind, rep)
res$source <- sprintf("SoupX %s autoEstCont on %s, %s",
                      as.character(utils::packageVersion("SoupX")),
                      basename(RUN), format(Sys.Date()))
out <- file.path(RUN, "qc", "soupx_rho.csv")
write.csv(res, out, row.names = FALSE)
msg("wrote %s", out)
cat("\n=== correction summary ===\n")
print(res[, c("sample", "rho", "cells", "pct_counts_removed")], row.names = FALSE)
cat(sprintf("\nCorrected matrices: %s/<sample>_sx/filter_matrix\n", dest))
cat("SOUPX_ADJUST_DONE\n")
