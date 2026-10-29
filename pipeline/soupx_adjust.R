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
# Usage: Rscript pipeline/soupx_adjust.R [run_dir] [sample_base_dir] [--rho-max=0.3]
#   run_dir          defaults to $SCRNA_RESULTS_DIR (set by run_pipeline.sh, which
#                    exposes this script as step "00": `run_pipeline.sh cm <samples> 00`
#                    after a finished 01-05 run)
#   sample_base_dir  defaults to the parent of $SCRNA_SAMPLE1, else Samples/
#   --rho-max        samples whose estimated rho exceeds this are REPORTED but not
#                    corrected (accepted = FALSE in soupx_rho.csv). Default 0.3: above
#                    that the "cells" are mostly soup and subtracting it does not
#                    rescue them (H1D20_1, rho 0.795, is the motivating case).
#
# autoEstCont runs WITHOUT forceAccept. When SoupX's own sanity check rejects the
# estimate the sample is recorded with accepted = FALSE and rho = NA instead of
# silently taking a value SoupX itself did not trust.
# =============================================================================
suppressPackageStartupMessages({
  library(Matrix); library(SoupX); library(Seurat)
})

args    <- commandArgs(trailingOnly = TRUE)
flag    <- function(k, d) { v <- sub(paste0("^--", k, "="), "", grep(paste0("^--", k, "="), args, value = TRUE)); if (length(v)) v[1] else d }
pos     <- args[!grepl("^--", args)]
RUN     <- if (length(pos) >= 1) pos[1] else Sys.getenv("SCRNA_RESULTS_DIR", "")
if (!nzchar(RUN)) stop("usage: soupx_adjust.R <run_dir> [sample_base_dir] [--rho-max=0.3]")
base    <- if (length(pos) >= 2) pos[2] else {
  s1 <- Sys.getenv("SCRNA_SAMPLE1", ""); if (nzchar(s1)) dirname(s1) else "Samples" }
RHO_MAX <- as.numeric(flag("rho-max", "0.3"))
dest    <- "Samples_soupx"
dir.create(dest, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(RUN, "qc"), showWarnings = FALSE, recursive = TRUE)

rds <- file.path(RUN, "integrated", "integrated_annotated.rds")
if (!file.exists(rds)) stop("no integrated_annotated.rds in ", RUN,
                            " — SoupX needs this run's cluster assignments")
msg <- function(fmt, ...) message(sprintf(paste0("[soupx] ", fmt), ...))

s    <- readRDS(rds)
SAMP <- sort(unique(as.character(s$sample)))
# The run's own QC floors, so "cells passing QC after correction" uses the cut these
# cells already passed. Correction only lowers counts, so the floors are what can bite.
thr_f <- file.path(RUN, "qc", "qc_thresholds.csv")
THR   <- if (file.exists(thr_f)) read.csv(thr_f, stringsAsFactors = FALSE) else NULL
# Older runs predate qc_thresholds.csv: fall back to the species config's fixed floors
# (config.R pinned to this run dir, species from SCRNA_SPECIES, as 08b does).
if (is.null(THR)) {
  Sys.setenv(SCRNA_RESULTS_DIR = normalizePath(RUN))
  # Top level and global, as every step does: config.R finds its helpers through
  # sys.frame(1)$ofile and its species overlays read its globals. None of its names
  # collide with this script's.
  source(file.path(dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "config.R"))
  THR <- data.frame(Sample = SAMP, min_features = QC$min_features, min_counts = QC$min_counts)
  msg("no qc_thresholds.csv — using config QC floors (%d genes, %d UMIs)", THR$min_features[1], THR$min_counts[1])
}
n_qc  <- function(m, sm) {
  t <- THR[THR$Sample == sm, , drop = FALSE]
  if (!nrow(t)) return(NA_integer_)
  sum(colSums(m > 0) >= t$min_features[1] & colSums(m) >= t$min_counts[1])
}
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
  est <- tryCatch(autoEstCont(sc, doPlot = FALSE), error = function(e) e)
  before <- sum(toc[g, ])
  # Per-cell medians before/after correction: SoupX removes counts, never cells, so the
  # customer report shows what each cell keeps rather than a cell count that cannot move.
  row <- data.frame(sample = sm, rho = NA_real_, accepted = FALSE, cells = ncol(toc),
                    counts_before = before, counts_after = NA_real_,
                    pct_counts_removed = NA_real_,
                    med_umi_before = median(colSums(toc[g, ])), med_umi_after = NA_real_,
                    med_genes_before = median(colSums(toc[g, ] > 0)), med_genes_after = NA_real_,
                    cells_qc_before = n_qc(toc[g, ], sm), cells_qc_after = NA_integer_,
                    note = "", stringsAsFactors = FALSE)
  if (inherits(est, "error")) {
    row$note <- paste("autoEstCont rejected:", conditionMessage(est))
    msg("  NOT corrected — %s", row$note); rep[[sm]] <- row; next
  }
  rho <- est$fit$rhoEst; row$rho <- round(rho, 4)
  if (rho > RHO_MAX) {
    row$note <- sprintf("rho %.3f > rho_max %.2f: sample is mostly ambient, left uncorrected", rho, RHO_MAX)
    msg("  NOT corrected — %s", row$note); rep[[sm]] <- row; next
  }
  adj <- adjustCounts(est, roundToInt = TRUE)
  after <- sum(adj)
  msg("  rho=%.3f  counts %.1fM -> %.1fM (removed %.2f%%)  cells=%d",
      rho, before/1e6, after/1e6, 100*(before-after)/before, ncol(adj))
  outdir <- file.path(dest, paste0(sm, "_sx"), "filter_matrix")
  ftk <- F$ft[match(rownames(adj),
                    make.unique(if (ncol(F$ft) >= 2) F$ft[[2]] else F$ft[[1]])), , drop = FALSE]
  write_mtx(adj, ftk, colnames(adj), outdir)
  row$accepted <- TRUE; row$counts_after <- after
  row$pct_counts_removed <- round(100*(before-after)/before, 2)
  row$med_umi_after <- median(colSums(adj)); row$med_genes_after <- median(colSums(adj > 0))
  row$cells_qc_after <- n_qc(adj, sm)
  rep[[sm]] <- row
  msg("  wrote %s", outdir)
}

if (!length(rep)) stop("no sample had raw_matrix/ + filtered matrix under ", base)
res <- do.call(rbind, rep)
res$source <- sprintf("SoupX %s autoEstCont (no forceAccept, rho_max %.2f) on %s, %s",
                      as.character(utils::packageVersion("SoupX")), RHO_MAX,
                      basename(RUN), format(Sys.Date()))
out <- file.path(RUN, "qc", "soupx_rho.csv")
write.csv(res, out, row.names = FALSE)
msg("wrote %s", out)
cat("\n=== correction summary ===\n")
print(res[, c("sample", "rho", "accepted", "cells", "pct_counts_removed", "note")], row.names = FALSE)
ok <- res$sample[res$accepted]
if (length(ok)) {
  cat(sprintf("\nCorrected matrices: %s/<sample>_sx/filter_matrix  (%d of %d samples)\n",
              dest, length(ok), nrow(res)))
  cat("Re-run on corrected data with:\n  bash pipeline/run_pipeline.sh", Sys.getenv("SCRNA_SPECIES", "human"),
      paste(file.path(dest, paste0(ok, "_sx")), collapse = " "), "\n")
}
cat("SOUPX_ADJUST_DONE\n")
