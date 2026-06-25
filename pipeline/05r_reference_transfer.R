# =============================================================================
# 05r_reference_transfer.R - Run-INDEPENDENT label transfer onto a finished run.
#
#   Classifies a run's cells against the FROZEN bat reference (built by
#   build_reference.R) so labels do not depend on the run's sample mix. ADDITIVE:
#   writes cell_type_ref alongside the de-novo cell_type; never replaces it.
#
#   Gated on REFERENCE_MODEL (config.R / SCRNA_REFERENCE_MODEL). No model -> skip.
#   Design: docs/frozen_reference_scope.md
#
#   Outputs (in <run_dir>/annotation/):
#     reference_transfer_cells.csv.gz   barcode,sample,cluster,cell_type_ref,ref_pruned,ref_delta
#     reference_transfer_composition.csv  per-sample cell_type_ref %  (+ ALL row)
#
#   Usage:
#     SCRNA_SPECIES=bat SCRNA_REFERENCE_MODEL=<model.rds> \
#       Rscript pipeline/05r_reference_transfer.R <run_dir> [--model=PATH] [--no-fine-tune]
# =============================================================================
.pipeline_dir <- local({
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(f)) dirname(f)
  else { a <- commandArgs(trailingOnly = FALSE)
    d <- sub("--file=", "", a[grep("--file=", a)]); if (length(d)) dirname(normalizePath(d)) else "." }
})
source(file.path(.pipeline_dir, "config.R"))
suppressPackageStartupMessages({ library(Seurat); library(SingleR); library(BiocParallel) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: 05r_reference_transfer.R <run_dir> [--model=PATH] [--no-fine-tune]")
RUN_DIR <- args[1]
getflag <- function(k, d) { m <- grep(paste0("^--", k, "="), args, value = TRUE)
  if (length(m)) sub(paste0("^--", k, "="), "", m[1]) else d }
MODEL <- getflag("model", if (exists("REFERENCE_MODEL")) REFERENCE_MODEL else "")
FT    <- !("--no-fine-tune" %in% args) && (if (exists("REF_FINE_TUNE")) isTRUE(REF_FINE_TUNE) else TRUE)

if (!nzchar(MODEL)) { message("05r: no REFERENCE_MODEL set — skipping label transfer."); quit(status = 0) }
if (!file.exists(MODEL)) { message("05r: REFERENCE_MODEL not found: ", MODEL, " — skipping."); quit(status = 0) }
rds <- file.path(RUN_DIR, "integrated", "integrated_annotated.rds")
if (!file.exists(rds)) { message("05r: no integrated_annotated.rds — skipping."); quit(status = 0) }

message("=== 05r reference transfer ===\n  run  : ", RUN_DIR, "\n  model: ", MODEL,
        "\n  fine.tune: ", FT)
b <- readRDS(MODEL)
o <- readRDS(rds); DefaultAssay(o) <- "RNA"; o <- NormalizeData(o, verbose = FALSE)
d <- GetAssayData(o, assay = "RNA", layer = "data")

# align test genes to the trained reference gene set (pre-trained model needs exact match)
g <- b$meta$genes
miss <- setdiff(g, rownames(d))
if (length(miss)) {
  z <- Matrix::Matrix(0, nrow = length(miss), ncol = ncol(d), sparse = TRUE,
                      dimnames = list(miss, colnames(d)))
  d <- rbind(d[intersect(g, rownames(d)), , drop = FALSE], z)
}
d <- d[g, , drop = FALSE]
message("  genes: ", length(g), " (", length(miss), " zero-filled)")

BP <- tryCatch(MulticoreParam(workers = min(8, max(1, parallel::detectCores() - 2))),
               error = function(e) SerialParam())
message("[", format(Sys.time(), "%H:%M:%S"), "] classifySingleR (fine.tune=", FT, ") ...")
p <- classifySingleR(test = d, trained = b$model, fine.tune = FT, prune = TRUE, BPPARAM = BP)
message("[", format(Sys.time(), "%H:%M:%S"), "] done.")

# per-cell confidence: top1-top2 score gap (delta)
gap <- apply(p$scores, 1, function(x){ s <- sort(x, decreasing = TRUE); if (length(s) >= 2) s[1]-s[2] else NA })
sample_col <- if ("sample" %in% colnames(o@meta.data)) "sample" else "orig.ident"
cells <- data.frame(
  barcode      = colnames(o),
  sample       = as.character(o@meta.data[[sample_col]]),
  cluster      = as.character(o$seurat_clusters),
  cell_type_ref= ifelse(is.na(p$pruned.labels), "Unassigned", p$pruned.labels),
  ref_label    = p$labels,
  ref_delta    = round(gap, 4),
  stringsAsFactors = FALSE)

dir.create(file.path(RUN_DIR, "annotation"), showWarnings = FALSE, recursive = TRUE)
cells_path <- file.path(RUN_DIR, "annotation", "reference_transfer_cells.csv.gz")
write.csv(cells, gzfile(cells_path), row.names = FALSE)

# per-sample composition (%) of cell_type_ref, plus an ALL row
comp_one <- function(v) round(100 * prop.table(table(factor(v))), 2)
samps <- sort(unique(cells$sample))
types <- sort(unique(cells$cell_type_ref))
comp <- do.call(rbind, lapply(c("ALL", samps), function(s) {
  v <- if (s == "ALL") cells$cell_type_ref else cells$cell_type_ref[cells$sample == s]
  row <- setNames(rep(0, length(types)), types); t <- comp_one(v); row[names(t)] <- as.numeric(t)
  data.frame(sample = s, n = length(v), t(row), check.names = FALSE)
}))
comp_path <- file.path(RUN_DIR, "annotation", "reference_transfer_composition.csv")
write.csv(comp, comp_path, row.names = FALSE)

message("\n  overall cell_type_ref composition (%):")
print(sort(comp_one(cells$cell_type_ref), decreasing = TRUE))
message("  pct Unassigned: ", round(100 * mean(cells$cell_type_ref == "Unassigned"), 2), "%")
message("\n=== DONE ===\n  cells      : ", cells_path, "\n  composition: ", comp_path)
