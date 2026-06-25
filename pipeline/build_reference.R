# =============================================================================
# build_reference.R (v2) - Build a FROZEN, MARKER-VALIDATED cell-type reference
#   for run-independent label transfer (cross-run anchor benchmarking).
#
#   Trains a SingleR classifier ONCE from a curated, high-quality run so future
#   runs get reproducible labels via classifySingleR (see 05r_reference_transfer.R).
#   The trained model IS the frozen artifact.
#
#   Design: docs/frozen_reference_scope.md
#   v2 curation (vs v1 "AUTO-clusters-only", which dropped RBC/Platelet/CD8 T):
#     - keep ALL populations, selected by PER-CELL MARKER VALIDATION (a cell joins
#       the reference only if its label's defining bat markers are expressed) so
#       every reference cell is auditable and labels are not blindly inherited.
#     - recover CD8 T from CD4 T by CD8A/CD8B (cell_type collapses CD8 into CD4).
#     - hold out the anchor samples (so the anchor benchmark is not circular).
#     - save model + provenance + baseline anchor proportions + validation table.
#
#   Usage (set SCRNA_SPECIES=bat so marker overrides apply):
#     SCRNA_SPECIES=bat Rscript pipeline/build_reference.R <ref_run_dir> \
#        [--holdout=Aksh1,ES332] [--out=Results/frozen_reference] \
#        [--min-cells=20] [--validate-thr=0.25] [--cd8-thr=0.5] [--label-col=cell_type]
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
source(file.path(.pipeline_dir, "config.R"))   # MARKERS (bat overrides), CELLTYPE_COLORS

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
})

# ---- args -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: build_reference.R <ref_run_dir> [--holdout=A,B] [--out=DIR] [--min-cells=N] [--validate-thr=X] [--cd8-thr=X] [--label-col=COL]")
RUN_DIR <- args[1]
getflag <- function(k, d) { m <- grep(paste0("^--", k, "="), args, value = TRUE)
  if (length(m)) sub(paste0("^--", k, "="), "", m[1]) else d }
HOLDOUT      <- strsplit(getflag("holdout", "Aksh1,ES332"), ",")[[1]]
OUT_DIR      <- getflag("out", file.path("Results", "frozen_reference"))
MIN_CELLS    <- as.integer(getflag("min-cells", "20"))
VALIDATE_THR <- as.numeric(getflag("validate-thr", "0.25"))
CD8_THR      <- as.numeric(getflag("cd8-thr", "0.5"))
LABEL_COL    <- getflag("label-col", "cell_type")
SPECIES      <- Sys.getenv("SCRNA_SPECIES", "human")

stopifnot(dir.exists(RUN_DIR))
rds <- file.path(RUN_DIR, "integrated", "integrated_annotated.rds")
stopifnot(file.exists(rds))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("=== build_reference v2 (marker-validated) ===")
message("  run_dir : ", RUN_DIR)
message("  species : ", SPECIES, "   holdout: ", paste(HOLDOUT, collapse = ", "))

# ---- load -------------------------------------------------------------------
o <- readRDS(rds)
DefaultAssay(o) <- "RNA"
o <- NormalizeData(o, verbose = FALSE)
d <- GetAssayData(o, assay = "RNA", layer = "data")
md <- o@meta.data
n0 <- ncol(o)

if (!LABEL_COL %in% colnames(md)) {
  for (c in c("cell_type", "celltype", "consensus_label", "singler_label_clean"))
    if (c %in% colnames(md)) { LABEL_COL <- c; break }
}
stopifnot(LABEL_COL %in% colnames(md))
SAMPLE_COL <- if ("sample" %in% colnames(md)) "sample" else "orig.ident"
message("  label col: ", LABEL_COL, "   sample col: ", SAMPLE_COL, "   total cells: ", n0)

# ---- defining markers per label (from bat-overridden MARKERS) ----------------
present  <- function(gs) intersect(gs, rownames(d))
label_markers <- function(L) {
  g <- if (grepl("^CD8", L))               c(MARKERS$T_pan, MARKERS$CD8_T)
       else if (grepl("^CD4|Treg", L))     c(MARKERS$T_pan, MARKERS$CD4_T)
       else if (grepl("γδ|gamma", L)) c(MARKERS$T_pan, MARKERS$gamma_delta_T)
       else if (grepl("^NK", L))           MARKERS$NK
       else if (grepl("B cell|Plasma", L)) MARKERS$B_cell
       else if (grepl("FCGR3A", L))        MARKERS$FCGR3A_mono
       else if (grepl("Mono|Macro", L))    MARKERS$CD14_mono
       else if (grepl("Neutrophil", L))    MARKERS$Neutrophil
       else if (grepl("DC", L))            MARKERS$DC
       else if (grepl("Platelet", L))      MARKERS$Platelet
       else if (grepl("RBC", L))           MARKERS$RBC
       else if (grepl("HSPC|Progenitor", L)) MARKERS$HSPC
       else if (grepl("Baso", L))          MARKERS$Basophil
       else if (grepl("Eosino", L))        MARKERS$Eosinophil
       else if (grepl("Mast", L))          MARKERS$Mast_cell
       else character(0)
  present(g)
}
# mean expression of a gene set across given cells (0 if no genes)
mexpr <- function(genes, cells) {
  if (!length(genes)) return(rep(0, length(cells)))
  if (length(genes) == 1) as.numeric(d[genes, cells]) else Matrix::colMeans(d[genes, cells, drop = FALSE])
}

lab <- as.character(md[[LABEL_COL]]); lab[is.na(lab)] <- "Unassigned"

# ---- (A) recover CD8 T from CD4 T (cell_type collapses CD8 into CD4) ---------
cd8g <- present(c("CD8A", "CD8B")); cd4g <- present(c("CD4", "IL7R"))
if (length(cd8g)) {
  isCD4 <- grepl("^CD4 T", lab)
  cd8e <- mexpr(cd8g, seq_len(n0)); cd4e <- if (length(cd4g)) mexpr(cd4g, seq_len(n0)) else 0
  rec  <- isCD4 & cd8e > CD8_THR & cd8e > cd4e
  lab[rec] <- "CD8 T"
  message("  (A) CD8 T recovered from CD4 T (CD8A/B>", CD8_THR, "): ", sum(rec))
}

# ---- (B) per-cell MARKER VALIDATION: keep only marker-supported exemplars ----
keep <- rep(FALSE, n0)
val_score <- rep(NA_real_, n0)
for (L in unique(lab)) {
  if (L == "Unassigned") next
  cells <- which(lab == L)
  mk <- label_markers(L)
  if (!length(mk)) { message("  no markers for '", L, "' -> cannot validate, dropped"); next }
  s <- mexpr(mk, cells)
  val_score[cells] <- s
  keep[cells] <- s >= VALIDATE_THR
}
is_anchor <- as.character(md[[SAMPLE_COL]]) %in% HOLDOUT
keep <- keep & !is_anchor
message("  (B) marker-validated cells (>= ", VALIDATE_THR, ", anchors excluded): ", sum(keep))

# ---- (C) drop tiny labels ----------------------------------------------------
tab <- table(lab[keep]); small <- names(tab)[tab < MIN_CELLS]
if (length(small)) { keep <- keep & !(lab %in% small)
  message("  (C) drop labels < ", MIN_CELLS, " cells: ", paste(small, collapse = ", ")) }

ref_labels <- lab[keep]
message("\n  => reference cells: ", sum(keep), " / ", n0)
message("  reference label composition:")
print(sort(table(ref_labels), decreasing = TRUE))

# ---- marker-validation audit table (per label: n, mean/median defining marker)
val_tab <- do.call(rbind, lapply(sort(unique(ref_labels)), function(L) {
  cs <- which(keep & lab == L)
  data.frame(label = L, n = length(cs),
             mean_marker = round(mean(val_score[cs]), 2),
             median_marker = round(median(val_score[cs]), 2),
             markers = paste(label_markers(L), collapse = ","))
}))
message("\n  marker-validation audit (defining-marker expression per label):")
print(val_tab[, c("label","n","mean_marker","median_marker")], row.names = FALSE)

# ---- train ------------------------------------------------------------------
ref_mat <- d[, keep, drop = FALSE]
message("\n[", format(Sys.time(), "%H:%M:%S"), "] trainSingleR (",
        nrow(ref_mat), " genes x ", ncol(ref_mat), " cells)...")
trained <- trainSingleR(ref = ref_mat, labels = ref_labels, de.method = "classic")
message("[", format(Sys.time(), "%H:%M:%S"), "] trained.")

# ---- baseline: classify the held-out anchors WITH the frozen model ----------
baseline <- list()
if (any(is_anchor)) {
  message("\nclassifying held-out anchors with the frozen model (baseline)...")
  for (s in HOLDOUT) {
    sel <- is_anchor & as.character(md[[SAMPLE_COL]]) == s
    if (!any(sel)) { message("  ", s, ": not present"); next }
    p <- classifySingleR(test = d[, sel, drop = FALSE], trained = trained, fine.tune = TRUE)
    prop <- sort(round(100 * table(p$pruned.labels, useNA = "no") / sum(!is.na(p$pruned.labels)), 2),
                 decreasing = TRUE)
    baseline[[s]] <- list(n = sum(sel),
                          pct_unassigned = round(100 * mean(is.na(p$pruned.labels)), 2),
                          proportions = prop)
    message("  ", s, " (n=", sum(sel), ", unassigned ", baseline[[s]]$pct_unassigned, "%):")
    print(prop)
  }
}

# ---- provenance + save ------------------------------------------------------
git_commit <- tryCatch(trimws(system("git rev-parse --short HEAD", intern = TRUE)),
                       error = function(e) NA_character_)
stamp <- format(Sys.Date())
meta <- list(
  built = as.character(Sys.time()), date = stamp, git = git_commit,
  method = "v2-marker-validated", source_run = RUN_DIR, species = SPECIES,
  label_col = LABEL_COL, sample_col = SAMPLE_COL, holdout = HOLDOUT,
  filters = list(min_cells = MIN_CELLS, validate_thr = VALIDATE_THR, cd8_thr = CD8_THR),
  n_total = n0, n_reference = sum(keep),
  labels = as.list(sort(table(ref_labels), decreasing = TRUE)),
  validation = val_tab, genes = rownames(ref_mat)
)
bundle <- list(model = trained, meta = meta, baseline = baseline)

base <- sprintf("frozen_ref_%s_%s_%s_v2",
                gsub("[^A-Za-z0-9]", "", basename(RUN_DIR)), SPECIES, stamp)
model_path <- file.path(OUT_DIR, paste0(base, ".rds"))
saveRDS(bundle, model_path)
writeLines(capture.output(str(meta, max.level = 2, list.len = 300)),
           file.path(OUT_DIR, paste0(base, ".meta.txt")))
write.csv(val_tab, file.path(OUT_DIR, paste0(base, ".validation.csv")), row.names = FALSE)
if (length(baseline)) {
  bdf <- do.call(rbind, lapply(names(baseline), function(s) {
    pr <- baseline[[s]]$proportions
    data.frame(sample = s, cell_type = names(pr), pct = as.numeric(pr))
  }))
  write.csv(bdf, file.path(OUT_DIR, paste0(base, ".baseline.csv")), row.names = FALSE)
}

message("\n=== DONE (v2) ===")
message("  model     : ", model_path)
message("  validation: ", file.path(OUT_DIR, paste0(base, ".validation.csv")))
message("  baseline  : ", file.path(OUT_DIR, paste0(base, ".baseline.csv")))
