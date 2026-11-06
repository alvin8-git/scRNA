# =============================================================================
# 08f_customer_export.R — the analysis data behind the report, for the customer
#
# The HTML report shows conclusions; this exports the data they rest on so a
# customer (or a bioinformatician they hire later) can verify any claim without
# re-running anything.
#
# Writes into <run>/export/:
#   per_cell_assignments.csv.gz   one row per cell — the verification backbone
#   pseudobulk_cpm.csv.gz         gene x sample expression
#   <run>_analysis.xlsx           every derived table, one sheet each
#
# Everything here is derived from the run's own object and CSVs, so it is
# deterministic and works for any validated sample type (human PBMC, human whole
# blood, cardiac, bat blood, bat wing). Nothing is hardcoded to a species: the
# marker panel comes from config, labels from the object.
#
# Also writes the analysis object as a SUPPLEMENT to the HTML report (too large to
# embed or email): <run>_seurat_slim.rds (counts + UMAP + labels, xz) and
# <run>.h5ad (AnnData, written by hand with hdf5r; X = raw counts, CSR).
# xz is ~8x slower than gzip, so on large runs this costs minutes: skip with
# SCRNA_NO_OBJECT_EXPORT=1 or --no-object.
#
# Usage: Rscript pipeline/08f_customer_export.R <run_dir> [--no-object]
# Env:   SCRNA_EXPORT_DIR  write export/ here instead of <run_dir>/export
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
pos     <- args[!startsWith(args, "--")]
run_dir <- if (length(pos)) pos[1] else stop("usage: 08f_customer_export.R <run_dir> [--no-object]")
no_obj  <- "--no-object" %in% args || nzchar(Sys.getenv("SCRNA_NO_OBJECT_EXPORT"))
Sys.setenv(SCRNA_RESULTS_DIR = run_dir)
source(file.path(.pipeline_dir, "config.R"))
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

msg <- function(fmt, ...) message(sprintf(paste0("[08f] ", fmt), ...))
rp  <- function(...) file.path(run_dir, ...)
exp_dir <- Sys.getenv("SCRNA_EXPORT_DIR", rp("export")); dir.create(exp_dir, showWarnings = FALSE, recursive = TRUE)
csv <- function(p) if (file.exists(p)) tryCatch(read.csv(p, stringsAsFactors = FALSE),
                                                error = function(e) NULL) else NULL

rds <- rp("integrated", "integrated_annotated.rds")
if (!file.exists(rds)) stop("no integrated_annotated.rds in ", run_dir)
msg("reading %s", rds)
obj <- readRDS(rds)
md  <- obj@meta.data
lab <- if ("cell_type_ref" %in% colnames(md) && !all(is.na(md$cell_type_ref)))
         "cell_type_ref" else "cell_type"
sl  <- {
  # brace-wrapped: a bare multi-line if/else at top level parses as two
  # statements once the true branch completes, giving "unexpected 'else'"
  if (exists("SAMPLE_NAMES")) {
    c(intersect(SAMPLE_NAMES, unique(as.character(md$sample))),
      sort(setdiff(unique(as.character(md$sample)), SAMPLE_NAMES)))
  } else {
    sort(unique(as.character(md$sample)))
  }
}
msg("%s cells | %d samples | labels from %s", format(nrow(md), big.mark = ","),
    length(sl), lab)

# ---- 1. per-cell assignments: the sheet that makes everything checkable -----
emb <- tryCatch(Embeddings(obj, "umap"), error = function(e) NULL)
cells <- data.frame(
  cell_barcode = rownames(md),
  sample       = as.character(md$sample),
  cluster      = as.character(md$seurat_clusters),
  cell_type    = as.character(md[[lab]]),
  n_umis       = md$nCount_RNA,
  n_genes      = md$nFeature_RNA,
  stringsAsFactors = FALSE)
for (c0 in c("percent.mt", "percent.ribo", "doublet_score"))
  if (c0 %in% colnames(md)) cells[[sub("\\.", "_", c0)]] <- round(md[[c0]], 4)
if ("cell_type_denovo" %in% colnames(md))
  cells$cell_type_denovo <- as.character(md$cell_type_denovo)
if (!is.null(emb)) { cells$UMAP_1 <- round(emb[, 1], 4); cells$UMAP_2 <- round(emb[, 2], 4) }
f1 <- file.path(exp_dir, "per_cell_assignments.csv.gz")
gz <- gzfile(f1, "w"); write.csv(cells, gz, row.names = FALSE); close(gz)
msg("per_cell_assignments.csv.gz  %.2f MB  (%s rows x %d cols)",
    file.info(f1)$size/1e6, format(nrow(cells), big.mark = ","), ncol(cells))

# ---- 2. pseudobulk CPM, gene x sample ---------------------------------------
cmat <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
pb   <- sapply(split(seq_len(ncol(cmat)), factor(cells$sample, levels = sl)),
               function(i) Matrix::rowSums(cmat[, i, drop = FALSE]))
cpm  <- round(t(t(pb) / colSums(pb)) * 1e6, 3)
pbdf <- data.frame(gene = rownames(cpm), cpm, check.names = FALSE)
pbdf <- pbdf[order(-rowMeans(cpm)), ]
f2 <- file.path(exp_dir, "pseudobulk_cpm.csv.gz")
gz <- gzfile(f2, "w"); write.csv(pbdf, gz, row.names = FALSE); close(gz)
msg("pseudobulk_cpm.csv.gz        %.2f MB  (%s genes)",
    file.info(f2)$size/1e6, format(nrow(pbdf), big.mark = ","))

# ---- 3. derived tables for the workbook ------------------------------------
mt <- intersect(c("percent.mt", "percent_mt"), colnames(md))
overview <- do.call(rbind, lapply(sl, function(s) {
  i <- cells$sample == s
  data.frame(Sample = s, Cells = sum(i),
             `Pct of run` = round(100*sum(i)/nrow(cells), 2),
             `Median genes` = round(median(cells$n_genes[i])),
             `Median UMIs` = round(median(cells$n_umis[i])),
             `Median pct mito` = if (length(mt)) round(median(md[[mt[1]]][i]), 2) else NA_real_,
             `Cell types` = length(unique(cells$cell_type[i])),
             check.names = FALSE)
}))

ct_counts <- as.data.frame(table(Sample = factor(cells$sample, levels = sl),
                                 `Cell type` = cells$cell_type), stringsAsFactors = FALSE)
names(ct_counts)[3] <- "Cells"
ct_counts$`Pct of sample` <- round(100*ct_counts$Cells /
                                   ave(ct_counts$Cells, ct_counts$Sample, FUN = sum), 2)

# marker panel from config — records which panel genes exist in the data, so the
# panel itself is auditable rather than implicit
# MARKERS is the NAMED list of panels; ALL_MARKERS is a flat character vector of
# every gene. Using the latter here gave names() == NULL and silently emptied
# both this table and the cluster-evidence one.
panel_list <- {
  if (exists("MARKERS") && is.list(MARKERS) && length(MARKERS)) {
    # MARKERS also carries config flags (e.g. compute_integrated = TRUE); keep
    # only character vectors or a logical leaks in as a one-gene "panel"
    MARKERS[vapply(MARKERS, is.character, logical(1))]
  } else if (exists("ALL_MARKERS") && length(ALL_MARKERS)) {
    list(`marker panel` = as.character(ALL_MARKERS))
  } else {
    list()
  }
}
panel <- if (length(panel_list))
  do.call(rbind, lapply(names(panel_list), function(p)
    data.frame(Panel = p, Gene = panel_list[[p]],
               `Present in data` = panel_list[[p]] %in% rownames(cmat),
               check.names = FALSE, stringsAsFactors = FALSE))) else NULL

# per-cluster evidence: label + the markers that decided it
evidence <- NULL
if (!is.null(panel)) {
  pg  <- intersect(unique(panel$Gene), rownames(cmat))
  sub <- cmat[pg, , drop = FALSE]
  ord <- order(suppressWarnings(as.integer(unique(cells$cluster))), unique(cells$cluster))
  # Detection % across the whole run, so a cluster's markers can be reported as
  # ENRICHMENT rather than raw prevalence. Ranking by raw detection surfaces the
  # most ubiquitous genes (KRT8, COL1A1) which discriminate nothing — the useful
  # markers are those detected far more in this cluster than elsewhere.
  pc_all <- round(100 * Matrix::rowSums(sub > 0) / ncol(sub), 1)
  evidence <- do.call(rbind, lapply(unique(cells$cluster)[ord], function(k) {
    i    <- which(cells$cluster == k)
    pc   <- round(100*Matrix::rowSums(sub[, i, drop = FALSE] > 0)/length(i), 0)
    lbl  <- names(sort(table(cells$cell_type[i]), decreasing = TRUE))[1]
    # markers of the panel matching the assigned label — the direct evidence
    own  <- if (!is.null(panel_list[[lbl]])) intersect(panel_list[[lbl]], names(pc)) else character(0)
    own_s <- if (length(own)) {
      o <- sort(pc[own], decreasing = TRUE)
      paste(sprintf("%s %d%%", names(o), as.integer(o)), collapse = ", ")
    } else NA_character_
    # most enriched genes relative to the rest of the run, regardless of panel
    enr  <- pc - pc_all[names(pc)]
    top  <- head(sort(enr, decreasing = TRUE), 6)
    data.frame(
      Cluster = k, Cells = length(i), `Cell type` = lbl,
      `Markers of the assigned type (pct of cluster cells)` = own_s,
      `Most enriched vs rest of run (pct points above run average)` =
        paste(sprintf("%s +%d", names(top), as.integer(round(top))), collapse = ", "),
      check.names = FALSE, stringsAsFactors = FALSE)
  }))
}

mk <- csv(rp("integrated", "integrated_cluster_markers.csv"))
if (!is.null(mk)) {
  sc <- intersect(c("avg_log2FC", "avg_logFC"), names(mk)); cl <- intersect("cluster", names(mk))
  if (length(sc) && length(cl)) {
    mk <- mk[order(mk[[cl]], -mk[[sc[1]]]), ]
    mk <- do.call(rbind, lapply(split(mk, mk[[cl]]), head, 25)); rownames(mk) <- NULL
  }
}

# Provenance MUST describe the run, not this export process. Re-sourcing config.R
# here gives the export's own state: with no SCRNA_SAMPLE* set, SAMPLE_NAMES
# falls back to H1/H2, no cluster-map guard matches, and the export would report
# "curated map applied: FALSE" for a run that demonstrably used one. The run's
# own step-05 log is the authoritative record of what actually happened.
log05 <- rp("logs", "05_05_annotate.log")
lg    <- if (file.exists(log05)) readLines(log05, warn = FALSE) else character(0)
from_log <- function(pat, default = NA_character_) {
  h <- grep(pat, lg, value = TRUE)
  if (!length(h)) return(default)
  trimws(sub("^[^:]*:\\s*", "", sub("\\s*$", "", h[1])))
}
params <- data.frame(Parameter = character(), Value = character(), stringsAsFactors = FALSE)
add <- function(k, v) params[nrow(params)+1, ] <<- c(k, as.character(v))
add("Run", basename(run_dir))
add("Species profile", {
  s <- grep("^\\[Species\\]", lg, value = TRUE)
  if (length(s)) trimws(sub("^\\[Species\\]\\s*", "", s[1])) else
    Sys.getenv("SCRNA_SPECIES", "human (assumed — not recorded in log)")
})
add("Cells analysed", format(nrow(cells), big.mark = ",")); add("Samples", length(sl))
add("Samples analysed", paste(sl, collapse = ", "))
add("Label source", lab); add("Clusters", length(unique(cells$cluster)))
add("QC min genes per cell", QC$min_features); add("QC max genes per cell", QC$max_features)
if (!is.null(QC$max_mt)) add("QC max pct mitochondrial", QC$max_mt)
add("SingleR reference", from_log("Reference:", get0("SINGLER_REF", ifnotfound = NA)))
add("Curated cluster map applied",
    if (!length(lg)) "unknown (no step-05 log found)"
    else any(grepl("Using manual CLUSTER_CELLTYPE_MAP", lg)))
add("Marker panels used", length(panel_list))
add("Exported", format(Sys.time(), "%Y-%m-%d %H:%M"))
for (p in c("R", "Seurat", "harmony", "SingleR", "scDblFinder", "SoupX"))
  add(paste0("Version: ", p), if (p == "R") paste(R.version$major, R.version$minor, sep = ".")
      else tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed"))

sheets <- list(
  Run_parameters   = params,
  Sample_overview  = overview,
  CellType_counts  = ct_counts,
  Cluster_evidence = evidence,
  Cluster_markers  = mk,
  Marker_panel     = panel,
  QC_cell_fate     = csv(rp("qc", "cell_fate.csv")),
  QC_summary       = csv(rp("qc", "qc_summary_table.csv")),
  Ambient_rho      = csv(rp("qc", "soupx_rho.csv")),
  Annotation_detail = csv(rp("annotation", "cluster_annotation_table.csv")))
sheets <- sheets[!vapply(sheets, function(x) is.null(x) || !NROW(x), logical(1))]
for (nm in names(sheets)) {
  f <- file.path(exp_dir, paste0(tolower(nm), ".csv"))
  write.csv(sheets[[nm]], f, row.names = FALSE)
}
msg("wrote %d derived tables as CSV", length(sheets))
# ---- 4. analysis object supplement (never fails the export) -----------------
if (no_obj) {
  msg("object export skipped (--no-object / SCRNA_NO_OBJECT_EXPORT)")
} else {
  tryCatch({
    source(file.path(.pipeline_dir, "object_export.R"))
    export_objects(obj, md, lab, exp_dir, basename(run_dir), msg)
  }, error = function(e) msg("object export skipped: %s", conditionMessage(e)))
}

saveRDS(list(sheets = sheets, run = basename(run_dir), exp_dir = exp_dir),
        file.path(tempdir(), "08f_sheets.rds"))
cat(sprintf("EXPORT_TABLES_DONE %d %s\n", length(sheets), exp_dir))
