# =============================================================================
# 08c_benchmark_concordance.R - Cross-run anchor benchmark + sort-status readout.
#
#   Compares the anchor samples' frozen-reference cell-type composition in THIS
#   run against the reference's stored baseline (anchors classified by the same
#   model when it was built). Identical input + same model => proportions should
#   match within DRIFT_FLAG_PP; large drift = run/pipeline artifact to investigate.
#
#   Also prints a per-batch RBC/Platelet/Neutrophil readout (whole-blood signature:
#   high => presort/whole blood; depleted => postsort PBMC).
#
#   Requires 05r_reference_transfer.R to have run first. Gated on REFERENCE_MODEL.
#   Design: docs/frozen_reference_scope.md
#
#   Outputs (in <run_dir>/benchmark/):
#     concordance.csv   anchor x cell_type: run% vs baseline% vs delta vs flag
#     benchmark_report.md
#
#   Usage:
#     SCRNA_REFERENCE_MODEL=<model.rds> Rscript pipeline/08c_benchmark_concordance.R \
#        <run_dir> [--model=PATH] [--anchors=A,B] [--drift=5]
# =============================================================================
.pipeline_dir <- local({
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(f)) dirname(f)
  else { a <- commandArgs(trailingOnly = FALSE)
    d <- sub("--file=", "", a[grep("--file=", a)]); if (length(d)) dirname(normalizePath(d)) else "." }
})
source(file.path(.pipeline_dir, "config.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: 08c_benchmark_concordance.R <run_dir> [--model=PATH] [--anchors=A,B] [--drift=5]")
RUN_DIR <- args[1]
getflag <- function(k, d) { m <- grep(paste0("^--", k, "="), args, value = TRUE)
  if (length(m)) sub(paste0("^--", k, "="), "", m[1]) else d }
MODEL   <- getflag("model", if (exists("REFERENCE_MODEL")) REFERENCE_MODEL else "")
ANCHORS <- { a <- getflag("anchors", paste(if (exists("ANCHOR_SAMPLES")) ANCHOR_SAMPLES else character(0), collapse = ","))
             if (nzchar(a)) strsplit(a, ",")[[1]] else character(0) }
DRIFT   <- as.numeric(getflag("drift", if (exists("DRIFT_FLAG_PP")) DRIFT_FLAG_PP else 5))

cells_path <- file.path(RUN_DIR, "annotation", "reference_transfer_cells.csv.gz")
if (!file.exists(cells_path)) { message("08c: no reference_transfer_cells.csv.gz (run 05r first) — skipping."); quit(status = 0) }
if (!nzchar(MODEL) || !file.exists(MODEL)) { message("08c: REFERENCE_MODEL missing — skipping."); quit(status = 0) }
if (!length(ANCHORS)) { message("08c: no ANCHOR_SAMPLES set — skipping."); quit(status = 0) }

cells <- read.csv(cells_path, stringsAsFactors = FALSE)
base  <- readRDS(MODEL)$baseline
dir.create(file.path(RUN_DIR, "benchmark"), showWarnings = FALSE, recursive = TRUE)

# proportions (%) over non-Unassigned cells (matches baseline computation)
prop_excl <- function(v) { v <- v[v != "Unassigned"]; round(100 * prop.table(table(v)), 2) }

message("=== 08c benchmark concordance ===\n  anchors: ", paste(ANCHORS, collapse = ", "),
        "   drift flag: ", DRIFT, " pp\n")

rows <- list(); max_drift <- 0; any_flag <- FALSE
for (s in ANCHORS) {
  if (!s %in% cells$sample) { message("  ", s, ": not in this run — skip"); next }
  if (is.null(base[[s]])) { message("  ", s, ": no baseline in model — skip"); next }
  run_p  <- prop_excl(cells$cell_type_ref[cells$sample == s])
  base_p <- base[[s]]$proportions
  types  <- sort(union(names(run_p), names(base_p)))
  rp <- setNames(rep(0, length(types)), types); rp[names(run_p)]  <- as.numeric(run_p)
  bp <- setNames(rep(0, length(types)), types); bp[names(base_p)] <- as.numeric(base_p)
  delta <- round(rp - bp, 2)
  flag  <- abs(delta) > DRIFT
  max_drift <- max(max_drift, max(abs(delta)))
  any_flag <- any_flag || any(flag)
  df <- data.frame(anchor = s, cell_type = types, run_pct = rp, baseline_pct = bp,
                   delta_pp = delta, flag = ifelse(flag, "DRIFT", ""), row.names = NULL)
  rows[[s]] <- df
  message("-- ", s, " (run n=", sum(cells$sample == s), ") --")
  print(df[order(-abs(df$delta_pp)), ], row.names = FALSE)
  fl <- df$cell_type[flag]
  message("   max |drift|: ", max(abs(delta)), " pp", if (length(fl)) paste0("  FLAGGED: ", paste(fl, collapse = ", ")) else "  (all within tolerance)")
  message("")
}
conc <- do.call(rbind, rows)
if (!is.null(conc)) write.csv(conc, file.path(RUN_DIR, "benchmark", "concordance.csv"), row.names = FALSE)

# ---- per-batch whole-blood / sort-status readout ----------------------------
wb <- c("Neutrophil", "RBC", "Platelet")
cells$batch <- ifelse(cells$sample %in% ANCHORS, "anchor", "main")
message("=== whole-blood signature per sample (%) — high => presort/whole blood ===")
sig <- do.call(rbind, lapply(sort(unique(cells$sample)), function(s) {
  v <- cells$cell_type_ref[cells$sample == s]
  data.frame(sample = s, batch = unique(cells$batch[cells$sample == s]), n = length(v),
             Neutrophil = round(100*mean(v=="Neutrophil"),1),
             RBC = round(100*mean(v=="RBC"),1), Platelet = round(100*mean(v=="Platelet"),1))
}))
print(sig, row.names = FALSE)
write.csv(sig, file.path(RUN_DIR, "benchmark", "wholeblood_signature.csv"), row.names = FALSE)

# ---- markdown report --------------------------------------------------------
verdict <- if (any_flag) {
  sprintf("DRIFT DETECTED (max %.1f pp > %.0f pp)", max_drift, DRIFT)
} else {
  sprintf("PASS (max anchor drift %.1f pp <= %.0f pp)", max_drift, DRIFT)
}
# Per-sample collection method (terminal cardiac puncture vs peripheral vein draw). Drives the
# leukogram: a terminal central bleed avoids restraint stress neutrophilia (and can pick up
# splenic/central pooled lymphocytes), so such samples read low-neutrophil / high-lymphoid.
COLLECTION <- c(Aksh1 = "cardiac puncture (terminal/sacrifice)")
coll_of <- function(s) if (s %in% names(COLLECTION)) COLLECTION[[s]] else "peripheral vein draw"
coll_lines <- vapply(sort(unique(cells$sample)), function(s)
  sprintf("- **%s** — %s", s, coll_of(s)), character(1))

md <- c(sprintf("# Benchmark concordance — %s", basename(RUN_DIR)),
        sprintf("Model: `%s`", basename(MODEL)),
        sprintf("Anchors: %s | drift flag: %g pp", paste(ANCHORS, collapse=", "), DRIFT),
        sprintf("\n**Verdict: %s**\n", verdict),
        "Anchor proportions are run-INDEPENDENT (frozen reference); identical input + same model",
        "should reproduce across runs. Drift = pipeline/run artifact, not biology.\n",
        "See concordance.csv and wholeblood_signature.csv.\n",
        "## Collection method", coll_lines, "",
        "## Caveats (biological interpretation)",
        "- **Aksh1 is a methodological outlier, not a biological baseline.** It was obtained by",
        "  terminal cardiac puncture, the others by conscious peripheral vein draw. Aksh1 has the",
        "  lowest neutrophils (~10%) and the highest lymphoid fractions (B cell ~15%, CD4 T ~29%).",
        "  Terminal central bleeds avoid restraint stress neutrophilia and may include splenic/central",
        "  pooled lymphocytes, so this profile reflects HOW it was collected. Aksh1 stays valid as a",
        "  reproducibility anchor (labels reproduce across runs); ES332 (vein) is the representative anchor.",
        "- **B cells (0.1–7% in vein samples) are normal-to-low for neutrophil-dominated whole blood.**",
        "  Proportions are a closed sum: high neutrophils suppress the lymphoid %. B cells are well",
        "  captured in droplet scRNA (trustworthy as relative values), UNLIKE neutrophils which suffer",
        "  granulocyte dropout — so a de-novo 0% neutrophil result is the artifact, not the high values.",
        "- **Draw volume does not change proportions** (fixed cell loading); collection site/stress does.",
        "- **Whole-blood signature (high neutrophil + RBC + platelet) => presort / whole blood.** High",
        "  neutrophils in captive bats are biologically credible (neutrophil-dominant adult pteropodids;",
        "  amplified by captivity + handling stress). See docs/bat_neutrophil_literature.md.")
writeLines(md, file.path(RUN_DIR, "benchmark", "benchmark_report.md"))
message("\n=== ", verdict, " ===")
message("  report: ", file.path(RUN_DIR, "benchmark", "benchmark_report.md"))
