#!/usr/bin/env Rscript
# test_regressions.R - guards against specific past bugs reintroducing.
# Run:  Rscript pipeline/tests/test_regressions.R   (exits 1 on any failure)
# Moved out of validate_config.R, which now covers config invariants only.

.dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  d <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(d) > 0) dirname(normalizePath(d[1])) else "."
})
.pipeline <- dirname(.dir)   # tests/ -> pipeline/
fails <- character(0)

# Regression T2: 10_rarefaction.R must not hardcode ES03_newkit as ground truth.
f <- file.path(.pipeline, "10_rarefaction.R")
if (file.exists(f) &&
    any(grepl('ref_sample\\s*<-\\s*"ES03_newkit"', readLines(f, warn = FALSE))))
  fails <- c(fails, "10_rarefaction.R hardcodes ES03_newkit as ground truth")

# Regression T3: run_pipeline.sh must export SCRNA_RESULTS_DIR. Without it, bash and
# config.R derive the run dir independently and disagree for >4 samples (bash collapses
# to <first>_<N>samples, config.R does not) — the R steps then write to one directory
# while the 08b/08c paths bash hands them point at another, silently skipping both.
f <- file.path(.pipeline, "run_pipeline.sh")
if (file.exists(f) &&
    !any(grepl("^\\s*export SCRNA_RESULTS_DIR=", readLines(f, warn = FALSE))))
  fails <- c(fails, "run_pipeline.sh does not export SCRNA_RESULTS_DIR")

# Regression T4: setup_env.sh must not pin r-base below 4.4. config.R:17 calls `%||%`
# hundreds of lines before defining it, so it depends on base R's `%||%` (added in
# 4.4.0). On 4.3.x that line errors and nothing in the pipeline loads.
f <- file.path(.pipeline, "setup_env.sh")
if (file.exists(f)) {
  .pin <- grep("^\\s*r-base=", readLines(f, warn = FALSE), value = TRUE)
  .ver <- sub("^\\s*r-base=([0-9.]+).*$", "\\1", .pin)
  if (length(.ver) && any(package_version(.ver) < package_version("4.4.0")))
    fails <- c(fails, paste0("setup_env.sh pins r-base=", .ver, " (< 4.4.0; config.R needs base %||%)"))
}

# Regression T5: config.R's options(error=) handler must not dispatch conditionMessage()
# on an unguarded .Last.error. Under Rscript .Last.error is often NULL; the failed
# dispatch overwrites the error buffer, so the geterrmessage() fallback returns the
# dispatch error and every real error is reported as "no applicable method for
# 'conditionMessage' applied to an object of class NULL". This masked the step-04
# future.globals.maxSize overflow on the 8-sample bat run.
f <- file.path(.pipeline, "config.R")
if (file.exists(f)) {
  .l <- readLines(f, warn = FALSE)
  .h <- paste(.l[seq_len(min(30L, length(.l)))], collapse = "\n")
  if (grepl("conditionMessage", .h, fixed = TRUE) &&
      !grepl('inherits(last, "condition")', .h, fixed = TRUE))
    fails <- c(fails, "config.R error handler calls conditionMessage() without an inherits() guard")
}

# Regression T6: the bat (whole blood) Neutrophil panel must not contain MPO or ELANE.
# They are promyelocyte/marrow azurophilic-granule genes, absent from mature circulating
# neutrophils: measured 0.0% and 0.2% detection in the 8-sample ES49 cohort (92,864
# cells, 2026-07-30). scType scores by absolute mean, so inert genes dilute the panel and
# suppressed the neutrophil call entirely (0.0% reported; 15,748 real neutrophils were
# labelled CD14+ Mono). They were originally added on a reviewer's "conserved across
# mammals" argument — true of the genome, false of expression in blood.
f <- file.path(.pipeline, "config_species_bat.R")
if (file.exists(f)) {
  .l <- readLines(f, warn = FALSE)
  .n <- grep('^\\s*MARKERS\\$Neutrophil\\s*<-', .l, value = TRUE)
  if (length(.n) && any(grepl('"MPO"|"ELANE"', .n[1])))
    fails <- c(fails, "bat Neutrophil panel re-introduces MPO/ELANE (0.0% expressed in mature bat blood neutrophils)")
}

# Regression T7: 05_annotate.R must not assume a blood MARKERS layout when building the
# scType gene sets. It used to construct .gs_pos unconditionally from named slots
# (MARKERS$T_pan, MARKERS$CD4_T, ...). A non-blood overlay such as SCRNA_SPECIES=cm
# replaces MARKERS wholesale with cell-type-named entries, leaving every slot NULL, so
# .gs_pos collapsed to an empty list and step 05 died with the unhelpful "attempt to set
# 'colnames' on an object with less than two dimensions".
f <- file.path(.pipeline, "05_annotate.R")
if (file.exists(f)) {
  .l <- paste(readLines(f, warn = FALSE), collapse = "\n")
  if (grepl("MARKERS\\$T_pan", .l) && !grepl(".blood_slots", .l, fixed = TRUE))
    fails <- c(fails, "05_annotate.R builds scType gene sets from blood MARKERS slots with no non-blood fallback")
}

# Regression T8: a curated CLUSTER_CELLTYPE_MAP guarded on a QC value must live in the
# species overlay, not in config.R. config.R defines the map blocks well before it
# sources the overlays, so a guard term reading QC$... there sees the BASE value and the
# map leaks into a run whose QC (and therefore cluster numbering) differs. Observed:
# the cm map guarded on QC$min_features == 200 still applied to a min_features 1000 run.
# Two ways to fail this: put a QC-guarded map back in config.R, or drop the QC term from
# the overlay's guard.
f <- file.path(.pipeline, "config.R")
if (file.exists(f)) {
  .l <- readLines(f, warn = FALSE)
  .map_start <- grep("^\\s*CLUSTER_CELLTYPE_MAP\\s*<-\\s*c\\(", .l)
  .src <- grep("source\\(file\\.path\\(PIPELINE_DIR", .l)
  if (length(.map_start) && length(.src) &&
      any(.map_start < max(.src)) &&
      any(grepl("QC\\$", .l[seq_len(max(.src))]) &
          seq_len(max(.src)) %in% unlist(lapply(.map_start, function(i) max(1, i - 12):i))))
    fails <- c(fails, "config.R has a QC-guarded CLUSTER_CELLTYPE_MAP before the species overlay is sourced (QC is not final there)")
}
f <- file.path(.pipeline, "config_species_cm.R")
if (file.exists(f)) {
  # Inspect the `if (...)` CODE line, not the whole file: the explanatory comment above
  # the map also contains the string "QC$min_features == 200", so a whole-file grep
  # passes even after the guard term is deleted from the condition.
  .l <- readLines(f, warn = FALSE)
  .code <- .l[!grepl("^\\s*#", .l)]
  if (any(grepl("CLUSTER_CELLTYPE_MAP <- c(", .code, fixed = TRUE)) &&
      !any(grepl("^\\s*if \\(.*min_features", .code)))
    fails <- c(fails, "cm CLUSTER_CELLTYPE_MAP guard is missing its QC$min_features term (cluster numbering depends on the QC floor)")
}

# Regression T9: 08b_html_report.R must take CELLTYPE_COLORS and QC thresholds from config.R,
# not local copies. The copies drifted silently: cm cell types (Pluripotent, Epicardial,
# Proepicardial, ...) got fallback colours in the HTML while the PDFs used the palette, and the
# QC scatters drew human thresholds on a cm run (min_features 200 vs the overlay's 1000).
f <- file.path(.pipeline, "08b_html_report.R")
if (file.exists(f)) {
  .code <- grep("^\\s*#", readLines(f, warn = FALSE), value = TRUE, invert = TRUE)
  if (any(grepl("CELLTYPE_COLORS\\s*<-|QC_THRESH\\s*<-\\s*list\\(", .code)) ||
      !any(grepl('source\\(.*"config\\.R"', .code)))
    fails <- c(fails, "08b_html_report.R re-declares CELLTYPE_COLORS/QC_THRESH instead of sourcing config.R")
}

# Regression T10: .combine_pdfs() is defined only in pdf_helpers.R (sourced by config.R). The
# bat_wing project steps 11-14 each carried a private copy that rasterised every page to a
# 150-dpi image via magick, shadowing the shared helper, which merges PDFs losslessly.
.def <- Filter(function(f) basename(f) != "pdf_helpers.R" &&
                 any(grepl("^\\s*\\.combine_pdfs\\s*<-\\s*function", readLines(f, warn = FALSE))),
               list.files(.pipeline, pattern = "\\.R$", recursive = TRUE, full.names = TRUE))
if (length(.def))
  fails <- c(fails, paste0(".combine_pdfs() redefined outside pdf_helpers.R: ", paste(basename(unlist(.def)), collapse = ", ")))

# Regression T11: 06_visualize.R must not sort() the sample list. It used to do
#   SAMPLE_NAMES <- sort(unique(as.character(merged$sample)))
# which orders samples alphabetically and silently scrambles every figure that iterates them:
# umap_split_by_sample.pdf, the composition panels and the proportion bars. On the PBMC loading
# series PBMC_5K plotted last (after PBMC_30K_S_K); a D0..D30 time course would run D0 after D11.
# Order must follow the command-line order (config SAMPLE_NAMES) or the object's factor levels.
f <- file.path(.pipeline, "06_visualize.R")
if (file.exists(f)) {
  .code <- grep("^\\s*#", readLines(f, warn = FALSE), value = TRUE, invert = TRUE)
  if (any(grepl("SAMPLE_NAMES\\s*<-\\s*sort\\(", .code)))
    fails <- c(fails, "06_visualize.R sorts SAMPLE_NAMES alphabetically (scrambles per-sample figure order)")
}

if (length(fails) > 0) {
  for (x in fails) message("FAIL: ", x)
  quit(status = 1)
}
message("All regression checks passed.")
