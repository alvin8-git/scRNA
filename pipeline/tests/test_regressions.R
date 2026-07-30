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

if (length(fails) > 0) {
  for (x in fails) message("FAIL: ", x)
  quit(status = 1)
}
message("All regression checks passed.")
