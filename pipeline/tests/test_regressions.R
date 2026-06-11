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

if (length(fails) > 0) {
  for (x in fails) message("FAIL: ", x)
  quit(status = 1)
}
message("All regression checks passed.")
