#!/usr/bin/env Rscript
# validate_config.R - Runs before any pipeline step. Exit code 1 on failure.
# Robust path resolution: works whether sourced OR run as top-level Rscript.
.this_dir <- local({
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(f)) dirname(f)
  else {
    a <- commandArgs(trailingOnly = FALSE)
    d <- sub("--file=", "", a[grep("--file=", a)])
    if (length(d) > 0) dirname(normalizePath(d)) else "."
  }
})
source(file.path(.this_dir, "config.R"))
errors <- character(0)

# Check 1: CELLTYPE_COLORS covers all types in CLUSTER_CELLTYPE_MAP
if (!is.null(CLUSTER_CELLTYPE_MAP)) {
  missing_colors <- setdiff(unique(CLUSTER_CELLTYPE_MAP), names(CELLTYPE_COLORS))
  if (length(missing_colors) > 0)
    errors <- c(errors, paste0("CLUSTER_CELLTYPE_MAP types missing from CELLTYPE_COLORS: ",
                               paste(missing_colors, collapse = ", ")))
}

# Check 2: CLUSTER_CELLTYPE_MAP keys are quoted integers
if (!is.null(CLUSTER_CELLTYPE_MAP)) {
  bad_keys <- names(CLUSTER_CELLTYPE_MAP)[!grepl("^[0-9]+$", names(CLUSTER_CELLTYPE_MAP))]
  if (length(bad_keys) > 0)
    errors <- c(errors, paste0("CLUSTER_CELLTYPE_MAP keys must be cluster numbers: ",
                               paste(bad_keys, collapse = ", ")))
}

# Check 3: SAMPLE_PATHS point to existing directories
for (nm in names(SAMPLE_PATHS)) {
  if (!dir.exists(SAMPLE_PATHS[[nm]]))
    errors <- c(errors, paste0("Sample path does not exist: ", SAMPLE_PATHS[[nm]]))
}

# Regression T2: step 10 must not hardcode ES03_newkit
step10_lines <- readLines(file.path(.this_dir, "10_rarefaction.R"))
if (any(grepl('ref_sample\\s*<-\\s*"ES03_newkit"', step10_lines)))
  errors <- c(errors, "REGRESSION: 10_rarefaction.R hardcodes ES03_newkit as ground truth")

if (length(errors) > 0) {
  message("CONFIG VALIDATION FAILED:")
  for (e in errors) message("  ERROR: ", e)
  quit(status = 1)
}
message("Config validation passed.")
