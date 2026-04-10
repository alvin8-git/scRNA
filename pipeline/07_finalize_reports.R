# =============================================================================
# 07_finalize_reports.R - Merge per-step partial PDFs into 5 category reports:
#   report_qc.pdf           (01_load_qc only)
#   report_doublets.pdf     (02_doublets)
#   report_individual.pdf   (01_load_qc + 03_individual)
#   report_annotation.pdf   (05_annotate)
#   report_integrated.pdf   (04_integrate + 06_visualize)
# =============================================================================
source("/data/alvin/scRNA/pipeline/config.R")

combine <- function(inputs, output) {
  existing <- inputs[file.exists(inputs)]
  if (length(existing) == 0) {
    message("  Skipped (no inputs found): ", basename(output))
    return(invisible(NULL))
  }
  .combine_pdfs(existing, output)
  message("  -> ", output)
}

message("\n=== Finalising PDF reports ===")

combine(
  c(file.path(DIRS$qc, "qc_report.pdf")),
  file.path(DIRS$reports, "report_qc.pdf")
)

combine(
  c(file.path(DIRS$doublets, "doublets_report.pdf")),
  file.path(DIRS$reports, "report_doublets.pdf")
)

combine(
  c(file.path(DIRS$qc,         "qc_report.pdf"),
    file.path(DIRS$individual,  "individual_report.pdf")),
  file.path(DIRS$reports, "report_individual.pdf")
)

combine(
  c(file.path(DIRS$annotation, "annotation_report.pdf")),
  file.path(DIRS$reports, "report_annotation.pdf")
)

combine(
  c(file.path(DIRS$integrated, "integration_report.pdf"),
    file.path(DIRS$integrated, "visualization_report.pdf")),
  file.path(DIRS$reports, "report_integrated.pdf")
)

message("\nFinal reports in: ", DIRS$reports)
message("  report_qc.pdf")
message("  report_doublets.pdf")
message("  report_individual.pdf")
message("  report_annotation.pdf")
message("  report_integrated.pdf")
message("\n07_finalize_reports.R complete.")
