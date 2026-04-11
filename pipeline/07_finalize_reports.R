# =============================================================================
# 07_finalize_reports.R - Merge per-step partial PDFs into 5 category reports
#   plus an Overall_report combining curated plots from all stages.
#
#   01-QC_report.pdf
#   02-Doublet_report.pdf
#   03-Individual_report.pdf
#   04-Annotation_report.pdf
#   05-Integrated_report.pdf
#   Overall_report.pdf
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
  file.path(DIRS$reports, "01-QC_report.pdf")
)

combine(
  c(file.path(DIRS$doublets, "doublets_report.pdf")),
  file.path(DIRS$reports, "02-Doublet_report.pdf")
)

combine(
  c(file.path(DIRS$qc,        "qc_report.pdf"),
    file.path(DIRS$individual, "individual_report.pdf")),
  file.path(DIRS$reports, "03-Individual_report.pdf")
)

combine(
  c(file.path(DIRS$annotation, "annotation_report.pdf")),
  file.path(DIRS$reports, "04-Annotation_report.pdf")
)

combine(
  c(file.path(DIRS$integrated, "integration_report.pdf"),
    file.path(DIRS$integrated, "visualization_report.pdf")),
  file.path(DIRS$reports, "05-Integrated_report.pdf")
)

# =============================================================================
# Overall_report.pdf  -  curated cross-stage summary
# =============================================================================
overall_inputs <- c(
  # QC (all pages)
  file.path(DIRS$qc, "qc_report.pdf"),
  # Doublets (all pages)
  file.path(DIRS$doublets, "doublets_report.pdf"),
  # Annotation: SingleR score heatmap only
  file.path(DIRS$annotation, "singler_scores_heatmap.pdf"),
  # Individual: HVG + top-5-markers dot plot per sample
  unlist(lapply(SAMPLE_NAMES, function(nm) c(
    file.path(DIRS$individual, nm, paste0(nm, "_hvg.pdf")),
    file.path(DIRS$individual, nm, paste0(nm, "_dotplot_markers.pdf"))
  ))),
  # Integrated
  file.path(DIRS$integrated, "harmony_before_after.pdf"),
  file.path(DIRS$integrated, "integrated_umap_sample.pdf"),
  file.path(DIRS$integrated, "integrated_umap_celltype.pdf"),
  file.path(DIRS$integrated, "umap_split_by_sample.pdf"),
  file.path(DIRS$integrated, "umap_triptych.pdf"),
  file.path(DIRS$integrated, "integrated_dotplot.pdf"),
  file.path(DIRS$integrated, "integrated_heatmap.pdf"),
  file.path(DIRS$integrated, "celltype_composition_combined.pdf"),
  file.path(DIRS$integrated, "violin_key_markers.pdf")
)
combine(overall_inputs, file.path(DIRS$reports, "Overall_report.pdf"))

message("\nFinal reports in: ", DIRS$reports)
message("  01-QC_report.pdf")
message("  02-Doublet_report.pdf")
message("  03-Individual_report.pdf")
message("  04-Annotation_report.pdf")
message("  05-Integrated_report.pdf")
message("  Overall_report.pdf")
message("\n07_finalize_reports.R complete.")
