# =============================================================================
# 13_cellchat.R — Cell-cell communication analysis (bat_wing mode)
# Requires: integrated_annotated.rds from step 04+05
# Outputs:  CellChat RDS per condition, chord diagrams, bubble plots,
#           differential signaling heatmap, cellchat_report.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(CellChat); library(Seurat); library(dplyr)
  library(ggplot2); library(patchwork); library(pdftools); library(magick)
})
.pipeline_dir <- {
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "."
}
source(file.path(.pipeline_dir, "config.R"))

if (length(CONDITION_LEVELS) < 2) {
  message("Only one condition — skipping 13_cellchat.R")
  quit(save = "no", status = 0)
}
dir.create(DIRS$cellchat, showWarnings = FALSE, recursive = TRUE)

.combine_pdfs <- function(paths, out) {
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) return(invisible(NULL))
  imgs <- lapply(paths, function(p) {
    n <- tryCatch(pdftools::pdf_length(p), error = function(e) 1L)
    lapply(seq_len(n), function(i) magick::image_read_pdf(p, pages = i, density = 150))
  })
  magick::image_write(magick::image_join(unlist(imgs, recursive = FALSE)), out, format = "pdf")
}

rds_path <- file.path(DIRS$integrated, "integrated_annotated.rds")
if (!file.exists(rds_path))
  stop("Missing: ", rds_path, " — run steps 04+05 first")
merged <- readRDS(rds_path)
merged$condition <- unname(SAMPLE_CONDITIONS[merged$sample])

# CellChat uses human LR database (bat gene names are human-compatible)
CellChatDB     <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB,
                            search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))

# Build one CellChat object per condition
.build_cellchat <- function(seurat_obj, condition_label) {
  message("Building CellChat for: ", condition_label)
  sub <- subset(seurat_obj, condition == condition_label)
  Idents(sub) <- "cell_type"
  keep_types  <- setdiff(unique(sub$cell_type), CONTAMINATION_TYPES)
  sub <- subset(sub, idents = keep_types)

  cc <- createCellChat(object = sub, group.by = "cell_type", assay = "RNA")
  cc@DB <- CellChatDB.use
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean", nboot = 100)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

cc_list <- lapply(CONDITION_LEVELS, function(cond) {
  tryCatch(.build_cellchat(merged, cond),
           error = function(e) { message("CellChat failed for ", cond, ": ", e$message); NULL })
})
names(cc_list) <- CONDITION_LEVELS
cc_list <- Filter(Negate(is.null), cc_list)

if (length(cc_list) == 0) {
  message("No CellChat objects built — skipping outputs")
  quit(save = "no", status = 0)
}

# Save individual objects
for (cond in names(cc_list)) {
  saveRDS(cc_list[[cond]],
          file.path(DIRS$cellchat, paste0("cellchat_", cond, ".rds")))
}

cond_a <- CONDITION_LEVELS[1]
cond_b <- CONDITION_LEVELS[2]

# =============================================================================
# PART 1: Per-condition chord diagrams (interaction count)
# =============================================================================
chord_pdfs <- character(0)
for (cond in names(cc_list)) {
  tmp <- tempfile(fileext = ".pdf")
  tryCatch({
    pdf(tmp, width = 8, height = 8)
    netVisual_circle(cc_list[[cond]]@net$count,
                     vertex.weight = groupSize(cc_list[[cond]]),
                     weight.scale = TRUE, label.edge = FALSE,
                     title.name = paste0("Interaction count — ", cond))
    dev.off()
    chord_pdfs <- c(chord_pdfs, tmp)
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    message("chord failed for ", cond, ": ", e$message)
  })
}
.combine_pdfs(chord_pdfs, file.path(DIRS$cellchat, "chord_diagrams.pdf"))

# =============================================================================
# PART 2: Bubble plot — all significant interactions
# =============================================================================
bubble_pdfs <- character(0)
for (cond in names(cc_list)) {
  tmp <- tempfile(fileext = ".pdf")
  tryCatch({
    p <- netVisual_bubble(cc_list[[cond]], remove.isolate = FALSE,
                          angle.x = 45, font.size = 8) +
      ggtitle(paste0("Interactions — ", cond))
    ggsave(tmp, p, width = 12,
           height = max(6, nrow(cc_list[[cond]]@LR$LRsig) * 0.1))
    bubble_pdfs <- c(bubble_pdfs, tmp)
  }, error = function(e) message("bubble failed for ", cond, ": ", e$message))
}
.combine_pdfs(bubble_pdfs, file.path(DIRS$cellchat, "bubble_plots.pdf"))

# =============================================================================
# PART 3: Differential signaling (merged comparison) — only if both conditions built
# =============================================================================
all_cc_pdfs <- c(chord_pdfs, bubble_pdfs)

if (length(cc_list) >= 2) {
  cc_merged <- tryCatch(
    mergeCellChat(cc_list, add.names = names(cc_list)),
    error = function(e) { message("mergeCellChat failed: ", e$message); NULL }
  )

  if (!is.null(cc_merged)) {
    # Differential interaction count/weight bar plot
    diff_pdf <- tempfile(fileext = ".pdf")
    tryCatch({
      pdf(diff_pdf, width = 10, height = 5)
      compareInteractions(cc_merged, show.legend = FALSE, group = names(cc_list))
      dev.off()
      all_cc_pdfs <- c(all_cc_pdfs, diff_pdf)
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      message("compareInteractions failed: ", e$message)
    })

    # Per-pathway strength comparison
    shared_pathways <- Reduce(intersect,
      lapply(cc_list, function(cc) cc@netP$pathways))
    top_paths <- head(shared_pathways, 12)

    for (pw in top_paths) {
      tmp <- tempfile(fileext = ".pdf")
      tryCatch({
        pdf(tmp, width = 10, height = 5)
        rankNet(cc_merged, mode = "comparison", measure = "weight",
                sources.use = NULL, targets.use = NULL,
                stacked = TRUE, do.stat = TRUE,
                color.use = c("steelblue", "tomato"),
                title = paste0("Pathway strength: ", pw))
        dev.off()
        all_cc_pdfs <- c(all_cc_pdfs, tmp)
      }, error = function(e) {
        if (dev.cur() > 1) dev.off()
        message("rankNet failed for ", pw, ": ", e$message)
      })
    }

    # Interaction summary table
    tryCatch({
      df_net <- subsetCommunication(cc_merged)
      write.csv(df_net,
                file.path(DIRS$cellchat, "interactions_summary.csv"), row.names = FALSE)
    }, error = function(e) message("Could not export interaction table: ", e$message))
  }
}

.combine_pdfs(all_cc_pdfs, file.path(DIRS$cellchat, "cellchat_report.pdf"))
message("Step 13 complete — CellChat results in: ", DIRS$cellchat)
