# =============================================================================
# 12_pathways.R — GO/KEGG pathway enrichment and GSEA per cell type
# Requires: all_DEG_combined.csv from step 11
# Outputs:  enrichment CSVs, dot plots, GSEA plots → DIRS$pathways
# Note:     Bat genome uses human gene names — org.Hs.eg.db and hsa KEGG apply
# =============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
  library(ggplot2); library(dplyr); library(pdftools); library(magick)
})

.pipeline_dir <- {
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "."
}
source(file.path(.pipeline_dir, "config.R"))

dir.create(DIRS$pathways, showWarnings = FALSE, recursive = TRUE)

.combine_pdfs <- function(paths, out) {
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) return(invisible(NULL))
  imgs <- lapply(paths, function(p) {
    n <- tryCatch(pdftools::pdf_length(p), error = function(e) 1L)
    lapply(seq_len(n), function(i) magick::image_read_pdf(p, pages = i, density = 150))
  })
  magick::image_write(magick::image_join(unlist(imgs, recursive = FALSE)),
                      out, format = "pdf")
}

de_path <- file.path(DIRS$differential, "all_DEG_combined.csv")
if (!file.exists(de_path)) stop("Missing: ", de_path, " — run step 11 first")
all_de <- read.csv(de_path)

# Silence clusterProfiler verbose output
options(clusterProfiler.download.method = "wget")

# Helper: gene symbols → Entrez IDs (bat uses human gene names)
.to_entrez <- function(genes) {
  res <- suppressMessages(
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db, drop = TRUE)
  )
  res$ENTREZID
}

cell_types  <- unique(all_de$cell_type)
all_enrich  <- list()
all_pdfs    <- character(0)

for (ct in cell_types) {
  de_ct <- all_de %>% dplyr::filter(cell_type == ct, p_val_adj < 0.05)
  if (nrow(de_ct) < 5) {
    message("Skipping ", ct, " — fewer than 5 significant DEGs")
    next
  }
  message("Pathways: ", ct, " (", nrow(de_ct), " sig DEGs)")

  up_genes   <- de_ct %>% dplyr::filter(avg_log2FC >  0.5) %>% dplyr::pull(gene)
  down_genes <- de_ct %>% dplyr::filter(avg_log2FC < -0.5) %>% dplyr::pull(gene)
  universe   <- all_de %>% dplyr::filter(cell_type == ct) %>% dplyr::pull(gene) %>% unique()

  ct_label <- gsub("[/ ()]", "_", ct)
  ct_dir   <- file.path(DIRS$pathways, ct_label)
  dir.create(ct_dir, showWarnings = FALSE)

  # --- GO enrichment (up-regulated) -----------------------------------------
  if (length(up_genes) >= 5) {
    ego_up <- tryCatch(
      suppressMessages(
        enrichGO(gene          = .to_entrez(up_genes),
                 universe      = .to_entrez(universe),
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.2,
                 readable      = TRUE)
      ),
      error = function(e) { message("  GO up failed: ", e$message); NULL }
    )
    if (!is.null(ego_up) && nrow(ego_up) > 0) {
      write.csv(as.data.frame(ego_up),
                file.path(ct_dir, "GO_BP_up.csv"), row.names = FALSE)
      p        <- dotplot(ego_up, showCategory = 15) +
        labs(title = paste0(ct, " — GO BP (up in recovering)")) +
        theme_classic(base_size = 9)
      pdf_file <- file.path(ct_dir, "GO_BP_up.pdf")
      ggsave(pdf_file, p, width = 8, height = 6)
      all_pdfs <- c(all_pdfs, pdf_file)
      all_enrich[[paste0(ct, "_GO_up")]] <- as.data.frame(ego_up) %>%
        dplyr::mutate(cell_type = ct, direction = "up")
    }
  }

  # --- GO enrichment (down-regulated) ----------------------------------------
  if (length(down_genes) >= 5) {
    ego_dn <- tryCatch(
      suppressMessages(
        enrichGO(gene          = .to_entrez(down_genes),
                 universe      = .to_entrez(universe),
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.2,
                 readable      = TRUE)
      ),
      error = function(e) { message("  GO down failed: ", e$message); NULL }
    )
    if (!is.null(ego_dn) && nrow(ego_dn) > 0) {
      write.csv(as.data.frame(ego_dn),
                file.path(ct_dir, "GO_BP_down.csv"), row.names = FALSE)
      p        <- dotplot(ego_dn, showCategory = 15) +
        labs(title = paste0(ct, " — GO BP (up in healthy)")) +
        theme_classic(base_size = 9)
      pdf_file <- file.path(ct_dir, "GO_BP_down.pdf")
      ggsave(pdf_file, p, width = 8, height = 6)
      all_pdfs <- c(all_pdfs, pdf_file)
    }
  }

  # --- KEGG enrichment (up-regulated) ----------------------------------------
  if (length(up_genes) >= 5) {
    ekegg <- tryCatch(
      suppressMessages(
        enrichKEGG(gene          = .to_entrez(up_genes),
                   universe      = .to_entrez(universe),
                   organism      = "hsa",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)
      ),
      error = function(e) { message("  KEGG failed: ", e$message); NULL }
    )
    if (!is.null(ekegg) && nrow(ekegg) > 0) {
      ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      write.csv(as.data.frame(ekegg),
                file.path(ct_dir, "KEGG_up.csv"), row.names = FALSE)
      p        <- dotplot(ekegg, showCategory = 15) +
        labs(title = paste0(ct, " — KEGG (up in recovering)")) +
        theme_classic(base_size = 9)
      pdf_file <- file.path(ct_dir, "KEGG_up.pdf")
      ggsave(pdf_file, p, width = 8, height = 6)
      all_pdfs <- c(all_pdfs, pdf_file)
    }
  }

  # --- GSEA on full ranked gene list -----------------------------------------
  de_all_ct <- all_de %>% dplyr::filter(cell_type == ct) %>%
    dplyr::arrange(dplyr::desc(avg_log2FC))
  ranked_sym <- setNames(de_all_ct$avg_log2FC, de_all_ct$gene)

  map_df <- tryCatch(
    suppressMessages(
      bitr(names(ranked_sym), fromType = "SYMBOL", toType = "ENTREZID",
           OrgDb = org.Hs.eg.db, drop = TRUE)
    ),
    error = function(e) NULL
  )

  if (!is.null(map_df) && nrow(map_df) >= 20) {
    ranked_entrez <- setNames(ranked_sym[map_df$SYMBOL], map_df$ENTREZID)
    ranked_entrez <- sort(ranked_entrez, decreasing = TRUE)

    gsea_res <- tryCatch(
      suppressMessages(
        gseGO(geneList     = ranked_entrez,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 10,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
      ),
      error = function(e) { message("  GSEA failed: ", e$message); NULL }
    )
    if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
      write.csv(as.data.frame(gsea_res),
                file.path(ct_dir, "GSEA_GO_BP.csv"), row.names = FALSE)
      top_ids <- head(
        gsea_res@result$ID[order(abs(gsea_res@result$NES), decreasing = TRUE)], 5
      )
      pdf_gsea <- tryCatch({
        p        <- gseaplot2(gsea_res, geneSetID = top_ids,
                              title = paste0(ct, " — GSEA top pathways"))
        pdf_file <- file.path(ct_dir, "GSEA_GO_BP.pdf")
        ggsave(pdf_file, p, width = 10, height = 6)
        pdf_file
      }, error = function(e) { message("  gseaplot2 failed: ", e$message); NA_character_ })
      if (!is.na(pdf_gsea)) all_pdfs <- c(all_pdfs, pdf_gsea)
    }
  }
}

# Combined summary: top 3 pathways per cell type (by adjusted p)
if (length(all_enrich) > 0) {
  summary_df <- dplyr::bind_rows(all_enrich) %>%
    dplyr::group_by(cell_type, direction) %>%
    dplyr::slice_min(p.adjust, n = 3) %>%
    dplyr::ungroup()
  write.csv(summary_df,
            file.path(DIRS$pathways, "pathway_summary.csv"), row.names = FALSE)
  message("Saved pathway_summary.csv (", nrow(summary_df), " rows)")
}

# Combined PDF report
if (length(all_pdfs) > 0) {
  .combine_pdfs(all_pdfs, file.path(DIRS$pathways, "pathways_report.pdf"))
  message("Saved pathways_report.pdf")
}

message("Step 12 complete — pathway results in: ", DIRS$pathways)
