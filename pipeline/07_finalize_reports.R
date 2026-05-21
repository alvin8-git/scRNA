# =============================================================================
# 07_finalize_reports.R - Merge per-step partial PDFs into 5 category reports
#   plus an Overall_report.pdf where every page is normalized to A4 size with
#   a bold title banner and description / Good-Bad caption.
#
#   01-QC_report.pdf
#   02-Doublet_report.pdf
#   03-Individual_report.pdf
#   04-Annotation_report.pdf
#   05-Integrated_report.pdf
#   Overall_report.pdf
# =============================================================================
source("/data/alvin/scRNA/pipeline/config.R")

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0) x else y

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
  c(file.path(DIRS$annotation, "annotation_report.pdf"),
    file.path(DIRS$annotation, "contamination_summary.pdf")),
  file.path(DIRS$reports, "04-Annotation_report.pdf")
)

combine(
  c(file.path(DIRS$integrated, "integration_report.pdf"),
    file.path(DIRS$integrated, "visualization_report.pdf")),
  file.path(DIRS$reports, "05-Integrated_report.pdf")
)

# =============================================================================
# Overall_report.pdf — every page normalized to A4 portrait or landscape,
#   with a bold title banner at top and a description/interpretation strip
#   at the bottom.  Requires the 'pdftools' R package to rasterize source PDFs.
# =============================================================================

make_overall_report <- function(output_path) {
  suppressPackageStartupMessages({
    library(cowplot)
    library(ggplot2)
  })

  if (!requireNamespace("pdftools", quietly = TRUE)) {
    message("  pdftools not available — using simple merge for Overall_report.")
    message("  Install with: install.packages('pdftools')")
    return(FALSE)
  }

  # A4 dimensions in inches
  A4P <- c(w = 8.27,  h = 11.69)   # portrait
  A4L <- c(w = 11.69, h = 8.27)    # landscape

  # ------------------------------------------------------------------
  # Render one page of a PDF to a nativeRaster (returns NULL on error)
  # ------------------------------------------------------------------
  .render <- function(path, page = 1L, dpi = 150L) {
    if (is.null(path) || !file.exists(path)) return(NULL)
    n <- tryCatch(pdftools::pdf_length(path), error = function(e) 0L)
    if (n < page) return(NULL)
    tryCatch(
      pdftools::pdf_render_page(path, page = page, dpi = dpi, numeric = FALSE),
      error = function(e) NULL
    )
  }

  # ------------------------------------------------------------------
  # Build a single A4 page: image grid + title banner + caption strip.
  # imgs   : list of nativeRaster objects (NULLs are dropped)
  # ncols  : columns in the image grid (rows calculated automatically)
  # Returns list(plot, w, h) or NULL if no images
  # ------------------------------------------------------------------
  .build_page <- function(imgs, title = NULL, caption = NULL,
                           landscape = FALSE, ncols = 1L) {
    imgs <- Filter(Negate(is.null), imgs)
    if (length(imgs) == 0) return(NULL)

    dims <- if (landscape) A4L else A4P
    has_title   <- !is.null(title)   && nchar(trimws(title))   > 0
    has_caption <- !is.null(caption) && nchar(trimws(caption)) > 0

    th <- if (has_title)   0.055 else 0
    ch <- if (has_caption) 0.14  else 0
    ih <- 1 - th - ch

    # Image grid
    img_plots <- lapply(imgs, function(img) ggdraw() + draw_image(img))
    nrows <- ceiling(length(img_plots) / ncols)
    img_grid <- if (length(img_plots) == 1L) {
      img_plots[[1]]
    } else {
      plot_grid(plotlist = img_plots, ncol = ncols, nrow = nrows)
    }

    # Row list and relative heights
    rows <- list(); rel_h <- numeric(0)

    if (has_title) {
      rows[[length(rows) + 1L]] <- ggdraw() +
        draw_label(title, fontface = "bold", size = 13,
                   x = 0.03, y = 0.5, hjust = 0, vjust = 0.5)
      rel_h <- c(rel_h, th)
    }
    rows[[length(rows) + 1L]] <- img_grid
    rel_h <- c(rel_h, ih)
    if (has_caption) {
      rows[[length(rows) + 1L]] <- ggdraw() +
        draw_label(caption, size = 7.5, color = "#444444",
                   x = 0.03, y = 0.92, hjust = 0, vjust = 1, lineheight = 1.25)
      rel_h <- c(rel_h, ch)
    }

    pg <- if (length(rows) == 1L) rows[[1L]] else
      plot_grid(plotlist = rows, ncol = 1L, rel_heights = rel_h)

    list(plot = pg, w = unname(dims["w"]), h = unname(dims["h"]))
  }

  # Write one page spec to a temp PDF; returns path or NULL
  .save_page <- function(spec) {
    if (is.null(spec)) return(NULL)
    tf <- tempfile(fileext = ".pdf")
    pdf(tf, width = spec$w, height = spec$h)
    tryCatch(
      { print(spec$plot); dev.off(); tf },
      error = function(e) { try(dev.off(), silent = TRUE); NULL }
    )
  }

  # Convenience: build + save in one call; appends path to `pages` in parent
  pages <- character(0)
  .add <- function(imgs, title, caption = NULL, landscape = FALSE, ncols = 1L) {
    cap <- caption %||% .get_caption(title)
    spec <- .build_page(imgs, title, cap, landscape, ncols)
    tf <- .save_page(spec)
    if (!is.null(tf)) pages <<- c(pages, tf)
  }

  # Combined caption helper: paste two captions, skip empty ones
  .two_caps <- function(key1, key2) {
    c1 <- .get_caption(key1) %||% ""
    c2 <- .get_caption(key2) %||% ""
    trimws(paste(c1, c2, sep = if (nchar(c1) > 0 && nchar(c2) > 0) "\n" else ""))
  }

  # Paired-image page: title | img1 | cap1 (directly below) | img2 | cap2
  # Used when each image needs its own caption directly beneath it.
  .add_paired <- function(img1, img2, title, key1, key2, landscape = FALSE) {
    img1 <- if (is.character(img1)) .render(img1) else img1
    img2 <- if (is.character(img2)) .render(img2) else img2
    if (is.null(img1) && is.null(img2)) return(invisible(NULL))

    dims <- if (landscape) A4L else A4P
    cap1 <- .get_caption(key1) %||% ""
    cap2 <- .get_caption(key2) %||% ""

    has_title <- !is.null(title) && nchar(trimws(title)) > 0
    has_c1    <- nchar(trimws(cap1)) > 0
    has_c2    <- nchar(trimws(cap2)) > 0 && !is.null(img2)

    rows <- list(); rel_h <- numeric(0)

    if (has_title) {
      rows[[length(rows)+1L]] <- ggdraw() +
        draw_label(title, fontface = "bold", size = 13,
                   x = 0.03, y = 0.5, hjust = 0, vjust = 0.5)
      rel_h <- c(rel_h, 0.05)
    }

    if (!is.null(img1)) {
      rows[[length(rows)+1L]] <- ggdraw() + draw_image(img1)
      rel_h <- c(rel_h, if (is.null(img2)) 0.88 else 0.40)
      if (has_c1) {
        rows[[length(rows)+1L]] <- ggdraw() +
          draw_label(cap1, size = 7, color = "#444444",
                     x = 0.03, y = 0.92, hjust = 0, vjust = 1, lineheight = 1.2)
        rel_h <- c(rel_h, 0.07)
      }
    }

    if (!is.null(img2)) {
      rows[[length(rows)+1L]] <- ggdraw() + draw_image(img2)
      rel_h <- c(rel_h, 0.38)
      if (has_c2) {
        rows[[length(rows)+1L]] <- ggdraw() +
          draw_label(cap2, size = 7, color = "#444444",
                     x = 0.03, y = 0.92, hjust = 0, vjust = 1, lineheight = 1.2)
        rel_h <- c(rel_h, 0.10)
      }
    }

    rel_h <- rel_h / sum(rel_h)  # normalise to 1
    pg <- plot_grid(plotlist = rows, ncol = 1L, rel_heights = rel_h)
    tf <- .save_page(list(plot = pg, w = unname(dims["w"]), h = unname(dims["h"])))
    if (!is.null(tf)) pages <<- c(pages, tf)
  }

  message("  Building Overall_report.pdf pages (A4 normalized)...")

  # ------------------------------------------------------------------ #
  # 1. QC: violin + scatter stacked — one portrait page per sample      #
  # ------------------------------------------------------------------ #
  for (nm in SAMPLE_NAMES) {
    .add_paired(
      img1      = file.path(DIRS$qc, paste0(nm, "_violin_qc.pdf")),
      img2      = file.path(DIRS$qc, paste0(nm, "_scatter_qc.pdf")),
      title     = paste0(nm, " — QC Violin & Scatter Plots"),
      key1      = "QC Violin",
      key2      = "QC Scatter",
      landscape = FALSE
    )
  }

  # ------------------------------------------------------------------ #
  # 2. Doublet UMAP + Score Histogram stacked — one portrait per sample #
  # ------------------------------------------------------------------ #
  for (nm in SAMPLE_NAMES) {
    .add_paired(
      img1      = file.path(DIRS$doublets, paste0(nm, "_doublet_umap.pdf")),
      img2      = file.path(DIRS$doublets, paste0(nm, "_doublet_score_hist.pdf")),
      title     = paste0(nm, " — Doublet Detection UMAP & Score Distribution"),
      key1      = "Doublet.*UMAP",
      key2      = "Doublet.*Hist|Score Dist",
      landscape = FALSE
    )
  }

  # ------------------------------------------------------------------ #
  # 3. SingleR score heatmap — landscape                                #
  # ------------------------------------------------------------------ #
  .add(
    list(.render(file.path(DIRS$annotation, "singler_scores_heatmap.pdf"))),
    title     = "SingleR Reference Score Heatmap",
    landscape = TRUE
  )

  # ------------------------------------------------------------------ #
  # 3b. Contamination & rare cell type summary — portrait               #
  # ------------------------------------------------------------------ #
  {
    contam_pdf <- file.path(DIRS$annotation, "contamination_summary.pdf")
    if (file.exists(contam_pdf)) {
      n_contam <- tryCatch(pdftools::pdf_length(contam_pdf), error = function(e) 1L)
      for (pg in seq_len(n_contam)) {
        .lbl <- if (n_contam == 1L) "Contamination & Rare Cell Type Summary" else
          sprintf("Contamination & Rare Cell Type Summary (%d/%d)", pg, n_contam)
        .add(list(.render(contam_pdf, page = pg)), title = .lbl, landscape = FALSE)
      }
    }
  }

  # ------------------------------------------------------------------ #
  # 4. HVG + Top-5 Markers per sample — each caption below its image   #
  # ------------------------------------------------------------------ #
  for (nm in SAMPLE_NAMES) {
    .add_paired(
      img1      = file.path(DIRS$individual, nm, paste0(nm, "_hvg.pdf")),
      img2      = file.path(DIRS$individual, nm, paste0(nm, "_dotplot_markers.pdf")),
      title     = paste0(nm, " — Highly Variable Genes & Top 5 Markers per Cluster"),
      key1      = "Variable Genes|Highly Variable",
      key2      = "Top.*Marker.*Dot|Dot.*Top.*Marker",
      landscape = FALSE
    )
  }

  # ------------------------------------------------------------------ #
  # 5. Harmony before / after — landscape                               #
  # ------------------------------------------------------------------ #
  .add(
    list(.render(file.path(DIRS$integrated, "harmony_before_after.pdf"))),
    title     = "Harmony Batch Correction — UMAP Before & After",
    landscape = TRUE
  )

  # ------------------------------------------------------------------ #
  # 6. UMAP by sample (left) + UMAP cell type (right) — landscape 2-up #
  # ------------------------------------------------------------------ #
  .add(
    list(
      .render(file.path(DIRS$integrated, "integrated_umap_sample.pdf")),
      .render(file.path(DIRS$integrated, "integrated_umap_celltype.pdf"))
    ),
    title     = "Integrated UMAP — Sample Origin (left) & Cell Type (right)",
    landscape = TRUE, ncols = 2L
  )

  # ------------------------------------------------------------------ #
  # 7. UMAP split by sample — landscape (one page per 2 samples)       #
  # ------------------------------------------------------------------ #
  {
    split_pdf  <- file.path(DIRS$integrated, "umap_split_by_sample.pdf")
    n_split    <- tryCatch(pdftools::pdf_length(split_pdf), error = function(e) 0L)
    for (pg in seq_len(n_split)) {
      lbl <- if (n_split == 1L) "Integrated UMAP — Split by Sample" else
               sprintf("Integrated UMAP — Split by Sample (%d/%d)", pg, n_split)
      .add(list(.render(split_pdf, page = pg)), title = lbl, landscape = TRUE)
    }
  }

  # ------------------------------------------------------------------ #
  # 8. UMAP triptych — landscape                                        #
  # ------------------------------------------------------------------ #
  .add(
    list(.render(file.path(DIRS$integrated, "umap_triptych.pdf"))),
    title     = "UMAP Triptych — Clusters / Cell Types / Sample Origin",
    landscape = TRUE
  )

  # ------------------------------------------------------------------ #
  # 9. Canonical PBMC dot plot — landscape                              #
  # ------------------------------------------------------------------ #
  .add(
    list(.render(file.path(DIRS$integrated, "integrated_dotplot.pdf"))),
    title     = "Canonical PBMC Markers Dot Plot",
    landscape = TRUE
  )

  # ------------------------------------------------------------------ #
  # 10. Heatmap top-3 markers per cluster — landscape                  #
  # ------------------------------------------------------------------ #
  .add(
    list(.render(file.path(DIRS$integrated, "integrated_heatmap.pdf"))),
    title     = "Heatmap — Top 3 DE Markers per Cluster",
    landscape = TRUE
  )

  # ------------------------------------------------------------------ #
  # 11. Cell type composition — portrait (one page per ≤3 samples)    #
  # ------------------------------------------------------------------ #
  {
    comp_pdf  <- file.path(DIRS$integrated, "celltype_composition_combined.pdf")
    n_comp    <- tryCatch(pdftools::pdf_length(comp_pdf), error = function(e) 0L)
    for (pg in seq_len(n_comp)) {
      lbl <- if (n_comp == 1L) "Cell Type Composition — Proportion & Count per Sample" else
               sprintf("Cell Type Composition — Proportion & Count per Sample (%d/%d)", pg, n_comp)
      .add(list(.render(comp_pdf, page = pg)), title = lbl, landscape = FALSE)
    }
  }

  # ------------------------------------------------------------------ #
  # 12. Violin key lineage markers — landscape                          #
  # ------------------------------------------------------------------ #
  .add(
    list(.render(file.path(DIRS$integrated, "violin_key_markers.pdf"))),
    title     = "Key Lineage Marker Expression by Cell Type",
    landscape = TRUE
  )

  # Combine all temp pages into the final PDF
  valid <- Filter(file.exists, pages)
  if (length(valid) == 0) {
    message("  No pages rendered for Overall_report.pdf — check source files exist.")
    return(FALSE)
  }
  .combine_pdfs(valid, output_path)
  unlink(valid[startsWith(valid, tempdir())])
  message("  Overall_report.pdf: ", length(valid), " pages -> ", output_path)
  TRUE
}

# Run A4-formatted Overall_report; fall back to simple PDF merge if pdftools absent
if (!make_overall_report(file.path(DIRS$reports, "Overall_report.pdf"))) {
  message("  Falling back to simple PDF merge for Overall_report.pdf")
  combine(
    c(
      file.path(DIRS$qc,       "qc_report.pdf"),
      file.path(DIRS$doublets, "doublets_report.pdf"),
      file.path(DIRS$annotation, "singler_scores_heatmap.pdf"),
      file.path(DIRS$annotation, "contamination_summary.pdf"),
      unlist(lapply(SAMPLE_NAMES, function(nm) c(
        file.path(DIRS$individual, nm, paste0(nm, "_hvg.pdf")),
        file.path(DIRS$individual, nm, paste0(nm, "_dotplot_markers.pdf"))
      ))),
      file.path(DIRS$integrated, "harmony_before_after.pdf"),
      file.path(DIRS$integrated, "integrated_umap_sample.pdf"),
      file.path(DIRS$integrated, "integrated_umap_celltype.pdf"),
      file.path(DIRS$integrated, "umap_split_by_sample.pdf"),
      file.path(DIRS$integrated, "umap_triptych.pdf"),
      file.path(DIRS$integrated, "integrated_dotplot.pdf"),
      file.path(DIRS$integrated, "integrated_heatmap.pdf"),
      file.path(DIRS$integrated, "celltype_composition_combined.pdf"),
      file.path(DIRS$integrated, "violin_key_markers.pdf")
    ),
    file.path(DIRS$reports, "Overall_report.pdf")
  )
}

message("\nFinal reports in: ", DIRS$reports)
message("  01-QC_report.pdf")
message("  02-Doublet_report.pdf")
message("  03-Individual_report.pdf")
message("  04-Annotation_report.pdf")
message("  05-Integrated_report.pdf")
message("  Overall_report.pdf")
message("\n07_finalize_reports.R complete.")
