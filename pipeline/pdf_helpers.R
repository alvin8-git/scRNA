# pdf_helpers.R — shared PDF page-builder utilities for steps 07 and 08.
# Source this file after config.R has been sourced.

# A4 dimensions in inches
A4P <- c(w = 8.27,  h = 11.69)  # portrait
A4L <- c(w = 11.69, h = 8.27)   # landscape

# Combine a vector of PDF file paths into one PDF
.combine_pdfs <- function(files, output) {
  files <- files[file.exists(files)]
  if (length(files) == 0) { message("  No PDFs to combine for: ", output); return(invisible(NULL)) }
  if (length(files) == 1) { file.copy(files, output, overwrite = TRUE); return(invisible(output)) }
  if (requireNamespace("qpdf", quietly = TRUE)) {
    qpdf::pdf_combine(files, output)
  } else if (nchar(Sys.which("pdfunite")) > 0) {
    system2("pdfunite", c(shQuote(files), shQuote(output)))
  } else if (nchar(Sys.which("gs")) > 0) {
    system2("gs", c("-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=pdfwrite",
                    paste0("-sOutputFile=", shQuote(output)), shQuote(files)))
  } else {
    warning("Cannot combine PDFs — install R package 'qpdf' or system 'pdfunite'. Copying first file.")
    file.copy(files[1], output, overwrite = TRUE)
  }
  invisible(output)
}

# ------------------------------------------------------------------
# Render one page of a PDF to a nativeRaster (returns NULL on error)
# ------------------------------------------------------------------
.render <- function(path, page = 1L, dpi = 150L) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  n <- tryCatch(pdftools::pdf_length(path), error = function(e) 0L)
  if (n < page) return(NULL)
  suppressWarnings(tryCatch(
    pdftools::pdf_render_page(path, page = page, dpi = dpi, numeric = FALSE),
    error = function(e) { warning("PDF render failed (page ", page, " of ", path, "): ", conditionMessage(e)); NULL }
  ))
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
  suppressWarnings(tryCatch(
    { print(spec$plot); dev.off(); tf },
    error = function(e) { try(dev.off(), silent = TRUE); NULL }
  ))
}
