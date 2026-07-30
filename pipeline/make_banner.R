#!/usr/bin/env Rscript
# Build github_banner_collage.png from the 8-sample run: an annotated UMAP cluster
# atlas + the per-sample cell-type proportion bars, using the canonical palette.
suppressWarnings(suppressMessages({
  library(Seurat); library(ggplot2); library(patchwork)
}))
PIPE <- "/data/alvin/scRNA/pipeline"
# Source config.R BARE. It locates pdf_helpers.R via `sys.frame(1)$ofile`, so any
# wrapper call — suppressWarnings(), suppressMessages(), try() — adds a stack frame,
# ofile comes back NULL, the `%||% "."` fallback looks for ./pdf_helpers.R, and the
# source fails. Wrapped in try(silent=TRUE) that failure was invisible and the banner
# silently drew with the generic palette below instead of the canonical CELLTYPE_COLORS.
source(file.path(PIPE, "config.R"))
PAL_CANON <- if (exists("CELLTYPE_COLORS")) CELLTYPE_COLORS else NULL
if (is.null(PAL_CANON)) warning("CELLTYPE_COLORS unavailable - falling back to the generic palette")

R8  <- "/data/alvin/scRNA/Results/results_Aksh1-ES03-ES14-ES258-ES332-ES35-ES407-ES459_filtered"
rds <- file.path(R8, "integrated", "integrated_annotated.rds")
out <- "/data/alvin/scRNA/github_banner_collage.png"

cat("loading", rds, "...\n")
obj <- readRDS(rds)
md  <- obj@meta.data

col <- function(cands) { h <- intersect(cands, colnames(md)); if (length(h)) h[1] else NA }
ct_col <- col(c("cell_type","celltype","consensus_label","singler_label_clean","singler_label"))
sm_col <- col(c("sample","orig.ident","Sample"))
stopifnot(!is.na(ct_col), !is.na(sm_col))

reds <- Reductions(obj)
umap_name <- { u <- grep("umap", reds, value = TRUE, ignore.case = TRUE); if (length(u)) u[1] else reds[1] }
emb <- as.data.frame(Embeddings(obj, umap_name))[, 1:2]; colnames(emb) <- c("UMAP_1","UMAP_2")

ct <- as.character(md[[ct_col]]); ct[is.na(ct) | ct == ""] <- "Unassigned"
cells <- data.frame(cell_type = ct, sample = as.character(md[[sm_col]]),
                    UMAP_1 = emb[rownames(md),1], UMAP_2 = emb[rownames(md),2],
                    stringsAsFactors = FALSE)
rm(obj); gc()

cts  <- sort(unique(cells$cell_type))
base <- c("#4e79a7","#f28e2b","#e15759","#76b7b2","#59a14f","#edc948","#b07aa1",
          "#ff9da7","#9c755f","#bab0ac","#86bcb6","#d37295","#fabfd2","#8cd17d","#499894")
PAL <- if (!is.null(PAL_CANON)) {
  p <- unname(PAL_CANON[cts]); miss <- is.na(p); if (any(miss)) p[miss] <- rep(base, length.out = sum(miss)); setNames(p, cts)
} else setNames(rep(base, length.out = length(cts)), cts)

# downsample for the scatter (same 6000/sample cap as the report)
set.seed(1)
idx <- unlist(lapply(split(seq_len(nrow(cells)), cells$sample),
              function(ix) if (length(ix) > 6000) sample(ix, 6000) else ix), use.names = FALSE)
dd  <- cells[idx, ]
ctr <- aggregate(cbind(UMAP_1, UMAP_2) ~ cell_type, data = dd, FUN = median)

p_umap <- ggplot(dd, aes(UMAP_1, UMAP_2, color = cell_type)) +
  geom_point(size = 0.45, alpha = 0.85, stroke = 0) +
  ggrepel::geom_text_repel(data = ctr, aes(UMAP_1, UMAP_2, label = cell_type),
    inherit.aes = FALSE, size = 3.1, fontface = "bold", colour = "#16202e",
    bg.color = "white", bg.r = 0.14, max.overlaps = Inf, seed = 1,
    min.segment.length = 0, segment.size = 0.2, segment.colour = "grey60") +
  scale_color_manual(values = PAL, guide = "none") +
  labs(subtitle = "Integrated UMAP atlas · hover, zoom, toggle in the live report") +
  theme_void(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "#516074", size = 11, margin = margin(b = 6)),
        plot.margin = margin(8, 12, 8, 12))

tab  <- table(cells$cell_type, cells$sample)
prop <- as.data.frame(as.table(sweep(tab, 2, colSums(tab), "/") * 100), stringsAsFactors = FALSE)
colnames(prop) <- c("cell_type", "sample", "pct")
p_bar <- ggplot(prop, aes(sample, pct, fill = cell_type)) +
  geom_col(width = 0.82) +
  scale_fill_manual(values = PAL, guide = "none") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(subtitle = "Cell-type proportions per sample", y = NULL, x = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(colour = "#516074", size = 11, margin = margin(b = 6)),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 9, colour = "#516074"),
        axis.text.y = element_text(size = 8, colour = "#90a0b3"),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "#eef1f5"),
        plot.margin = margin(8, 14, 8, 6))

banner <- (p_umap | p_bar) + plot_layout(widths = c(1.35, 1)) +
  plot_annotation(
    title = "scRNA-seq Seurat pipeline — interactive HTML run report",
    subtitle = "Multi-sample QC → doublets → Harmony integration → annotation, compiled into one self-contained, emailable dashboard.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 19, colour = "#16202e"),
      plot.subtitle = element_text(size = 12, colour = "#516074", margin = margin(t = 2, b = 4)),
      plot.margin = margin(14, 16, 10, 16), plot.background = element_rect(fill = "white", colour = NA)))

dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
ggsave(out, banner, width = 12.6, height = 4.7, dpi = 132, bg = "white", device = dev)
cat("wrote", out, "(", round(file.size(out)/1024), "KB )\n")
