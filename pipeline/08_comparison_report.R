# =============================================================================
# 08_comparison_report.R
# Targeted comparison report for multi-sample analysis.
# Pulls existing plots that address 5 follow-up questions:
#   1. Sample quality & cell fate (MT%, capture rates, QC thresholds)
#   2. Doublet rates (expected higher for high-input samples)
#   3. Cell type composition (monocyte depletion, lymphocyte enrichment)
#   4. Integrated UMAP (split by sample, triptych)
#   5. Key marker expression (violins, dotplot)
# Output: Comparison_report.pdf in RESULTS_DIR
# =============================================================================
.pipeline_dir <- local({
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(f)) dirname(f)
  else {
    a <- commandArgs(trailingOnly = FALSE)
    d <- sub("--file=", "", a[grep("--file=", a)])
    if (length(d) > 0) dirname(normalizePath(d)) else "."
  }
})
source(file.path(.pipeline_dir, "config.R"))   # config.R also sources pdf_helpers.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(dplyr)
})

if (!requireNamespace("pdftools", quietly = TRUE))
  stop("pdftools required. Install: install.packages('pdftools')")

message("\n=== Building Comparison Report ===")

# =============================================================================
# HELPERS (step-08-specific; A4P/A4L/.render/.build_page/.save_page from pdf_helpers.R)
# =============================================================================

pages <- character(0)
.add <- function(imgs, title, caption = NULL, landscape = FALSE, ncols = 1L) {
  spec <- .build_page(imgs, title, caption, landscape, ncols)
  tf   <- .save_page(spec)
  if (!is.null(tf)) pages <<- c(pages, tf)
}

# Text-only page (section header or interpretation)
.add_text <- function(header, body, landscape = FALSE) {
  dims <- if (landscape) A4L else A4P
  pg <- ggdraw() +
    draw_label(header, fontface = "bold", size = 15,
               x = 0.05, y = 0.93, hjust = 0, vjust = 1) +
    draw_label(body, size = 9, color = "#333333",
               x = 0.05, y = 0.85, hjust = 0, vjust = 1, lineheight = 1.4) +
    theme(plot.background = element_rect(fill = "#FAFAFA", color = NA))
  spec <- list(plot = pg, w = unname(dims["w"]), h = unname(dims["h"]))
  tf <- .save_page(spec)
  if (!is.null(tf)) pages <<- c(pages, tf)
}

# Table page from a data.frame
.add_table <- function(df, title, caption = NULL, landscape = FALSE) {
  dims <- if (landscape) A4L else A4P
  # Render table as ggplot
  tbl <- tableGrob(df, rows = NULL,
                   theme = ttheme_minimal(base_size = 8,
                     core    = list(fg_params = list(hjust = 0, x = 0.05)),
                     colhead = list(fg_params = list(fontface = "bold"))))
  rows <- list(); rel_h <- numeric(0)
  rows[[1]] <- ggdraw() +
    draw_label(title, fontface = "bold", size = 13,
               x = 0.03, y = 0.5, hjust = 0, vjust = 0.5)
  rel_h <- c(rel_h, 0.06)
  rows[[2]] <- ggdraw() + draw_grob(tbl)
  rel_h <- c(rel_h, if (!is.null(caption)) 0.78 else 0.94)
  if (!is.null(caption)) {
    rows[[3]] <- ggdraw() +
      draw_label(caption, size = 7.5, color = "#444444",
                 x = 0.03, y = 0.95, hjust = 0, vjust = 1, lineheight = 1.3)
    rel_h <- c(rel_h, 0.16)
  }
  pg   <- plot_grid(plotlist = rows, ncol = 1L, rel_heights = rel_h)
  spec <- list(plot = pg, w = unname(dims["w"]), h = unname(dims["h"]))
  tf   <- .save_page(spec)
  if (!is.null(tf)) pages <<- c(pages, tf)
}

# =============================================================================
# COVER PAGE
# =============================================================================

sample_list <- paste(SAMPLE_NAMES, collapse = ", ")
cover_body <- paste0(
  "Samples analysed: ", sample_list, "\n\n",
  "This report addresses five follow-up comparison questions:\n\n",
  "  1.  Sample Quality & Cell Fate\n",
  "      MT% distributions, QC thresholds, capture rates, and cell losses\n",
  "      per step for each sample.\n\n",
  "  2.  Doublet Rates\n",
  "      Doublet UMAP and score histogram per sample. High-input samples\n",
  "      (e.g. H2 loaded at 30,000 cells) are expected to have elevated\n",
  "      doublet rates (10-20% vs typical 5-8%).\n\n",
  "  3.  Cell Type Composition\n",
  "      Proportion bar charts comparing cell type abundance across samples.\n",
  "      Key question: do post-sorted samples (H2, H3) show monocyte\n",
  "      depletion and lymphocyte (CD4 T, B cell) enrichment vs DemoScRNA?\n\n",
  "  4.  Integrated UMAP\n",
  "      Triptych (cluster / sample / cell type) and per-sample split UMAPs.\n",
  "      Assesses Harmony batch correction quality and sample-specific\n",
  "      cluster contributions.\n\n",
  "  5.  Key Marker Expression\n",
  "      Violin plots and dot plot of lineage markers across samples.\n",
  "      Confirms cell type identities and detects sample-specific\n",
  "      expression differences.\n\n",
  "Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M")
)
.add_text("Comparison Report — Multi-Sample scRNA-seq Analysis", cover_body)

# =============================================================================
# SECTION 1: SAMPLE QUALITY & CELL FATE
# =============================================================================

.add_text(
  "Section 1: Sample Quality & Cell Fate",
  paste0(
    "What to look for:\n\n",
    "  Capture rate  = cells retained / live cells loaded. Expected 40-60% for 3' 10x Genomics.\n",
    "  H1 (pre-sort, 48.4% viability): expect very low capture and possible stress signature\n",
    "  (elevated MT%). Dead cells lyse before encapsulation, and debris reduces chip efficiency.\n\n",
    "  H2 (post-sort, old kit, 30,000 input): capture may exceed expected if cells are healthy.\n",
    "  H3 (post-sort, new kit, 30,000 input): compare recovery to H2 to evaluate kit performance.\n\n",
    "  MT% violin: H1 should show broader or higher MT% distribution than DemoScRNA.\n",
    "  High MT% (>15%) indicates mitochondrial leakage from damaged/dying cells.\n",
    "  nFeature (genes/cell) and nCount (UMIs/cell) reflect library complexity and depth."
  )
)

# Cell fate table
fate_file <- file.path(DIRS$qc, "cell_fate.csv")
if (file.exists(fate_file)) {
  suppressPackageStartupMessages(library(gridExtra))
  fate <- read.csv(fate_file, check.names = FALSE)
  # Simplify column names for readability
  colnames(fate) <- gsub("_", " ", colnames(fate))
  .add_table(
    fate,
    title   = "Cell Fate Summary (all samples)",
    caption = paste0(
      "GEM barcodes: total droplets sequenced by CellRanger.  After CreateSeuratObject: cells passing min.features >= 200.\n",
      "After QC: cells passing nFeature, nCount, and MT% thresholds.  After doublet removal: final singlets retained.\n",
      "Pct retained: final singlets / GEM barcodes.  Low % in H1 reflects poor viability (48.4%) — dead cells do not\n",
      "encapsulate efficiently, and debris occludes microfluidic channels reducing droplet formation."
    )
  )
}

# QC violins — 2 per page (portrait)
qc_pairs <- list(
  c("DemoScRNA", "H1_pre_old"),
  c("H2_post_old", "H3_post_new")
)
for (pair in qc_pairs) {
  nms   <- pair[file.exists(file.path(DIRS$qc, paste0(pair, "_violin_qc.pdf")))]
  imgs  <- lapply(nms, function(nm) .render(file.path(DIRS$qc, paste0(nm, "_violin_qc.pdf"))))
  .add(imgs,
       title   = paste0("QC Violin Plots — ", paste(nms, collapse = " & ")),
       caption = paste0(
         "Each violin shows the distribution of nFeature_RNA (genes per cell), nCount_RNA (UMIs per cell),\n",
         "and percent.mt (mitochondrial gene %) before QC filtering.  Dashed lines mark applied thresholds.\n",
         "Good: compact violins within threshold lines, MT% median < 5%.  Warning: H1 may show a wider or\n",
         "higher MT% tail due to 48.4% input viability — cells with MT% > 20% are damaged/dying.\n",
         "Compare nFeature distributions between pre-sort (H1, DemoScRNA) and post-sort (H2, H3) to assess\n",
         "whether sorting enriches for cells with higher transcriptional complexity."
       ),
       ncols = length(imgs))
}

# QC scatter plots (nCount vs percent.mt) — 2 per page
for (pair in qc_pairs) {
  nms  <- pair[file.exists(file.path(DIRS$qc, paste0(pair, "_scatter_qc.pdf")))]
  imgs <- lapply(nms, function(nm) .render(file.path(DIRS$qc, paste0(nm, "_scatter_qc.pdf"))))
  .add(imgs,
       title   = paste0("QC Scatter Plots — ", paste(nms, collapse = " & ")),
       caption = paste0(
         "Scatter of nCount_RNA (x) vs percent.mt (y) and nCount_RNA vs nFeature_RNA.\n",
         "Good: cells cluster in the lower-left (low MT%, moderate counts).  Outliers top-left = damaged cells\n",
         "(high MT%, low counts).  Outliers far-right = potential doublets (high counts/features).\n",
         "H1 may show more top-left outliers reflecting stressed/dying cells from poor viability."
       ),
       ncols = length(imgs))
}

# =============================================================================
# SECTION 2: DOUBLET RATES
# =============================================================================

.add_text(
  "Section 2: Doublet Rates",
  paste0(
    "What to look for:\n\n",
    "  scDblFinder scores each cell for doublet probability based on simulated\n",
    "  artificial doublets. Cells above the threshold are labelled 'doublet' and removed.\n\n",
    "  Expected doublet rates for 10x Genomics 3' Chromium:\n",
    "    ~1% per 1,000 cells loaded (manufacturer guideline)\n",
    "    Typical PBMC run (~5,000-8,000 cells): 5-8% doublets\n",
    "    H2 loaded 30,000 cells  =>  expected 15-25% doublet rate\n",
    "    H1 loaded 19,481 cells with 48.4% viability  =>  ~10-15% expected\n\n",
    "  High doublet rates in H2/H3 would manifest as clusters sitting between two\n",
    "  cell type groups on the UMAP (e.g. a T+Monocyte doublet cluster).\n",
    "  These are removed before integration, but high rates reduce final cell yield."
  )
)

for (pair in qc_pairs) {
  umap_imgs <- lapply(pair, function(nm)
    .render(file.path(DIRS$doublets, paste0(nm, "_doublet_umap.pdf"))))
  hist_imgs <- lapply(pair, function(nm)
    .render(file.path(DIRS$doublets, paste0(nm, "_doublet_score_hist.pdf"))))

  valid_umap <- Filter(Negate(is.null), umap_imgs)
  valid_hist <- Filter(Negate(is.null), hist_imgs)

  if (length(valid_umap) > 0)
    .add(valid_umap,
         title   = paste0("Doublet UMAP — ", paste(pair, collapse = " & ")),
         caption = paste0(
           "UMAP coloured by doublet classification (red = doublet, grey = singlet).  Title shows count and\n",
           "percentage removed.  Good: doublets distributed around cluster edges or in inter-cluster space.\n",
           "Warning: large central doublet cluster suggests cell clumping during loading."
         ),
         ncols = length(valid_umap))

  if (length(valid_hist) > 0)
    .add(valid_hist,
         title   = paste0("Doublet Score Histogram — ", paste(pair, collapse = " & ")),
         caption = paste0(
           "Distribution of scDblFinder doublet scores (0 = singlet confidence, 1 = doublet confidence).\n",
           "Good: bimodal distribution with a clear singlet peak near 0.  Threshold (dashed line) separates\n",
           "singlets from doublets.  Broad or unimodal distributions indicate ambiguity in classification."
         ),
         ncols = length(valid_hist))
}

# =============================================================================
# SECTION 3: CELL TYPE COMPOSITION
# =============================================================================

.add_text(
  "Section 3: Cell Type Composition",
  paste0(
    "What to look for:\n\n",
    "  DemoScRNA (standard healthy PBMC reference):\n",
    "    CD4 T ~40-50%,  CD8 T ~20-25%,  NK ~5-15%,  B cell ~5-10%,  Monocyte ~10-15%\n\n",
    "  H1 (pre-sort PBMC):  composition should resemble DemoScRNA.  If monocyte % is\n",
    "  lower than DemoScRNA, this suggests some pre-processing depletion occurred.\n\n",
    "  H2 & H3 (post-sort):  elevated CD4 T + B cells relative to DemoScRNA is the key\n",
    "  observation.  Most likely explanations:\n",
    "    (a) Monocyte/granulocyte depletion during sort — proportionally enriches lymphocytes\n",
    "    (b) Lymphocyte positive selection (CD3+CD19 gate)\n",
    "    (c) Viability-based gating — monocytes are more fragile and die faster\n\n",
    "  H2 vs H3 kit comparison:  if cell type proportions differ between H2 (old kit) and\n",
    "  H3 (new kit) for the same sorted population, the kits affect cell type representation."
  )
)

comp_combined  <- file.path(DIRS$integrated, "celltype_composition_combined.pdf")
comp_prop      <- file.path(DIRS$integrated, "celltype_proportions_bar.pdf")
comp_counts    <- file.path(DIRS$integrated, "celltype_counts_bar.pdf")

.add(list(.render(comp_combined)),
     title   = "Cell Type Composition — All Samples Combined",
     caption = paste0(
       "Stacked bar chart showing % of each cell type per sample.  Each bar = one sample; colours = cell types.\n",
       "Key comparison: monocyte % (salmon/brown) in DemoScRNA vs H2/H3.  A drop in monocytes in post-sort\n",
       "samples confirms depletion during PBMC sorting.  Corresponding rise in CD4 T (blue) and B cell (green)\n",
       "reflects the proportional enrichment of lymphocytes after monocyte/granulocyte removal.\n",
       "H1 (pre-sort): expect composition similar to DemoScRNA but with possible shift due to poor viability."
     ))

.add(list(.render(comp_prop), .render(comp_counts)),
     title   = "Cell Type Proportions & Absolute Counts",
     caption = paste0(
       "Left: proportion (%) per cell type per sample — normalises for cell number differences between samples.\n",
       "Right: absolute cell counts per cell type — reflects actual recovery.\n",
       "H2 (22,489 cells) vs H3 (9,267 cells) from the same 30,000-cell input: old kit recovers ~2.4x more cells.\n",
       "Compare absolute NK and monocyte counts across samples: if H2/H3 counts are near-zero for monocytes\n",
       "this confirms active depletion rather than proportional shift."
     ),
     ncols = 2L)

# =============================================================================
# SECTION 4: INTEGRATED UMAP
# =============================================================================

.add_text(
  "Section 4: Integrated UMAP (Harmony Batch Correction)",
  paste0(
    "What to look for:\n\n",
    "  Triptych (cluster / sample / cell type):  after Harmony integration, cells from all\n",
    "  four samples should be intermixed within each cluster.  Sample-specific clusters\n",
    "  indicate biological differences (expected for sorted vs unsorted) rather than batch\n",
    "  effects — these are real biology.\n\n",
    "  Split by sample:  each panel shows one sample's cells.  Compare which UMAP regions\n",
    "  are dense (enriched) or sparse (depleted) in each sample vs DemoScRNA reference.\n\n",
    "  H1 will appear sparse overall (only 2,211 cells after QC/doublet removal).\n",
    "  H2/H3 should be dense in T cell and B cell regions, sparse in monocyte region.\n\n",
    "  Harmony mixing check:  if a cluster is >80% from one sample, Harmony did not fully\n",
    "  correct batch — but this may also reflect genuine biology (e.g. a cell type absent\n",
    "  in one sample due to depletion)."
  )
)

umap_triptych <- file.path(DIRS$integrated, "umap_triptych.pdf")
umap_split    <- file.path(DIRS$integrated, "umap_split_by_sample.pdf")
umap_sample   <- file.path(DIRS$integrated, "integrated_umap_sample.pdf")
n_split_pages <- tryCatch(pdftools::pdf_length(umap_split), error = function(e) 0L)

.add(list(.render(umap_triptych)),
     title     = "Integrated UMAP Triptych (Cluster / Sample / Cell Type)",
     caption   = paste0(
       "Left: Seurat clusters (Louvain, res=0.5).  Centre: sample of origin — good mixing after Harmony\n",
       "means no sample dominates a cluster.  Right: cell type annotation (SingleR + sub-type refinement).\n",
       "Compare cluster boundaries to cell type colours: if a cluster is split across two cell types, consider\n",
       "increasing resolution.  H1 cells (sparse) may not form their own clusters — they slot into existing\n",
       "DemoScRNA/H2/H3 clusters."
     ),
     landscape = TRUE)

for (pg in seq_len(n_split_pages)) {
  .add(list(.render(umap_split, page = pg)),
       title   = sprintf("Integrated UMAP — Split by Sample (page %d/%d)", pg, n_split_pages),
       caption = paste0(
         "Each panel shows one sample's cells on the shared Harmony UMAP coordinate space.\n",
         "Dense regions = enriched cell types in that sample.  Empty regions = absent/depleted cell types.\n",
         "Expected pattern: H2/H3 panels dense in CD4 T and B cell region (upper UMAP), sparse in\n",
         "monocyte region (right UMAP).  DemoScRNA panel balanced across all regions.\n",
         "H1 panel sparse everywhere — low cell number (2,211 cells) limits cluster representation."
       ),
       landscape = TRUE)
}

.add(list(.render(umap_sample)),
     title     = "Integrated UMAP — Coloured by Sample",
     caption   = paste0(
       "All cells coloured by sample of origin on the Harmony-integrated UMAP.\n",
       "Good mixing: each cluster contains a mixture of colours (samples).  Poor mixing: monochromatic\n",
       "clusters indicate batch effects not corrected by Harmony — increase HARMONY$theta in config.R.\n",
       "Biologically distinct regions (e.g. monocyte cluster dominated by DemoScRNA/H1) are expected\n",
       "if H2/H3 are depleted of monocytes by sorting."
     ),
     landscape = TRUE)

# =============================================================================
# SECTION 5: KEY MARKER EXPRESSION
# =============================================================================

.add_text(
  "Section 5: Key Marker Expression",
  paste0(
    "What to look for:\n\n",
    "  Violin plots of canonical lineage markers split by sample confirm that cell type\n",
    "  identities are consistent across samples and that sorting did not artefactually\n",
    "  activate marker expression.\n\n",
    "  Dot plot (integrated):  spot size = % cells expressing the gene in that cluster;\n",
    "  colour intensity = average scaled expression.  Confirms cell type annotation.\n\n",
    "  Sorting-related artefacts to watch for:\n",
    "    - CD4/CD8 T cells in sorted samples may show elevated activation markers\n",
    "      (CD69, HLA-DR) if the sort procedure was slow or used stimulatory antibodies.\n",
    "    - Monocytes in sorted samples should be absent or near-zero in marker violins.\n",
    "    - B cells (MS4A1, CD79A) should show clear, high expression in H2/H3 if\n",
    "      enriched by sorting."
  )
)

violin_markers <- file.path(DIRS$integrated, "violin_key_markers.pdf")
n_vln_pages    <- tryCatch(pdftools::pdf_length(violin_markers), error = function(e) 0L)
for (pg in seq_len(n_vln_pages)) {
  .add(list(.render(violin_markers, page = pg)),
       title   = sprintf("Key Lineage Marker Violins (page %d/%d)", pg, n_vln_pages),
       caption = paste0(
         "Violin + jitter plots of canonical marker genes, split by sample (colour) and cell type cluster.\n",
         "Marker–cell type associations shown in subplot titles (e.g. 'CD3D (T cells)').\n",
         "Good: high expression in expected cell type, low in all others.  Compare expression levels\n",
         "across samples within the same cell type — sorting should not alter marker expression.\n",
         "Warning: if monocyte markers (CD14, LYZ) are elevated in H2/H3 monocyte-depleted samples,\n",
         "it may indicate incomplete depletion or contaminating monocyte debris."
       ))
}

dotplot <- file.path(DIRS$integrated, "integrated_dotplot.pdf")
.add(list(.render(dotplot)),
     title     = "Integrated Dot Plot — Canonical Markers by Cluster",
     caption   = paste0(
       "Dot size = fraction of cells expressing the gene (> 0) in that cluster.\n",
       "Dot colour = average scaled expression (blue = high, grey = low).\n",
       "Use this to cross-validate cell type annotations: each cell type should show high\n",
       "expression of its canonical markers and low expression of off-target markers.\n",
       "Annotation errors appear as marker expression in unexpected clusters."
     ),
     landscape = TRUE)

# Annotation canonical marker dotplot
anno_dotplot <- file.path(DIRS$annotation, "canonical_markers_dotplot.pdf")
.add(list(.render(anno_dotplot)),
     title     = "Annotation Dot Plot — SingleR + Sub-type Refinement",
     caption   = paste0(
       "Dot plot used during cell type annotation (step 05).  Shows canonical marker expression\n",
       "across all clusters used to validate SingleR predictions.\n",
       "Compare to the integrated dot plot above: cluster labels should match marker expression patterns.\n",
       "Discrepancies (e.g. a CD4 T cluster also expressing CD14) indicate mixed clusters or misannotation."
     ),
     landscape = TRUE)

# =============================================================================
# ASSEMBLE FINAL PDF
# =============================================================================

output_path <- file.path(RESULTS_DIR, "Comparison_report.pdf")
valid_pages <- pages[file.exists(pages)]

if (length(valid_pages) == 0) {
  message("  No pages generated — check that analysis has been run (steps 01-06).")
} else {
  .combine_pdfs(valid_pages, output_path)
  message("\n  -> ", output_path)
  message("     Pages: ", length(valid_pages))
}

message("\n08_comparison_report.R complete.")
