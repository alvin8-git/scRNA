# Design: Interactive HTML run report

**Date:** 2026-06-12 · **Mode:** office-hours design doc (no implementation) · **Reviewer:** Fable 5
**Mockup:** [`docs/mockups/scrna-report-template.html`](mockups/scrna-report-template.html) — open in a browser.

---

## Problem

Each run drops ~40 PDFs into `Results/results_*/` (qc/, doublets/, annotation/, integrated/, reports/). The `Overall_report.pdf` staples them into one paginated file. The actual job a wetlab person does with it — "did the sort deplete B cells? are the samples comparable? which sample is junk?" — means **flipping between pages to eyeball one metric across samples.** The data is per-sample-then-per-plot; the question is per-metric-across-samples. The format fights the question.

Concrete case from this project: confirming **B-cell depletion in pre- vs post-sorted cells.** Today that's: find the pre-sort proportions plot, find the post-sort one, hold them side by side in your head. The pipeline already computes everything needed to answer it in one view; it just never assembles it that way.

## What Partek does that's worth stealing

The datasheet's Figure 5 is a clean "data viewer": **three UMAPs in a row** (Protein / Integrated / Gene), shared layout, left tool rail, tabbed. Figure 4 is literally a **"PBMC not sorting" vs "PBMC cell sorting"** comparison. Their pattern is *faceting by modality*. Your version: **facet the same panel by sample**, with a metric switcher. Browser-based, one artifact, no install — a wetlab scientist opens it from email.

We are not rebuilding Partek (a platform). We are stealing one interaction: *one metric, every sample, side by side.*

## The wedge (smallest thing worth shipping)

A single self-contained `<run>_report.html` with two views that the PDF pile can't do:

1. **Sample comparison** — pick a metric (UMAP by cell type / QC / doublets), see all samples as columns in one row, shared axes.
2. **Cell proportions** — stacked composition per sample + a **pre→post delta table** that names the change (B cell 18% → 3%, −15 pts) and reads it for you.

Plus an **Overview** with per-sample summary-stat cards + one sortable table. Everything else (DE, pathways, per-cell drilldown) is a later tab, not the wedge.

The mockup renders all three. The UMAP/QC/doublet plots there are placeholder dot-clouds — the point is layout and interaction, not the render.

## A key call: share one coordinate space

The honest side-by-side UMAP is **the integrated object split by sample** (`umap_split_by_sample.pdf` already exists), not independent per-sample UMAPs. In the integrated space every sample shares one set of coordinates, so a depleted cluster shows up as a *hole in the same place* across panels. Independent per-sample UMAPs have different coordinates and can't be compared by eye. The interactive report should default to the split-integrated view, with per-sample UMAPs available as a drilldown.

## Data sources — every panel already exists

| Report panel | Source (already produced) |
|---|---|
| Overview stat cards / table | `qc/cell_fate.csv` (loaded, post-QC, % retained, doublet %), per-sample `_seurat.rds` metadata (median genes/UMI/%MT), `annotation/cluster_annotation_table.csv` (n types, top type) |
| QC comparison | per-sample `_filtered.rds` metadata (`nFeature_RNA`, `nCount_RNA`, `percent.mt`) → re-plot interactively |
| Doublet comparison | `_singlets.rds` metadata (`scDblFinder.class`/`score`) + rates from `cell_fate.csv` |
| UMAP atlas (hero) | `integrated_annotated.rds` UMAP embedding + `cell_type`, split by `$sample` |
| Cell proportions + delta | `cell_type` × `sample` counts from `integrated_annotated.rds`; CIs from `proportions/bootstrap_proportions.csv` |
| Differential expression | `differential/de_*.csv` (step 06b) |

No new computation. The report is an **assembly + presentation layer** over existing outputs.

## Architecture — how to generate the HTML (the real decision)

| | A. Quarto / R Markdown + plotly | B. Static template + exported SVG/JSON | C. Shiny / JS SPA |
|---|---|---|---|
| Engine | `quarto`/`rmarkdown` renders one `.html` from the `.rds` + CSVs | pipeline exports per-panel SVG + `summary.json`; a hand-built template (the mockup) fills it | served app reads a run bundle |
| Interactivity | hover, zoom, legend-toggle, linked brushing (plotly + crosstalk) | layout-level (toggle samples, swap metric); plots static unless re-authored | full cross-filter / click-a-cluster-filters-all |
| New tech in an R pipeline | low (Quarto is R-native) | medium (a JS template to maintain) | high (Shiny server or JS build) |
| Self-contained file (emailable) | yes (`embed-resources: true`) | yes | no (needs a running server) |
| Effort (human / CC) | ~2 days / ~½ day | ~3 days / ~1 day | ~1–2 wks / ~2 days |
| Risk | low | medium | high |

**Recommendation: A (Quarto + plotly), with the mockup's comparison-grid as the layout spec.** It stays in R, emits one self-contained `<run>_report.html` you can email, and gives real interactivity (hover a cell for its label/QC, toggle a cell type off in the legend and watch every panel update) for cheap. B is the fallback if Quarto/pandoc is unavailable in the conda env. C is the wrong shape for a static per-run artifact — it spends a big innovation token to serve a file a wetlab person should just open.

### Can Quarto do the dropdown + interactive plots? Yes.

**Interactive plots (not PDFs).** In Quarto the same ggplot chunk becomes interactive by wrapping it in `plotly::ggplotly()` (or building native plotly) — it embeds as an htmlwidget with hover, zoom, pan, and click-legend-to-toggle, all inside the self-contained `.html` (`embed-resources: true`). For big UMAPs (10 samples × thousands of cells) render with WebGL (`plotly::toWebGL()` / `scattergl`) so it stays smooth; downsample to ~3–5k cells/sample for the overview panels and keep full resolution on the per-sample drilldown.

**Selective sample dropdown (for >10 samples).** Three Quarto-native options, pick by how much power you need:
- **OJS `Inputs.select()`** (Observable JS, first-class in Quarto) — a multi-select that chooses which samples populate the comparison columns; the grid re-renders to only the selected ones. Best fit for "pick 2–4 of 12 to compare side by side".
- **plotly `updatemenus`** — a dropdown baked into a single plotly figure that shows/hides each sample's trace. No extra JS.
- **crosstalk `filter_select`** — links a shared dataset so one dropdown filters cells by sample across linked plots.

Caveat at scale: dropdown-to-select is light at any sample count; true linked-brushing across *all* large plotly widgets gets heavy past a handful. So the design is **select-then-compare**, not all-panels-always-live. The mockup now demonstrates the selector (checkbox chips; a default subset shown of N).

Ships as **one new step** (e.g. `08b_html_report.R` or fold into step 07) that reads the run's outputs and writes `reports/<run>_report.html`. Keep the PDFs — this is additive.

## Build plan (phased, after you pick an architecture)

1. **Phase 0 — schema.** Define `summary.json` per run (the stat table + proportions + paths). One function over `integrated_annotated.rds` + `cell_fate.csv`. This is the contract; both A and B consume it.
2. **Phase 1 — the wedge.** Overview + Sample comparison (UMAP/QC/doublets) + Proportions delta, as the mockup. plotly for the UMAP/violins; the proportions delta is a plain table.
3. **Phase 2 — drilldowns.** Per-sample tab, DE tab (06b), marker dot plots. Reuse existing CSVs.
4. **Phase 3 — polish.** Export-to-PDF button (print stylesheet), shareable section anchors, a "flag a sample" annotation a reviewer can type and re-export.

## Decisions (2026-06-12) + the one open call

- **Per-run only** (confirmed). Cross-run comparison is out of scope — to compare samples across protocols, include them in the same run. `summary.json` stays per-run; no cohort layer.
- **Overview: no proportions** (confirmed). Overview is stat cards + the sortable summary table only; proportions live in their own tab.
- **Many samples: select-then-compare** (confirmed). The comparison view renders only selected samples (dropdown/checkboxes), never all 10+ at once.
- **Architecture: A (Quarto + plotly)** unless `scrna_seurat` can't install Quarto/pandoc. Quarto covers both the interactive plots and the selector (see above). Confirm and I'll spec step `08b_html_report.R`.

- **Landing view: Overview** (confirmed). The report opens on Overview — per-sample stat cards + the sortable summary table (which flags a junk sample at a glance). The UMAP atlas and Cell proportions are tabs one click away, not the landing page.

All layout and scope decisions are now locked; the design is ready to build.

## The assignment

Design is locked (architecture A, Overview landing, select-then-compare, per-run). Next concrete step: confirm Quarto + pandoc install into the `scrna_seurat` conda env (`mamba install -n scrna_seurat -c conda-forge quarto pandoc r-plotly`), then spec `08b_html_report.R` — the function that reads a run's outputs (`integrated_annotated.rds`, `cell_fate.csv`, the proportion/DE CSVs) and renders `reports/<run>_report.html`. The mockup is the layout contract for that step.
