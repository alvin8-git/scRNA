# scRNA customer report — product design notes

Audience: customers with raw/filtered matrices but no bioinformatics staff and no
compute. The HTML report is deliverable #1; the Excel workbook and others follow.

Status: discussion document, nothing implemented. Written 2026-10-24.

## 0. Decisions taken

### Scope of "sample-agnostic"

"Agnostic" means these four validated sample types, and no more:

| sample type | species keyword | status |
|---|---|---|
| human PBMC | `human` | validated |
| human whole blood | `human` + `SINGLER_REF <- "MonacoImmune"` by hand | validated |
| cardiac differentiation (hESC/iPSC → CM) | `cm` | validated |
| bat whole blood (*Eonycteris spelaea*) | `bat` | validated |
| bat wing tissue | `bat_wing` | validated |

**Anything else — mouse included — is out of scope until analysed in detail.**
The report machinery will run on other species because the marker panel, palette
and grouping all come from config rather than being hardcoded, but running is
not the same as validated: there is no curated marker panel, no reference
suitable for `SINGLER_REF`, and no one has checked whether the cell-type calls
mean anything. Treat output from an unvalidated sample type as unreviewed.

- **Format: a single self-contained HTML.** No zip, no linked `figures/` folder.
  The base64 ~33% penalty and the one-off JS library cost are both accepted.
- **Core content, always present** — these four earn their bytes first:
  1. Sample overview
  2. Interactive UMAPs
  3. Cell proportions and numbers (the "what is present and how has it changed"
     section — the most informative part of the deliverable)
  4. The cell markers used to define each cell type
- **Everything else is optional**, admitted only if the size budget allows after
  the core is built, and otherwise shipped as **supplementary HTML reports**:
  QC violins and scatters, doublet detection, HVG, knee/barcode-rank, saturation
  and rarefaction, ambient rho, integration quality, cell-cycle.

### What that scoping is worth, in measured bytes

Removing the non-core sections from the current 13.2 MB report:

| removed | size |
|---|---|
| Per-sample galleries (QC scatters ×2, HVG, doublet UMAP, doublet histogram) | 3.26 MB |
| QC violin widget | 1.87 MB |
| Doublet widget | 0.70 MB |
| **subtotal removed** | **5.83 MB** |

Core set, as currently built with plotly:

| kept | size |
|---|---|
| Head: plotly.js, CSS, jQuery, DataTables, fonts | 4.59 MB |
| 7 UMAP `scattergl` widgets | 2.34 MB |
| Stage-gallery PNGs (SingleR heatmap, lineage markers) | 0.32 MB |
| Composition, cell numbers, dotplot, DataTables | <0.10 MB |
| **core total** | **≈ 7.35 MB** |

So scoping alone lands the core under 10 MB — but **62% of what remains is the
plotly.js bundle plus UMAP widget JSON**, carrying only one of the four core
sections.

### BENCHMARKED 2026-10-24 — ggiraph is NOT lighter. Measured, self-contained:

| component | library | size |
|---|---|---|
| library floor (empty widget) | plotly | **3.78 MB** |
| library floor (empty widget) | ggiraph | **8.74 MB** |
| proportions stacked bar (77 marks) | plotly | 3.78 MB (marginal 0.01) |
| proportions stacked bar (77 marks) | ggiraph | 8.75 MB (marginal 0.01) |
| UMAP scattergl, 42,000 of 76,314 cells | plotly | 6.35 MB (marginal **2.58**) |
| UMAP raster, **all 76,314 cells**, png | scattermore | 0.12 MB |
| └ same, base64-embedded (×1.37) | scattermore | **0.17 MB** |
| UMAP annotation overlay (11 marks) | ggiraph | 8.73 MB (marginal 0.01) |

**The ggiraph recommendation was wrong.** Its floor is 8.74 MB — 2.3× plotly's —
because `gdtools`/`gfonts` bundle web fonts into every widget. And its marginal
cost per chart is 0.01 MB, identical to plotly's. It offers no advantage at
either end. Drop it; plotly is already installed and is the cheaper library here.

**The library was never the main lever — the UMAP strategy is.** Rasterising the
point cloud costs 0.17 MB embedded and plots *all* 76,314 cells; a plotly
`scattergl` widget costs 2.58 MB per 42,000 points, a **15× difference**, and
plots only 55% of them.

### Corrected core projection

My first projection (21.93 MB all-plotly) was wrong: it multiplied the cost of
one combined 42,000-point widget by seven. The report builds *per-sample*
widgets of ≤6,000 points each, measured at ~0.34 MB apiece — which agrees with
the single-widget measurement (2.58 ÷ 7 = 0.37). Correct figures:

| option | core total | cells plotted |
|---|---|---|
| plotly + 7 per-sample `scattergl` (today) | ≈ 6.4 MB | 42,000 of 76,314 |
| **plotly + 7 raster UMAPs** | **≈ 5.1 MB** | **all 76,314** |
| ggiraph + raster | ≈ 10.0 MB | all 76,314 |

**Recommendation: keep plotly, rasterise the UMAPs.** The core lands near 5 MB —
of which 3.78 MB is the unavoidable plotly floor — leaving roughly 4 MB of
headroom under a 9 MB ceiling for the optional sections, each of which costs
tens of KB.

If the 3.78 MB floor itself becomes the binding constraint, the next move is to
drop widgets entirely for the low-cardinality charts and emit plain inline SVG
with CSS tooltips — those charts only need 0.01 MB of marks, so the whole
JavaScript library is pure overhead for them. That is worth measuring only if
5 MB proves too large.

Recommended admission order for that headroom, by value per byte: sample QC
summary table (KBs), ambient rho per sample (KBs), proportion confidence
intervals (KBs), pseudobulk sample-correlation heatmap (KBs), knee plot,
saturation/rarefaction, then the per-sample QC and doublet galleries last since
they are the most expensive and the least informative to this audience.

## 1. What is actually there today (measured on H1D01sx_7samples_report.html)

**The report is already substantially interactive.** This corrects the premise
that it is mostly imported screenshots — it is not. Measured breakdown of the
13,232,379-byte file:

| component | size | share |
|---|---|---|
| Head: plotly.js (3.57 MB), CSS (0.64 MB), jQuery, DataTables, fonts | 4.59 MB | 35% |
| Widget JSON (plotly data payloads) | ~4.95 MB | 37% |
| base64 `data:` URIs — 47 of them, 38 PNG (2.67 MB decoded) | 4.00 MB | 30% |
| └ of which per-sample galleries | 3.26 MB | 25% |
| └ of which stage-gallery PNGs | 0.32 MB | 2% |

Already plotly widgets, no PDF involved: 7 UMAP `scattergl` (~334 KB each),
QC violin (**1.87 MB — the single largest widget**), doublets (0.70 MB),
composition / cell numbers / dotplot / DataTables (<100 KB combined).

**Only three PDFs are harvested**, not "many":

| figure | source | cost across 7 samples |
|---|---|---|
| Doublet UMAP | rasterised `<s>_doublet_umap.pdf` p1 @120 dpi (step 02) | 828 KB |
| Doublet score histogram | rasterised `<s>_doublet_score_hist.pdf` p1 @140 dpi (step 02) | 221 KB |
| SingleR score heatmap | rasterised `annotation/singler_scores_heatmap.pdf` p1 @130 dpi (step 05) | 188 KB |

Everything else in the per-sample galleries (QC UMIs vs genes 575 KB, UMIs vs
%MT 382 KB, HVG 429 KB) is **already rendered in-process from `ggsave`** — it is
a static PNG by choice, not because a PDF was in the way.

Nothing from step 07 is read at all.

**No plot objects are saved by steps 01–06b** — the only `saveRDS` calls write
Seurat objects. The doublet UMAP, doublet histogram and SingleR heatmap cannot
be rebuilt from those, which is *why* 08b rasterises them. That is the real
coupling to fix, and it affects exactly three figure types.

UMAPs are capped at 6,000 cells per sample, so the report plots 42,000 of the
cohort's 76,314 cells.

## 2. Size budget — and the tension in the request

Target: **≤ 10 MB self-contained**, ideally ≤ 8 MB. Current: 13.2 MB.

**The uncomfortable finding: fixing the "screenshots" will not get under 10 MB.**
All 38 PNGs together are 2.67 MB decoded. Meanwhile plotly.js plus widget JSON
is 9.54 MB — **72% of the file**. Adding more interactivity the current way makes
the problem worse, not better.

The size levers, largest first:

| lever | saving | cost |
|---|---|---|
| Drop plotly for `ggiraph` (SVG + d3, ~200 KB lib) | **−3.3 MB** fixed | loses WebGL; dense scatter must become raster |
| Replace the 1.87 MB QC violin widget with a ggiraph box/violin over summary stats | **−1.8 MB** | per-cell points no longer hoverable (no loss for this audience) |
| Raster UMAP (`scattermore`) + SVG annotation overlay, instead of 7 scattergl widgets | **−2.3 MB** JSON, +~0.5 MB raster | per-cell hover replaced by per-cluster hover |
| One UMAP widget with a sample selector instead of 7 | −1.7 MB | if staying on plotly |
| WebP instead of PNG for galleries | −1.0 to −1.6 MB | none worth noting |

Taking the first three: **13.2 MB → roughly 5 MB, with *more* interactivity on
the charts customers actually hover.** That is the plan worth pursuing.

So the rule: **interactivity is cheap for low-cardinality marks and expensive for
dense point clouds.** Spend it accordingly.

## 3. Proposed architecture — tier the charts by cardinality

### Tier A — dense point clouds: raster base + interactive overlay
UMAP (all variants), QC scatter (genes vs UMI), doublet-score UMAP.

`scattermore` is already installed and exists precisely for this: it rasterises
millions of points to a bitmap fast. Render the point cloud as a PNG, then
overlay *only* the annotations as SVG/HTML — cluster centroid labels, a legend
with click-to-toggle, and hover zones per cluster rather than per cell.

The customer gets: hover a cluster → cell type, n cells, % of sample, top 5
markers. They do not get per-cell hover, which is of no value to a non-specialist
anyway and costs 2–3 MB.

Optional upgrade: `plotly` `scattergl` (WebGL) on a **downsampled** 10–15k cell
layer, with the full-resolution raster underneath. Honest caveat to put in the
report: the interactive layer is a sample, the picture is complete.

### Tier B — low-cardinality summaries: fully interactive, ~20–80 KB each
Cell-type proportions (stacked + grouped), per-sample QC summary bars, HVG
top-N, marker dotplots, doublet rate per sample, SoupX rho per sample,
cell-cycle phase, rarefaction curves, barcode-rank/knee.

These have tens to low-hundreds of marks. Interactive tooltips here are nearly
free and are where customers actually hover.

**Library choice matters more than it looks.** `ggiraph` (SVG + d3) is ~200 KB
and renders ggplot objects directly with `tooltip=`/`data_id=` aesthetics — about
20× lighter than plotly and a near-drop-in for a ggplot-based pipeline. It is
**not currently installed**. If Tier A uses rasters rather than WebGL, ggiraph
alone covers everything and plotly can be dropped entirely, saving ~3.3 MB.
Recommendation: benchmark both on one real report before committing.

### Tier C — tables: paged, or pushed to the Excel deliverable
`DT` is installed. A paged DataTable of the top 20 markers per cluster is small;
the full marker table (thousands of rows × 7 columns) belongs in Excel, not in
the HTML. Link, don't embed.

## 4. Pipeline change — stop harvesting PDFs

The structural fix, and the reason the pipeline feels "chunky":

Narrower than it first appeared: **three figure types**, not the whole report.

**Today:** step 02 writes `<s>_doublet_umap.pdf` and `<s>_doublet_score_hist.pdf`
and step 05 writes `annotation/singler_scores_heatmap.pdf`; 08b rasterises page 1
of each. It must, because no step saves the plot objects and they cannot be
rebuilt from the Seurat `.rds` alone.

**Proposed:** every plot-producing step saves its figures as plot **objects**
(`.rds` of the ggplot, or a small list with the data + a builder function) into
`figures/<step>/<name>.rds`, alongside continuing to write the PDFs for anyone
who wants print output. Step 08b then consumes objects and renders each figure
*for the web* at web dimensions, interactive where Tier B applies.

Benefits:
- The PDF and the HTML stop fighting over one layout. Print wants 2-up A4;
  the web wants one chart, full width, responsive.
- No `magick`/poppler dependency in the HTML path, and no dpi guessing.
- Figures become individually addressable, so the report gains a manifest and
  missing figures degrade to a named placeholder rather than a silent gap.
  (This is the class of bug that silently dropped 18–24 pages per report before
  commit `a7e8020`.)
- A figure can be regenerated without re-running the analysis step.

Migration can be incremental: add object-saving to steps 03/05/06 first, switch
08b section by section, and keep the PDF-rasterising fallback until every
section has an object source. That keeps the report working throughout.

## 5. Charts worth adding (customer value, low size)

Ranked by value-per-byte for a non-specialist audience:

1. **Barcode-rank / knee plot** per sample — the standard "did cell calling
   work" figure. Customers coming from CellRanger/DNBelab reports expect it and
   notice its absence.
2. **Per-sample QC summary table** — cells, median genes, median UMI, % mito,
   doublet %, saturation, in one scannable grid. Highest information density in
   the whole report.
3. **Sequencing saturation + rarefaction** — answers "did you sequence deep
   enough, should I have paid for more reads". `pipeline/10_rarefaction.R`
   already computes this and it is not in the report.
4. **Ambient RNA (SoupX rho) per sample** — a genuine differentiator. Most
   vendors do not report contamination at all. This cohort ranged rho
   0.053–0.795, and the 0.795 library was unusable; a customer would want to
   know that before drawing conclusions.
5. **Pseudobulk sample–sample correlation heatmap** — catches swapped or
   mislabelled samples and confirms replicate concordance. Cheap (n×n marks),
   and it is the check that settled the H1D11_1 identity question.
6. **Proportion confidence intervals** — `pipeline/09_bootstrap_proportions.R`
   exists and is not in the report. Without CIs a customer will over-read a
   2-point difference between samples.
7. **Integration quality** — UMAP coloured by sample before/after Harmony, plus
   a per-cluster sample-entropy bar. Shows batch effect was handled rather than
   asserting it.
8. **Cell-cycle phase per cluster** — explains away "why is this cluster
   separate" questions before they are asked.
9. **Top-N marker heatmap per cell type** — the evidence behind the labels.
10. **Cell-cycle / mito vs complexity scatter** — already present as QC scatter.

Deliberately **not** recommended for the HTML: per-gene interactive lookup
across all ~30k genes. It is the single largest size driver and belongs in a
separate deliverable (below).

### Plain-language annotation is the actual product
For this audience, every section needs one or two sentences of "what this shows
and what good looks like", with the threshold stated. That is worth more than any
additional chart, and costs kilobytes. Example: "Doublets are two cells captured
in one droplet. We removed N (X%). Above ~10% suggests the sample was loaded too
densely."

## 6. Other deliverables worth offering

| deliverable | effort | customer value |
|---|---|---|
| **Excel: markers + proportions + QC** | low — CSVs exist | high. Already partly built: `cardiac_expression.xlsx` is gene × sample pseudobulk. The natural companion is gene × **cell type**. |
| **Methods & parameters appendix** | low | required for publication; QC thresholds, software versions, reference used |
| **Seurat `.rds` / `.h5ad` export** | low | lets the customer hire a bioinformatician later without paying for reprocessing. Strong goodwill, near-zero cost |
| **Per-sample 1-page PDF** | low | for printing and for pasting into slide decks |
| **Interactive cell browser** (cellxgene / Vitessce bundle) | medium–high | the per-gene lookup customers eventually ask for, delivered as its own artefact rather than bloating the report |
| **Provenance manifest** | low | input checksums, run parameters, pipeline commit — makes the report defensible |

## 7. Notes on the four core sections

**Sample overview.** The per-sample QC summary grid is the highest information
density in the report and costs kilobytes: cells, median genes, median UMI,
% mito, doublet %, saturation, and — worth including since most vendors omit it
— ambient rho. One scannable table, no figure.

**Interactive UMAPs.** Currently 7 `scattergl` widgets capped at 6,000 cells per
sample, so 42,000 of 76,314 cells are plotted. Switching to a `scattermore`
raster plots *all* cells at lower byte cost, with hover moving from per-cell to
per-cluster (cell type, n, % of sample, top markers). For this audience that is
a better interaction as well as a cheaper one. Needed views: coloured by cell
type, by sample, and by timepoint/condition.

**Proportions and numbers.** The user's framing — "what is present and how it has
changed" — means a change view is required, not just composition. Three charts:
stacked proportion per sample, absolute cell counts per type, and a **delta view**
across condition or timepoint. Add bootstrap confidence intervals
(`09_bootstrap_proportions.R` already computes them): without CIs a customer will
over-read a 2-point difference. All low-cardinality, so fully interactive for
tens of KB.

**Cell markers used to define cell types.** This is the evidence behind the
labels, and the curated maps in `config_species_cm.R` already record it as
per-cluster detection percentages. Surface exactly that: a cell type × marker
detection-% matrix as an interactive heatmap or table (22 clusters × ~45 markers
is under 1,000 values — trivial in size), alongside the existing dotplot. State
the rule the curation used, since it is the defensible part: identity
transcription factors decide a call, secreted cargo never does.

## 7b. Value-per-byte analysis for further additions

Measured baseline: the customer report is **10.47 MB** against a 20 MB ceiling, so
**≈9.5 MB of headroom**. Costs below are measured where a comparable figure
already exists in a built report, estimated only where marked.

| candidate | cost | value to a non-specialist | verdict |
|---|---|---|---|
| **Proportion confidence intervals** | 0.01 MB | High — without them a customer reads a 2-point difference as real. Already computed by `09_bootstrap_proportions.R` | **add first** |
| **Saturation / rarefaction** | ~0.05 MB | High — answers "should I have paid for more reads", the question customers actually ask. `10_rarefaction.R` already computes it | **add** |
| **Knee / barcode-rank** | ~0.05 MB | High — the standard "did cell calling work" figure; its absence is noticed by anyone who has seen a CellRanger report. Needs the raw matrices | **add** |
| **Integration before/after UMAP** | ~0.3 MB (2 rasters) | High — demonstrates the batch effect was handled instead of asserting it | **add** |
| **HVG per sample** | 0.43 MB (measured, 7 samples) | Low for this audience — variable-gene selection is a method detail, not a result. Keep in the internal report | supplementary |
| **Doublet UMAP + score histogram per sample** | 1.05 MB (measured) | Moderate — builds QC credibility, but the run-level Doublets tab already covers the finding. **Blocked**: needs the step-02 per-cell fix or it stays a rasterised PDF | supplementary, after the step-02 fix |
| **Cell-cycle phase per cluster** | ~0.05 MB | Moderate — pre-empts "why is this cluster separate". Needs scoring to be run | optional |
| **Top-N marker heatmap** | ~0.1 MB | Low — the dot plot and the new per-cluster evidence table already carry this | skip |
| **Per-gene expression lookup** | >3 MB | High demand, wrong vehicle — belongs in a separate browser deliverable | separate deliverable |

**The four "add" rows total ~0.4 MB**, taking the report to roughly 10.9 MB. Every
one answers a question a customer asks unprompted. The two "supplementary" rows
cost 1.5 MB between them and answer questions only a bioinformatician asks —
which is the right dividing line for this product.

**The size lesson holds**: charts are nearly free, libraries are not. plotly.js
plus DataTables and jQuery is ~4.6 MB, roughly 44% of the current report, and no
amount of chart curation changes that. Adding every row in the table above
would still leave the report under 13 MB.

## 8. Open questions

1. **Interactivity library:** install `ggiraph` and drop plotly (−3.3 MB, SVG
   only, dense scatter must become raster), or stay on plotly for WebGL? This is
   now the only decision blocking implementation. Recommendation: ggiraph plus
   `scattermore` rasters, benchmarked on one real report before committing.
2. **Scope:** cohort-agnostic product template (bat, human PBMC, cardiac, wing)
   or cardiac-specific first? The species overlays already differ in markers and
   QC thresholds, so a template needs the marker panel driven entirely from
   config rather than hardcoded.
3. **Per-gene lookup:** confirmed out of the main report; offer as a separate
   browser deliverable if customers ask.
