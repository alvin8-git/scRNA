# Pipeline Documentation

## Overview

This pipeline takes raw single-cell RNA-seq (scRNA-seq) count matrices and produces annotated cell populations, per-sample composition statistics, and publication-quality figures, using Seurat v5 in R. It handles **human PBMC** and **bat (*Eonycteris spelaea*) whole blood**, from both **10x Genomics** and **DNB C4** (dnbc4tools) platforms.

### New to single-cell analysis? The arc in one paragraph

scRNA-seq measures gene expression in thousands of individual cells at once: each cell is a column, each gene a row, each value a UMI (unique transcript) count. The raw matrix is noisy — empty droplets, dying cells, and **doublets** (two cells captured as one droplet) all masquerade as real cells. The pipeline removes that noise, then finds and names the biology:

1. **QC** (step 01) — drop low-quality cells: too few genes detected, or too many mitochondrial reads (a dying-cell signature).
2. **Doublet detection** (step 02) — flag and remove droplets that captured two cells.
3. **Per-sample processing** (step 03) — normalise counts, pick the most informative genes (highly variable genes, HVGs), reduce dimensions (PCA → UMAP), and **cluster** cells by expression similarity.
4. **Integration** (step 04) — with more than one sample, **Harmony** removes batch effects so the same cell type from different samples overlaps instead of splitting by sample.
5. **Annotation** (step 05) — give each cluster a biological identity (CD4 T, NK, Monocyte, …) using a reference atlas (SingleR), marker genes (scType), a consensus of the two, and an optional manual map.
6. **Quantify & report** (steps 06–10) — figures, differential expression, bootstrap proportion confidence intervals, and rarefaction (how many cells you need for a stable estimate).

Every step writes an `.rds` checkpoint, so you can re-run from any point. A pre-flight validator and a per-sample cache (see Configuration Reference) keep re-runs fast and safe.

---

## Data Format Notes (DNB C4 vs 10X)

The DNB C4 output is compatible with Seurat's `Read10X()` but has one key difference:

| | 10X Genomics | DNB C4 |
|--|--|--|
| `features.tsv.gz` columns | 2 (ENSEMBL ID, gene symbol) | 3 (ENSEMBL ID, gene symbol, feature type) |
| Barcode format | `ACGTACGT-1` | `CELL1_N3` |
| Matrix format | MTX (gzip) | MTX (gzip) — identical |

**Critical:** always pass `gene.column = 2` to `Read10X()` to use gene symbols instead of ENSEMBL IDs with version numbers (e.g. `ENSG00000243485.5`).

---

## Pipeline Steps

The pipeline is a sequence of R scripts in `pipeline/`. Core steps (run by default):

| Step | Script | What it does |
|------|--------|--------------|
| 00 | `validate_config.R` | Pre-flight config validation (see Configuration Reference) |
| 01 | `01_load_qc.R` | Load matrices, compute QC metrics, filter cells |
| 02 | `02_doublets.R` | scDblFinder doublet detection + removal |
| 03 | `03_individual.R` | Per-sample normalise → HVG → PCA → UMAP → cluster → markers |
| 04 | `04_integrate.R` | Merge + Harmony batch integration (multi-sample) |
| 05 | `05_annotate.R` | SingleR + scType + consensus + manual cell-type annotation |
| 05r | `05r_reference_transfer.R` | Run-independent `cell_type_ref` labels from a frozen reference (optional; gated on `REFERENCE_MODEL`) |
| 06 | `06_visualize.R` | Publication figures |
| 06b | `06b_differential.R` | Differential expression between samples per cell type (multi-sample) |
| 07 | `07_finalize_reports.R` | Merge per-step PDFs into 5 category reports + `Overall_report.pdf` |
| 08 | `08_comparison_report.R` | Standalone cross-sample comparison report |
| 08c | `08c_benchmark_concordance.R` | Cross-run anchor benchmark vs frozen baseline + whole-blood sort readout (optional; gated on `REFERENCE_MODEL`) |
| 09 | `09_bootstrap_proportions.R` | Bootstrap-normalised proportion CIs + pairwise chi-squared |
| 10 | `10_rarefaction.R` | Minimum-capture-depth rarefaction analysis |

Bat-wing-specific downstream analysis (steps 11–14: wing DEGs, pathway enrichment, CellChat, trajectory) lives under `pipeline/projects/bat_wing/` and runs only in `bat_wing` species mode.

The sections below detail each core step.

### Step 01 — Load & QC (`01_load_qc.R`)

**Input:** `H1/filter_matrix/`, `H2/filter_matrix/`  
**Output:** `results/qc/*.pdf`, `results/individual/H{1,2}/H{1,2}_filtered.rds`

Loads each sample with `Read10X(gene.column = 2)`, adds sample metadata, computes:
- `nFeature_RNA` — number of detected genes per cell
- `nCount_RNA` — total UMI count per cell
- `percent.mt` — percentage of reads mapping to mitochondrial genes (`^MT-`)
- `percent.ribo` — percentage mapping to ribosomal genes (`^RP[SL]`) (diagnostic only)

Default QC thresholds (adjust in `config.R`):

| Metric | Min | Max | Rationale |
|--------|-----|-----|-----------|
| nFeature_RNA | 200 | 5000 | Removes empty droplets and potential doublets |
| nCount_RNA | 500 | 25000 | Removes debris and outlier high-UMI cells |
| percent.mt | — | 20% | PBMC can legitimately reach 15-20% |

**Inspect** `H1_violin_qc.pdf` and `H2_violin_qc.pdf` before accepting defaults. For H1 (only 386 input cells), relax `max_features` to 6000 if a genuine population falls above 5000.

Barcodes are prefixed with the sample name (`RenameCells(add.cell.id = sample_name)`) so `H1_CELL1_N3` and `H2_CELL1_N3` remain distinct after merging.

---

### Step 02 — Doublet Detection (`02_doublets.R`)

**Input:** `H{1,2}_filtered.rds`  
**Output:** `results/doublets/*.pdf`, `H{1,2}_singlets.rds`

Uses **scDblFinder** (Bioconductor), which is natively compatible with Seurat v5. The tool:
1. Converts the Seurat object to `SingleCellExperiment`
2. Generates artificial doublets by random cell fusion
3. Trains a gradient boosting classifier to score each cell
4. Labels cells as `singlet` or `doublet`

Expected doublet rates used:
- H1 (386 cells): **3.1%** → ~12 expected doublets
- H2 (967 cells): **7.7%** → ~75 expected doublets

Adjust rates in `config.R → DOUBLET$doublet_rate`.

**Outputs per sample:**
- `*_doublet_umap.pdf` — UMAP with doublets highlighted in red
- `*_doublet_score_hist.pdf` — score distribution (bimodal = clean separation)

---

### Step 03 — Individual Analysis (`03_individual.R`)

**Input:** `H{1,2}_singlets.rds`  
**Output:** `results/individual/H{1,2}/*.pdf`, `H{1,2}_seurat.rds`

Full per-sample Seurat workflow:

1. **Normalisation** — `LogNormalize` with scale factor 10,000
2. **HVG selection** — `FindVariableFeatures(nfeatures = 2000, method = "vst")`
3. **Scaling** — `ScaleData(vars.to.regress = "percent.mt")` removes mitochondrial confounding
4. **PCA** — 30 PCs computed; elbow plot + PC heatmaps saved
5. **UMAP** — `dims = 1:20`, `seed = 42`
6. **Clustering** — `FindNeighbors` + `FindClusters` at resolutions 0.3, 0.4, 0.5, 0.6, 0.8; default `res = 0.5`
7. **Marker genes** — `FindAllMarkers(only.pos = TRUE, test = "wilcox", min.pct = 0.1, logfc.threshold = 0.25)`

For H1 (~350 cells), consider using `resolution = 0.4` to avoid over-clustering.

---

### Step 04 — Harmony Integration (`04_integrate.R`)

**Input:** `H{1,2}_seurat.rds`  
**Output:** `results/integrated/*.pdf`, `integrated_seurat.rds`

1. **Merge** — `merge(seu_H1, seu_H2)` then `JoinLayers()` (required in Seurat v5 after merge)
2. **Re-normalise** on merged object (normalisation must be done after merge)
3. **PCA** on merged HVGs
4. **Harmony** — runs on PCA embedding, returns a `harmony` reduction
   - `theta = 2`: diversity penalty; increase to 3-4 if samples don't intermix
   - `lambda = 1`: ridge regression penalty
   - `nclust = 50`: k-means clusters (~nCells/25 rule of thumb)
5. **UMAP** on Harmony embedding
6. **Clustering** on Harmony-corrected neighbours

**Key diagnostic:** `harmony_before_after.pdf` — H1/H2 should intermix after correction while cell type clusters remain separated. If cell types blur together, reduce `theta`.

---

### Step 05 — Annotation (`05_annotate.R`)

**Input:** `integrated_seurat.rds`  
**Output:** `results/annotation/*.pdf`, `cluster_annotation_table.csv`, `integrated_annotated.rds`

**Part 1: SingleR automated annotation**
- Reference controlled by `SINGLER_REF` in `config.R`:
  - `"HumanPrimaryCellAtlas"` (default) — broad PBMC cell types (~30 types)
  - `"MonacoImmune"` — blood-optimised, resolves CD4/CD8/Treg/γδ T, monocyte subtypes, pDC/mDC (29 types); set automatically by `bat` species override
- Run at single-cell level with `fine.tune = TRUE`
- Low-confidence calls pruned (`prune = TRUE`)
- Delta score (confidence margin) saved as `singler_delta`
- Raw labels normalised via `SINGLER_NORM` lookup — covers HumanPrimaryCell granulocyte/erythroid lineages + full Monaco label.main set

**Part 1b: scType marker-based scoring**
- Runs alongside SingleR with no reference transcriptome — uses `MARKERS` gene lists directly
- Genes weighted by inverse frequency (markers found in fewer cell types score higher, reducing false positives from shared markers e.g. FCGR3B shared between Neutrophil and FCGR3A_mono)
- Filters genes against `rownames(scale.data)` — HVG-safe (scale.data only contains the 2000 HVGs)
- Outputs: `sctype_labels_umap.pdf`, `singler_vs_sctype_comparison.csv`
- `sctype_label` added to metadata for cross-reference with SingleR

**Part 1c: Consensus auto-annotation**
- Fuses the two per-cluster calls (SingleR majority + scType) into a per-cluster **AUTO/REVIEW** decision: clusters where both independent methods agree are auto-accepted; where they disagree, flagged for review.
- Writes `consensus_annotation.csv` (per cluster: SingleR label, scType label, agreement, mean delta, decision) and prints a **pre-filled `CLUSTER_CELLTYPE_MAP`** to the log with `REVIEW` markers on the ambiguous clusters.
- Advisory only — it does not change the assignment logic below. It collapses "hand-build the whole map" into "resolve the few clusters where the two methods disagree". Best for PBMC / whole blood, where the cell-type vocabulary is small and known; novel tissue still needs manual review.

**Part 2: Canonical marker dot plot**
- Dot plot of all canonical PBMC markers across clusters
- Use this to manually assign `CLUSTER_CELLTYPE_MAP` in `config.R`

**Part 3: Label assignment**
- If `CLUSTER_CELLTYPE_MAP` is set: uses manual labels
- Otherwise: uses SingleR majority-vote per cluster as fallback

**Part 4 (contamination override)**
- After majority-vote assignment, per-cell SingleR labels for types listed in `CONTAMINATION_TYPES` (Neutrophil, RBC, HSPC, Platelet, Basophil, Eosinophil, Mast cell) override the cluster label
- This ensures scattered contamination cells appear on UMAP even if they are a minority in their cluster
- Override count logged per type: `[CONTAM] Preserved N <type> cell(s)`
- `contamination_summary.pdf` — stacked bar of contamination counts per sample

**Canonical marker genes used:**

| Cell Type | Markers |
|-----------|---------|
| T cell (pan) | CD3D, CD3E |
| CD4+ T | CD4, IL7R, CCR7 |
| CD8+ T | CD8A, CD8B, GZMK |
| NK | NKG7, GNLY, KLRD1 |
| B cell | MS4A1, CD79A, CD19 |
| CD14+ Mono | CD14, LYZ, CST3, S100A8 |
| FCGR3A+ Mono | FCGR3A, MS4A7 |
| DC | FCER1A, CLEC9A |
| Platelet | PPBP, PF4 |

---

### Step 05r — Frozen-Reference Label Transfer (`05r_reference_transfer.R`)

**Input:** `integrated/integrated_annotated.rds`, a frozen model bundle (`REFERENCE_MODEL`)
**Output:** `annotation/reference_transfer_cells.csv.gz`, `annotation/reference_transfer_composition.csv`

De-novo annotation (step 05) is run-relative: clustering and the SingleR/scType consensus are
recomputed over whatever samples are in the run, so the same sample can get different proportions
in different runs, and cross-species runs mislabel (bat neutrophils collapse into CD14+ Mono at
the human Neutrophil/Monocyte decision boundary). Step 05r removes that dependency.

It loads a **frozen reference model** built once by `pipeline/build_reference.R` (a
`trainSingleR` classifier saved with its gene set and per-anchor baseline), aligns the run's genes
to the model's gene set (zero-filling any missing), and runs `classifySingleR(fine.tune = TRUE)`
to give every cell a label that depends only on the fixed model, not the run. The result is
written as a second column, **`cell_type_ref`**, alongside the de-novo `cell_type` (additive,
never overwrites). Per-cell confidence is the top1-top2 score gap (`ref_delta`).

- Gated on `REFERENCE_MODEL` (`config.R` / `SCRNA_REFERENCE_MODEL`). No model → the script
  quits cleanly and the whole feature is a no-op, so default human runs are unaffected.
- `fine.tune` is per-cell inference against a fixed model; it sharpens calls but cannot overfit
  the query and does not fix batch effects.
- Build a reference: `SCRNA_SPECIES=bat Rscript pipeline/build_reference.R <run_dir> --holdout=Aksh1,ES332`.
- Full walk-through (build → transfer → benchmark → re-render): [`docs/howto-frozen-reference.md`](docs/howto-frozen-reference.md).

---

### Step 06 — Visualisation (`06_visualize.R`)

> **Note on label source:** when a 05r transfer exists, `06` calls `apply_reference_labels()`
> and plots `cell_type_ref` (run-independent); the proportion-bar subtitle reads
> `Labels: frozen reference`. Without a transfer it plots the de-novo `cell_type`
> (`Labels: de-novo`). The de-novo call is preserved as `cell_type_denovo`.

> **Note on split-by-sample UMAP:** panels are produced at ≤2 per page. A shared cell-type legend (extracted from the full merged object via `cowplot::get_legend()`) is placed as a right column; each per-sample panel uses `NoLegend()` to prevent duplicate legends.

**Input:** `integrated_annotated.rds`  
**Output:** `results/integrated/*.pdf`

All figures saved at 300 DPI as PDF.

| File | Description |
|------|-------------|
| `umap_triptych.pdf` | 3-panel: cluster / sample / cell type |
| `integrated_umap_cluster.pdf` | UMAP coloured by Seurat cluster |
| `integrated_umap_sample.pdf` | UMAP coloured by H1 vs H2 |
| `integrated_umap_celltype.pdf` | UMAP coloured by cell type |
| `umap_split_by_sample.pdf` | Side-by-side H1 vs H2 UMAP |
| `harmony_before_after.pdf` | Batch effect before/after Harmony |
| `feature_T_cells.pdf` | Feature plots: T cell markers |
| `feature_NK.pdf` | Feature plots: NK markers |
| `feature_B_cells.pdf` | Feature plots: B cell markers |
| `feature_monocytes.pdf` | Feature plots: monocyte markers |
| `feature_DC_platelet.pdf` | Feature plots: DC + platelet markers |
| `feature_all_canonical_markers.pdf` | All markers in one composite |
| `integrated_dotplot.pdf` | Dot plot: cell type × canonical markers |
| `integrated_heatmap.pdf` | Heatmap: top 3 markers per cluster |
| `celltype_proportions_bar.pdf` | Cell type proportions per sample |
| `celltype_counts_bar.pdf` | Cell type counts per sample |
| `violin_key_markers.pdf` | Violin plots for 8 lineage markers |

---

### Step 09 — Bootstrap Proportion Comparison (`09_bootstrap_proportions.R`)

**Input:** `integrated_annotated.rds`
**Output:** `reports/bootstrap_proportions_report.pdf`, `reports/bootstrap_summary.csv`

Normalises cell type proportion estimates across samples with different capture depths:

1. **Observed proportions** — per sample with multinomial 95% CI (Wilson-score per cell type)
2. **Bootstrap downsampling** — 1000 draws at n_min (smallest sample size) without replacement; removes capture-depth bias
3. **Per-cell-type CI plot** — dot + error bar per sample × cell type (faceted)
4. **Pairwise chi-squared** — all sample pairs tested; tile heatmap with p-values

> Note: bootstrap mean ≈ observed proportion (same point estimate, different uncertainty). The value of bootstrap is the CI at a normalised depth, not a shifted mean.

Run: `bash pipeline/run_pipeline.sh <samples> 09`

---

### Step 10 — Rarefaction Analysis (`10_rarefaction.R`)

**Input:** `integrated_annotated.rds`
**Output:** `reports/rarefaction_report.pdf`, `reports/rarefaction_summary.csv`

Empirically determines minimum capture depth for stable proportion estimates:

1. Uses the **largest sample** as ground truth
2. Subsamples at increasing depths (100, 250, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000)
3. 1000 bootstrap draws per depth; computes 95% CI width and RMSE vs ground truth
4. Fits theoretical CI ~ a/√n per cell type
5. Reports minimum n where CI ≤ 2× asymptotic value

**Findings for bat whole blood (ES03_newkit ground truth, n=6,520):**

| Precision target | Min cells needed |
|---|---|
| ±5% (crude) | ~500 |
| ±2% (moderate) | ~2,000 |
| ±1% (good) | ~4,000–5,000 |
| ±0.5% (high) | 6,000+ |

**Recommendation: ≥5,000 cells per bat whole blood sample** for full ±1% precision across all cell types.

Run: `bash pipeline/run_pipeline.sh <samples> 10`

---

### Step 08 — Comparison Report (`08_comparison_report.R`)

**Input:** existing output files in `RESULTS_DIR` (QC CSVs, doublet PDFs, composition CSVs, integrated PDFs, DE files)  
**Output:** `Comparison_report.pdf` in `RESULTS_DIR`

Standalone cross-sample summary report with sections:

| Section | Source files |
|---------|-------------|
| Cover page | pipeline metadata |
| Sample Quality & Cell Fate | `qc/cell_fate.csv`, QC violins/scatter PDFs |
| Doublet Rates | doublet UMAP + histogram PDFs per sample |
| Cell Type Composition | `integrated/celltype_composition_combined.pdf` |
| Integrated UMAP | all pages of `integrated/umap_split_by_sample.pdf` |
| Key Marker Expression | `integrated/integrated_dotplot.pdf`, violin PDFs |

Run it alone: `bash run_pipeline.sh <samples> 08`

---

### Step 08c — Cross-Run Benchmark (`08c_benchmark_concordance.R`)

**Input:** `annotation/reference_transfer_cells.csv.gz` (from 05r), the frozen model's baseline
**Output:** `benchmark/concordance.csv`, `benchmark/wholeblood_signature.csv`, `benchmark/benchmark_report.md`

Once labels are run-independent (05r), the **anchor** samples carried into every run become a
real reproducibility control. Step 08c compares each anchor's `cell_type_ref` composition in THIS
run against the model's stored baseline (the same anchors classified when the model was built).
Identical input + same model should reproduce, so any drift beyond `DRIFT_FLAG_PP` (default 5 pp)
is a pipeline/run artifact, not biology. In practice the anchors reproduce within ~2.3–2.4 pp
across the ES03-batch and ES17-batch integrations, which is what resolved the original cross-run
proportion discrepancy.

It also prints a per-sample **whole-blood signature** (Neutrophil / RBC / Platelet %): high values
indicate presort / whole blood, depleted values indicate postsort PBMC. The generated
`benchmark_report.md` includes a collection-method table and biological caveats (e.g. the Aksh1
cardiac-puncture outlier; why a de-novo 0% neutrophil result is a scRNA artifact, not a sort).

> **Statistical note:** two anchor samples are a valid reproducibility *control* (a fixed
> reference point), not a biological population *benchmark* (no variance, cannot generalise).

Gated on `REFERENCE_MODEL`; self-skips if 05r has not run. Design: `docs/frozen_reference_scope.md`.

---

## Configuration Reference (`config.R`)

All parameters are centralised in `config.R`. Every script sources this file at startup.

### Key parameters to adjust

```r
# QC thresholds
QC$min_features   <- 200    # lower = keep more cells
QC$max_features   <- 5000   # raise to 6000 if legitimate cells exceed this
QC$max_percent_mt <- 20     # lower to 15 for stricter filtering

# Clustering resolution
CLUSTER$default_res <- 0.5  # raise for more clusters, lower for fewer

# Harmony diversity penalty
HARMONY$theta <- 2           # raise to 3-4 if samples don't mix; lower to 1 if over-integrated

# Manual annotation (fill after inspecting canonical_markers_dotplot.pdf)
CLUSTER_CELLTYPE_MAP <- c(
  "0" = "CD4 T",
  "1" = "CD14+ Mono",
  # ...
)
```

### Config validation (`validate_config.R`)

Runs before step 01 (and standalone via `Rscript pipeline/validate_config.R`). Hard-errors on `CLUSTER_CELLTYPE_MAP` colour/key problems; **warns** (does not halt) on missing `SAMPLE_PATHS`, so an unmounted NAS doesn't block config-only validation — pass `--strict-paths` to escalate. Full check table: [`docs/reference-config.md`](docs/reference-config.md).

### Cache invalidation engine

Steps 01–03 cache their per-sample output under `sample_cache/<name>/` with a `.hash` sidecar. A step recomputes only when `cache_hash(nm, step)` changes:

```
cache_hash = md5( cumulative_params(step) + resolved_input_path + matrix_fingerprint )
```

- **cumulative_params** are per-step and nested (01: `QC` + species; 02: + `DOUBLET`; 03: + `NORM`/`DIM`/`CLUSTER`), so changing a threshold vector invalidates only the steps that depend on it — a clustering tweak does not bust the QC cache.
- **resolved_input_path** keys on the absolute sample path, so two experiments that share a folder name never collide in `sample_cache/`.
- **matrix_fingerprint** (size + mtime of the 10x files) busts the cache when the input matrix changes, even if config is byte-identical.

Steps 04+ are not cached, so `HARMONY` / `MARKERS` / `SINGLER_REF` changes take effect on every run. Delete `sample_cache/<name>/` to force a full recompute for one sample.

### Parallelism

Two RAM-governed worker pools: `PARALLEL$workers` for the per-sample phase (01–03, 8 GB/worker) and `PARALLEL$merge_workers` for the merged-object phase (04–06b, 16 GB/worker → fewer workers, since each holds a copy of the full merged object). A governor reads `MemAvailable` and caps both so `workers × budget` leaves ~20% headroom. OMP/OpenBLAS/MKL/BLAS are pinned to 1 thread per worker to stop CPU thrashing.

### Frozen-reference labels & benchmark

| Parameter | Env | Default | Effect |
|---|---|---|---|
| `REFERENCE_MODEL` | `SCRNA_REFERENCE_MODEL` | `""` (off) | Path to the model bundle from `build_reference.R`. Set it to turn on `05r` (label transfer) and `08c` (benchmark). Empty = both self-skip. |
| `ANCHOR_SAMPLES` | `SCRNA_ANCHORS` | `Aksh1,ES332` | Control samples carried into every run; their drift vs baseline is the benchmark metric. |
| `DRIFT_FLAG_PP` | `SCRNA_DRIFT_PP` | `5` | Anchor drift (percentage points) above which 08c flags a run/pipeline artifact. |
| `REF_FINE_TUNE` | — | `TRUE` | Per-cell `classifySingleR` fine-tuning in 05r (sharpens calls; cannot overfit the query). |

`apply_reference_labels(obj, run_dir)` (in `config.R`) is the shared helper that steps 06 and 09
call after loading the object: if `<run_dir>/annotation/reference_transfer_cells.csv.gz` exists it
swaps `cell_type` to `cell_type_ref` (keeping the de-novo call as `cell_type_denovo`) and records
`options(scrna.label_source)` so plots can label which view they show. Idempotent; a no-op when no
transfer exists, so non-bat runs are unchanged.

`SCRNA_RESULTS_DIR` points any step at a finished run directory (overrides the
samples-derived `RESULTS_DIR`), so you can re-render its PDFs/reports without reconstructing
every `SCRNA_SAMPLE*` var. Steps 06 and 07 recover the true sample list from the object and the
run-dir name rather than the config fallback. Example:

```bash
SCRNA_SPECIES=bat SCRNA_RESULTS_DIR=Results/results_<run>_filtered \
  Rscript pipeline/06_visualize.R && \
SCRNA_SPECIES=bat SCRNA_RESULTS_DIR=Results/results_<run>_filtered \
  Rscript pipeline/09_bootstrap_proportions.R
```

---

## Bat (*Eonycteris spelaea*) Whole-Blood Support

Pass `bat` as the first argument to `run_pipeline.sh` to activate species-specific overrides:

```bash
bash pipeline/run_pipeline.sh bat /path/ES03 /path/ES12
```

The `bat` keyword exports `SCRNA_SPECIES=bat`; `config.R` reads it and automatically applies:

| Setting | Human default | Bat override |
|---------|--------------|-------------|
| `SINGLER_REF` | `HumanPrimaryCellAtlas` | `MonacoImmune` |
| `CLUSTER$default_res` | 0.5 | 1.0 |
| `CLUSTER$resolutions` | 0.3–0.8 | 0.3, 0.5, 0.8, 1.0 |
| `MARKERS$CD14_mono` | CD14, LYZ, CST3, S100A8 | CD14, LYZ, S100A8, S100A9, **CSF1R** (CST3 absent; CSF1R added) |
| `MARKERS$FCGR3A_mono` | FCGR3A, MS4A7 | FCGR2A, FCGR3B, MS4A7 (FCGR3A absent) |
| `MARKERS$Neutrophil` | FCGR3B, CSF3R, CXCR2, CEACAM8 | FCGR3B, CSF3R, CXCR2, CEACAM6, **IDO1** (CEACAM8 absent; IDO1 bat-specific) |
| `MARKERS$gamma_delta_T` | — | TRDC, TRGC1, TRGC2 |
| `CONTAMINATION_TYPES` | Neutrophil, RBC, HSPC, Platelet, … | Basophil, Eosinophil, Mast cell only |
| B cell subtype markers | Uses IGHD, IGHM, IGHG1 | Uses TCL1A, IL4R, CD24, FCER2 (isotypes absent) |

**Bat reference GTF:** `Bat/ESpe_merged_fullref.gtf` — 30,180 genes (30,011 original + 169 new TCR loci merged from `ESpe_Peaks2UTRed_genome_tcellgenes_ZF_Dec2024_v2.gtf`). The new GTF has only 207 `gene`-level feature lines (TCR genes only); the remaining ~30k genes are represented in transcript/exon/CDS records — this does not affect GEM quantification since STAR counts via exon features.

**Expected bat whole-blood cell types** (4-sample run Sample6/Sample7/10/ES03_newkit, res=0.5, MonacoImmune):

| Cell type | Typical range | Notes |
|---|---|---|
| Monocyte | 23–80% | Dominant; bat whole blood is constitutively myeloid-high |
| CD4 T | 7–44% | Wide inter-individual variation; Sample10 outlier (44%) |
| CD8 T | 0.4–20% | Confirmed by CD3E 90.7%, NCAM1 0.5% — SingleR calls NK incorrectly |
| Neutrophil | 1–10% | CSF3R 82%, CXCR2 50%, IDO1 44% — genuine; human references under-call |
| DC | 1–3% | FLT3 73%, HLA-DR 98% — SingleR absorbs into Monocyte without override |
| B cell | 0.7–4% | Low in bat whole blood; consistent across samples |
| HSPC | 1–4% | Haematopoietic progenitors; normal in bat whole blood |
| NK (true) | 0–1% | NCAM1/CD56+; minor population |

**Mandatory CLUSTER_CELLTYPE_MAP overrides for bat whole blood** (verify cluster numbers per run):
- NK-labelled clusters with CD3E >80% → **CD8 T** (CTL, not NK)
- Monocyte-labelled clusters with FLT3 >40% and HLA-DR >70% → **DC**
- HSPC-labelled clusters with S100A8/A9 100% and CSF3R >30% → **Neutrophil**

**Minimum recommended capture depth:** ≥5,000 cells per sample for ±1% CI on all populations (rarefaction analysis using ES03_newkit as ground truth).

### Sample type: unsorted whole blood (presort)

The bat ES cohorts are **unsorted whole blood** (presort): not FACS-sorted and not
density-gradient (Ficoll/PBMC) separated, so granulocytes, RBC, and platelets are retained
(typically RBC-lysed, not density-separated). Confirmed for both the ES03 and ES17 batches via the
frozen-reference whole-blood signature (high Neutrophil + RBC + Platelet move together).

This is why the de-novo annotation's **0% neutrophils** in the ES17 cohort was an artifact, not
biology: the 24,737-cell neutrophil cluster was mislabeled CD14+ Mono (cross-species SingleR vote
flip), and droplet scRNA also systematically under-captures granulocytes. The frozen bat reference
recovers them (ES18 63%, ES171 59%). High neutrophil fractions in captive bat whole blood are
biologically credible: adult pteropodids are neutrophil-dominant, and captivity plus handling
stress amplify the count. Full citations: [`docs/bat_neutrophil_literature.md`](docs/bat_neutrophil_literature.md).

**Collection-method caveat:** `Aksh1` was drawn by terminal cardiac puncture, the other samples by
conscious vein draw. The terminal bleed avoids stress neutrophilia and may pool splenic/central
lymphocytes, so Aksh1 reads lowest-neutrophil / highest-lymphoid (B ~15%). It is valid as a
reproducibility anchor but not as a biological baseline; `ES332` is the representative anchor. Draw
volume does not change scRNA proportions (fixed cell loading); collection site and stress do. B-cell
fractions of 0.1–7% in the vein-draw samples are normal-to-low for neutrophil-dominated whole blood
(closed-sum: abundant neutrophils suppress the lymphoid percentage).

---

## Re-running Individual Steps

```bash
conda activate scrna_seurat

# Re-run a single step
bash /data/alvin/scRNA/pipeline/run_pipeline.sh 05

# Re-run annotation and visualisation only (most common after manual labelling)
bash /data/alvin/scRNA/pipeline/run_pipeline.sh /path/S1 /path/S2 05 06 07

# Re-run from doublet detection onwards
bash /data/alvin/scRNA/pipeline/run_pipeline.sh /path/S1 /path/S2 02 03 04 05 06 07

# Re-generate comparison report only (all earlier steps already done)
bash /data/alvin/scRNA/pipeline/run_pipeline.sh /path/S1 /path/S2 08

# Re-run visualisation + both reports
bash /data/alvin/scRNA/pipeline/run_pipeline.sh /path/S1 /path/S2 06 07 08
```

**Run any single step directly (no wrapper).** Every step resolves `config.R` relative to its own location, so you can run one stage with `Rscript` and your own env vars — handy when iterating on a single stage:

```bash
SCRNA_SAMPLE1=/path/S1 SCRNA_SAMPLE2=/path/S2 SCRNA_SPECIES=human \
  Rscript pipeline/05_annotate.R
```

If you omit `SCRNA_SAMPLE*` / `SCRNA_SPECIES`, `config.R` logs a warning and falls back to the hardcoded H1/H2 + `human` defaults — it never silently runs on the wrong inputs.

> Steps 01–03 skip automatically via the per-sample, per-step cache (`sample_cache/<name>/`; see **Cache invalidation engine**). Editing the input matrix, or any parameter a cached step depends on, auto-invalidates the relevant caches. Delete `sample_cache/<name>/` to force a full reprocess.

---

## Saved R Objects

| File | Content | Size |
|------|---------|------|
| `results/individual/H1/H1_filtered.rds` | H1 post-QC Seurat object | ~5 MB |
| `results/individual/H1/H1_singlets.rds` | H1 post-doublet-removal | ~5 MB |
| `results/individual/H1/H1_seurat.rds` | H1 fully processed (UMAP, clusters) | ~15 MB |
| `results/individual/H2/H2_*.rds` | Same for H2 | — |
| `results/integrated/integrated_seurat.rds` | Merged + Harmony integrated | ~30 MB |
| `results/integrated/integrated_annotated.rds` | Final annotated object | ~30 MB |

All paths above are relative to the run's results directory, `Results/results_<samples>_<filtered|raw>/` (e.g. `Results/results_H1-H2_filtered/`). Step 05 also writes `annotation/consensus_annotation.csv` and `annotation/cluster_annotation_table.csv`.

Load any object for interactive exploration:
```r
library(Seurat)
seu <- readRDS("Results/results_H1-H2_filtered/annotation/integrated_annotated.rds")
DimPlot(seu, group.by = "cell_type")
```
