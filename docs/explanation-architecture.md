# Pipeline Architecture

This document explains why the pipeline is structured the way it is — the decisions behind the modular design, the sample cache, parallel execution, and the two-pass report system.

---

## The problem

Single-cell RNA-seq analysis involves 10+ sequential stages, each taking minutes to hours. Running everything in a single script means:

1. A failure at step 07 forces you to rerun hours of QC, clustering, and annotation.
2. Iterating on annotation (the most manual step) requires rerunning visualisation every time.
3. Adding one new sample to an existing cohort forces reprocessing of samples you have already QCed.

---

## The approach

```
run_pipeline.sh
      │
      ├── 01_load_qc.R      ──writes──► qc/<sample>_seurat.rds
      ├── 02_doublets.R     ──writes──► doublets/<sample>_seurat.rds
      ├── 03_individual.R   ──writes──► individual/<sample>_seurat.rds
      │                                sample_cache/<sample>_seurat.rds  ← shared cache
      ├── 04_integrate.R    ──reads─── sample_cache/*  ──writes──► integrated_seurat.rds
      ├── 05_annotate.R     ──writes──► integrated_annotated.rds
      ├── 06_visualize.R    ──writes──► visualization_report.pdf
      └── 07_finalize.R     ──writes──► Overall_report.pdf
```

Each step writes a `.rds` checkpoint. If you re-run starting at step 05, steps 01–04 are skipped entirely.

---

## The sample cache

Steps 01–03 are per-sample. Their outputs are written to both the run-specific directory and `sample_cache/<sample>/`. When you run a new sample combination, step 04 checks the cache first:

```
Run: ES03 + ES12 → processes both, caches both
Run: ES03 + ES14 → ES03 loaded from cache; only ES14 processed
```

This makes adding a new sample to a cohort fast: only the new sample is processed; integration and annotation re-run on the full set.

**Trade-off:** The cache assumes QC parameters did not change between runs. If you lower `min_counts`, delete the sample's cache subfolder to force reprocessing with the new thresholds.

---

## Harmony for batch correction

When multiple samples are provided, UMAP clusters by sample rather than biology without correction. Harmony corrects this by embedding cells into a shared low-dimensional space where the sample dimension is minimised.

```
Without Harmony:          With Harmony:
  UMAP clusters by         UMAP clusters by
  sample identity          cell type
  ────────────             ────────────
  ●●●● (ES03)              ●○●○ (CD4 T)
  ○○○○ (ES12)              ●○●○ (NK)
```

**Trade-off:** Harmony can over-correct if `theta` is too high, merging biologically distinct populations. The default `theta=2` is appropriate for samples from the same tissue type. Raise to 3–5 only for strong batch effects (different sequencing platforms or tissue sources).

The pipeline saves a UMAP before and after Harmony (`integration_umap_before.pdf`, `integration_umap_after.pdf`) so you can verify the correction visually.

---

## Two annotation routes

Step 05 assigns cell types via two independent routes:

**Route 1 — SingleR (reference-based):**
Scores every cell against a curated reference transcriptome. Fast and unbiased, but constrained to cell types that exist in the reference. Good for identifying major immune populations.

**Route 2 — Cluster majority vote + manual override:**
Aggregates per-cell SingleR labels to a per-cluster label. The majority label wins. A `CLUSTER_CELLTYPE_MAP` lets you correct obvious errors without rerunning SingleR.

**Why both?** SingleR at the per-cell level is noisy for rare populations (a Neutrophil cluster containing 3% NK contamination would be mislabelled as NK by majority vote). The contamination-types list addresses this: cells whose per-cell SingleR label matches a contamination type always keep that label, bypassing majority vote.

See [Annotation Strategy](explanation-annotation.md) for a deeper treatment.

---

## Per-plot captions in `Overall_report.pdf`

Every figure in the Overall report has a `Good: … | Bad: …` interpretation caption. These are matched by pattern from `PLOT_CAPTIONS` in `config.R`. The goal: a reader who has never seen this pipeline can open the report and immediately understand what a good vs bad result looks like for each figure, without consulting documentation.

**Trade-off:** Captions are static text defined at pipeline write-time. They do not adapt to per-run statistics. A QC violin that looks good for PBMC might look different for whole blood — the captions provide guidance, not verdicts.

---

## Results directory naming

The results directory is named from the sample names: `results_ES03-ES12_filtered/`. This means:

- Multiple runs with different sample combinations do not overwrite each other.
- You can compare `results_ES03-ES12_filtered/` and `results_ES03-ES12-ES14_filtered/` side-by-side.
- The `_filtered` or `_raw` suffix records whether filtered or raw matrices were used.

**Trade-off:** Long sample lists produce long directory names. The pipeline truncates names beyond 8 samples to avoid filesystem limits.

---

## Parallelism

Steps that process multiple samples (03, 05) use `future` multicore with a worker cap:
```r
WORKERS <- min(parallel::detectCores() - 2, 8)
```

`BiocParallel` is used for scDblFinder (step 02) and SingleR (step 05) since these libraries have their own parallelism that conflicts with `future`.

`future_mem_gb = 16` limits each worker's RAM allocation. On a 64 GB machine with 8 workers this caps total pipeline memory at 128 GB (enough for large datasets). On 32 GB machines, lower `WORKERS` to 4 or reduce `future_mem_gb`.

---

## Related

- [Configuration Reference](reference-config.md) — tuning all the parameters discussed here
- [Pipeline Steps Reference](reference-pipeline-steps.md) — what each step reads and writes
- [Annotation Strategy](explanation-annotation.md) — deep dive into the two-route annotation
