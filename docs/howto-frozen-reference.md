# How to build and use a frozen reference

This guide takes you from a finished, well-annotated run to **run-independent cell-type labels**
and a **cross-run benchmark**. You will build a frozen reference once, point new runs at it, and
read the concordance report that tells you whether shared samples reproduce across runs.

Use this when de-novo annotation is run-relative: the same sample gets different cell-type
proportions in different runs, or a cross-species run mislabels a population (on bat whole blood,
neutrophils collapse into CD14+ Mono at the human Neutrophil/Monocyte boundary, so a sample can
read 0% neutrophils in one run and 25% in another). The frozen reference removes that dependency
by classifying every cell against a fixed model instead of the run's own clustering.

For the design rationale and locked decisions, see
[`frozen_reference_scope.md`](frozen_reference_scope.md). For the biology behind high bat-blood
neutrophils, see [`bat_neutrophil_literature.md`](bat_neutrophil_literature.md).

## Prerequisites

- The conda env from `pipeline/setup_env.sh`, activated: `conda activate scrna_seurat`.
- A finished run to serve as the reference source: a `Results/results_*_filtered/` directory that
  has `integrated/integrated_annotated.rds` and a cell-type annotation you trust. For bat whole
  blood the canonical source is the ES03 8-sample run.
- One or more samples that appear in both the reference run and your later runs. These are the
  **anchors** (default `Aksh1,ES332`); they are the control that the benchmark checks.

## Step 1: Build the reference (once)

`build_reference.R` trains a SingleR classifier from the annotated run, keeping only cells whose
label is corroborated by its defining markers, and saves a model bundle.

```bash
SCRNA_SPECIES=bat Rscript pipeline/build_reference.R \
  Results/results_Aksh1-ES03-ES14-ES258-ES332-ES35-ES407-ES459_filtered \
  --holdout=Aksh1,ES332
```

`--holdout` excludes the anchor samples from training so their baseline is an honest, held-out
estimate. Other flags: `--out=DIR` (output dir), `--min-cells=N` (drop labels with fewer than N
cells, default 20), `--validate-thr=X` (marker-validation cutoff, default 0.25),
`--cd8-thr=X` (CD8 recovery threshold, default 0.5), `--label-col=COL` (which metadata column to
train on).

This writes a model bundle plus audit files under `Results/frozen_reference/`:

```
frozen_ref_<run>_bat_<date>_v2.rds    # the model + meta (gene set) + per-anchor baseline (~96 MB)
frozen_ref_<run>_bat_<date>_v2.meta.txt
frozen_ref_<run>_bat_<date>_v2.validation.csv   # per-label mean defining-marker expression
frozen_ref_<run>_bat_<date>_v2.baseline.csv     # anchor proportions under the model
```

The `.rds` is gitignored (too large to commit); keep it on disk and reference it by path.

## Step 2: Point runs at the reference

Set `SCRNA_REFERENCE_MODEL` to the bundle path. On a full pipeline run, `05r` then runs right
after annotation (so the PDFs and HTML pick up the run-independent labels) and `08c` runs at the
end:

```bash
SCRNA_SPECIES=bat \
SCRNA_REFERENCE_MODEL=Results/frozen_reference/frozen_ref_<run>_bat_<date>_v2.rds \
  bash pipeline/run_pipeline.sh bat /path/A /path/B
```

Both steps self-skip when `SCRNA_REFERENCE_MODEL` is unset, so a default human run never touches
them. Optional knobs: `SCRNA_ANCHORS` (control samples, default `Aksh1,ES332`), `SCRNA_DRIFT_PP`
(drift-flag threshold, default 5).

## Step 3: Transfer labels onto an already-finished run

To add frozen-reference labels to a run that already completed, run `05r` then `08c` against its
directory. No clustering or integration is repeated:

```bash
MODEL=Results/frozen_reference/frozen_ref_<run>_bat_<date>_v2.rds

SCRNA_SPECIES=bat SCRNA_REFERENCE_MODEL=$MODEL \
  Rscript pipeline/05r_reference_transfer.R Results/results_<your-run>_filtered

SCRNA_SPECIES=bat SCRNA_REFERENCE_MODEL=$MODEL \
  Rscript pipeline/08c_benchmark_concordance.R Results/results_<your-run>_filtered
```

`05r` accepts `--model=PATH` and `--no-fine-tune`; `08c` accepts `--model=PATH`, `--anchors=A,B`,
and `--drift=N`. `05r` with fine-tune takes roughly 30 to 40 minutes on ~70k cells.

## Step 4: Re-render the PDFs and HTML with the new labels

The transfer writes a second label column; the figures need a re-render to show it. Use
`SCRNA_RESULTS_DIR` to point the rendering steps at the finished run without re-listing samples:

```bash
SCRNA_SPECIES=bat SCRNA_RESULTS_DIR=$PWD/Results/results_<your-run>_filtered \
  Rscript pipeline/06_visualize.R

SCRNA_SPECIES=bat SCRNA_RESULTS_DIR=$PWD/Results/results_<your-run>_filtered \
  Rscript pipeline/09_bootstrap_proportions.R

# finalized 05-Integrated_report.pdf + Overall_report.pdf
SCRNA_SPECIES=bat SCRNA_RESULTS_DIR=$PWD/Results/results_<your-run>_filtered \
  Rscript pipeline/07_finalize_reports.R

# interactive HTML report
Rscript pipeline/08b_html_report.R Results/results_<your-run>_filtered
```

The proportion PDFs now print `Labels: frozen reference` in the subtitle. The HTML report defaults
its bars and UMAP to `cell_type_ref` with a **Frozen reference vs De-novo** toggle and adds a
**Frozen-reference benchmark** section.

## Verification

Confirm the transfer and benchmark produced their outputs:

```bash
ls Results/results_<your-run>_filtered/annotation/reference_transfer_cells.csv.gz
cat Results/results_<your-run>_filtered/benchmark/benchmark_report.md
```

The report's verdict line should read `PASS (max anchor drift X pp <= 5 pp)`. A PASS means the
shared anchor samples reproduce their frozen-reference composition within tolerance across runs,
which is the whole point: the same input plus the same model gives the same proportions. If you
re-rendered the figures, the proportion-bar subtitle should say `Labels: frozen reference`, and
populations the de-novo view missed (e.g. neutrophils on cross-species data) should now appear.

## Troubleshooting

| Symptom | Cause / fix |
|---------|-------------|
| `05r` prints "no REFERENCE_MODEL set" and exits | `SCRNA_REFERENCE_MODEL` is empty or the path is wrong. Pass the full `.rds` path. |
| `05r` errors on `rownames(test)` vs `rownames(ref)` | Gene mismatch. `05r` aligns to the model's stored gene set and zero-fills missing genes automatically; if you see this, the model bundle is stale or hand-edited. Rebuild it with Step 1. |
| `08c` prints "no reference_transfer_cells.csv.gz" and skips | Run `05r` first; `08c` reads its output. |
| `08c` prints "no ANCHOR_SAMPLES set" | Set `SCRNA_ANCHORS` (or `--anchors=A,B`) to samples present in this run. |
| Benchmark verdict is DRIFT (> 5 pp) | A real run/pipeline difference, not biology, since input and model are fixed. Check the flagged cell types in `concordance.csv`; inspect whether the anchor was processed differently this run. |
| PDFs still show de-novo labels after `05r` | The figures were rendered before the transfer existed. Re-run Step 4. |
| `06`/`07` crash with "No cells found" under `SCRNA_RESULTS_DIR` | You are on an old build; `06`/`07` now recover the sample list from the object and run-dir name. Pull latest. |
| cDC2 looks inflated vs the de-novo call | Known residual: the cDC2 reference profile is weak and acts as a partial attractor on transfer. Anchors still pass; treat minor-type ref proportions with care. |

## Related

- [`frozen_reference_scope.md`](frozen_reference_scope.md) — design, locked decisions, acceptance test
- [`bat_neutrophil_literature.md`](bat_neutrophil_literature.md) — why high bat-blood neutrophils are credible
- [DOCUMENTATION.md](../DOCUMENTATION.md) — Step 05r, Step 08c, and the config reference
- [ReportGuide.md](../ReportGuide.md) — reading `cell_type_ref` and the benchmark report
