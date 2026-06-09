# How to Re-run from a Specific Step

Each pipeline step writes checkpoint `.rds` files so you can restart from any point without reprocessing earlier stages.

## Prerequisites

- Conda environment active: `conda activate scrna_seurat`
- A completed run under `Results/results_<samples>_filtered/`
- The working directory is the repository root (`/data/alvin/scRNA` or wherever you cloned it)

---

## Steps

### 1. Identify which step to restart from

| Reason to re-run | Start from step |
|-----------------|-----------------|
| Changed `QC` thresholds | `01` |
| Changed `CLUSTER$default_res` or `DIMS` | `03` or `04` |
| Corrected `CLUSTER_CELLTYPE_MAP` | `05` |
| Changed visualisation parameters (`PLOT`, `CELLTYPE_COLORS`) | `06` |
| Added/changed plot captions | `07` |
| Regenerating reports only | `07` |

### 2. Edit `pipeline/config.R`

Open `pipeline/config.R` and make your changes. Save the file.

### 3. Re-run from the target step

Pass your sample paths followed by the step numbers to run:

```bash
# Re-run steps 05, 06, 07 on an existing bat run
bash pipeline/run_pipeline.sh bat /path/to/ES03 /path/to/ES12 05 06 07

# Re-run only step 07 (reports) — no samples needed if config already set
bash pipeline/run_pipeline.sh 07

# Re-run steps 04 through 07
bash pipeline/run_pipeline.sh /path/to/H1 /path/to/H2 04 05 06 07
```

The runner exports `SCRNA_SAMPLE1`, `SCRNA_SAMPLE2`, … and the R scripts resolve the results directory from those env vars — so `config.R` always writes to the same `Results/results_<samples>_filtered/` directory.

### 4. Check the log

```bash
tail -f Results/results_ES03-ES12_filtered/logs/05_annotate.log
```

Each step appends to its log file. The last line will read `DONE Step XX finished in NNs` when successful.

---

## Verification

After re-running step 05, confirm `integrated_annotated.rds` was updated:

```bash
ls -lh Results/results_ES03-ES12_filtered/annotation/integrated_annotated.rds
```

The timestamp should match the current run. Open the regenerated `reports/05-Integrated_report.pdf` to verify the changes are reflected.

---

## Troubleshooting

**Error: `No seurat object for '<sample>'. Run steps 01-03 first.`**
Step 04 could not find the sample in `sample_cache/`. Either run steps 01–03 first, or check that the sample path passed to `run_pipeline.sh` matches what was used in the original run (the cache subfolder is named from the sample path basename).

**Results directory mismatch**
If you pass different sample paths than the original run, the results directory name changes. Pass the same paths in the same order to write to the same directory.

**Step 05 produces different clusters after re-running step 03**
Cluster numbers are non-deterministic between runs (depends on `set.seed` in Seurat). If you changed `CLUSTER$default_res` and re-ran from step 03, the cluster numbers in any existing `CLUSTER_CELLTYPE_MAP` are now wrong. Set `CLUSTER_CELLTYPE_MAP = NULL` and re-annotate from the new log.

---

## Related

- [How to Override Annotations](howto-override-annotations.md) — the most common reason to restart at step 05
- [Pipeline Steps Reference](reference-pipeline-steps.md) — what each step reads and writes
- [Configuration Reference](reference-config.md) — all tunable parameters
