# How to Add or Replace a Sample

Adding a new sample to an existing cohort only processes the new sample through steps 01–03; it loads existing samples from cache and re-runs integration and annotation on the full set.

## Prerequisites

- Conda environment active: `conda activate scrna_seurat`
- Existing completed run under `Results/results_<samples>_filtered/`
- New sample folder with a `filter_matrix/` or `filtered_feature_bc_matrix/` subfolder

---

## Steps

### 1. Confirm the new sample folder structure

```bash
ls /path/to/NewSample/
# Expected: filter_matrix/ (or filtered_feature_bc_matrix/)

ls /path/to/NewSample/filter_matrix/
# Expected: barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz
```

### 2. Add a color for the new sample (optional)

Open `pipeline/config.R` and add an entry to `SAMPLE_COLORS`:

```r
SAMPLE_COLORS <- c(
  H1        = "#E64B35",
  H2        = "#4DBBD5",
  NewSample = "#7E6148"   # add this line
)
```

If you skip this, the pipeline auto-assigns a colour from `hue_pal()`.

### 3. Run the full pipeline with all samples including the new one

```bash
# Original run was H1 + H2; now adding NewSample
bash pipeline/run_pipeline.sh /path/to/H1 /path/to/H2 /path/to/NewSample
```

The runner exports `SCRNA_SAMPLE1=/path/to/H1`, `SCRNA_SAMPLE2=/path/to/H2`, `SCRNA_SAMPLE3=/path/to/NewSample`. The results directory becomes `Results/results_H1-H2-NewSample_filtered/`.

**H1 and H2 are loaded from `sample_cache/` — only NewSample runs steps 01–03.** The pipeline then re-runs integration (step 04) and all downstream steps with all three samples.

### 4. Set `CLUSTER_CELLTYPE_MAP = NULL` before re-running

Adding a new sample changes the integrated embedding and may change cluster numbers. Open `pipeline/config.R` and set:

```r
CLUSTER_CELLTYPE_MAP <- NULL
```

After step 05, read the new log to get an updated map:

```bash
grep -A 30 "CLUSTER_CELLTYPE_MAP" Results/results_H1-H2-NewSample_filtered/logs/05_annotate.log
```

Then add corrections and re-run steps 05–07 as needed.

---

## Replacing a sample

To swap out a sample (e.g. replace a low-quality sample with a resequenced version):

1. Delete the old sample's cache to prevent it from being loaded:
   ```bash
   rm -rf sample_cache/OldSample/
   ```
2. Run with the replacement:
   ```bash
   bash pipeline/run_pipeline.sh /path/to/H1 /path/to/NewSample
   ```

The new sample is processed fresh through steps 01–03.

---

## Verification

```bash
# Confirm new results directory was created
ls Results/

# Check new sample appears in integrated object
conda run -n scrna_seurat Rscript -e '
  source("pipeline/config.R")
  seu <- readRDS(file.path(DIRS$annotation, "integrated_annotated.rds"))
  print(table(seu$sample))
'
```

The output should list all samples including the new one with cell counts.

---

## Troubleshooting

**New sample is not in the integrated object**
The sample cache for the new sample was not created. Check `logs/03_individual.log` for errors. Common cause: the matrix path was wrong or the `filter_matrix/` subfolder was missing.

**Old results directory overwritten accidentally**
The old run used `H1 + H2` → `results_H1-H2_filtered/`. The new run with three samples creates `results_H1-H2-NewSample_filtered/`. They are separate directories. Nothing is overwritten.

**Cluster numbers changed from the previous run**
Expected — adding a sample changes the Harmony embedding. Set `CLUSTER_CELLTYPE_MAP = NULL` and re-annotate from the new log.

---

## Related

- [How to Override Annotations](howto-override-annotations.md) — updating `CLUSTER_CELLTYPE_MAP` after adding a sample
- [How to Re-run from a Specific Step](howto-rerun-steps.md) — skipping early steps when only downstream changes are needed
- [Pipeline Architecture](explanation-architecture.md) — how the sample cache works
