# Getting Started: Your First scRNA-seq Analysis

You will run a complete two-sample scRNA-seq analysis — from raw 10x Genomics matrices to a full annotated UMAP and PDF reports. By the end you will have a `Results/` directory containing cell type labels, composition plots, and a curated summary report.

Expected time: **15–30 minutes** (depending on CPU cores and dataset size).

---

## What you'll need

- Linux system (Ubuntu 20.04 or later)
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed
- ~8 GB RAM (32 GB recommended for datasets > 10,000 cells)
- Two 10x Genomics sample folders, each containing a `filter_matrix/` or `filtered_feature_bc_matrix/` subfolder with `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`
- This repository cloned to your machine

---

## Step 1: Install the environment

```bash
cd /path/to/scRNA
bash pipeline/setup_env.sh
```

This installs R 4.3.3 and all required packages into a conda environment named `scrna_seurat`. Takes **10–20 minutes** on first install.

When it finishes you will see:

```
✓ Environment 'scrna_seurat' created successfully.
```

Activate the environment:

```bash
conda activate scrna_seurat
```

---

## Step 2: Verify your sample folders

```bash
ls /path/to/SampleA/filter_matrix/
# barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz

ls /path/to/SampleB/filter_matrix/
# barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz
```

Both folders must have all three files. If your folder is named `filtered_feature_bc_matrix/` instead of `filter_matrix/`, the pipeline detects both automatically.

---

## Step 3: Run the pipeline

```bash
bash pipeline/run_pipeline.sh /path/to/SampleA /path/to/SampleB
```

For bat whole-blood samples, add the species flag first:

```bash
bash pipeline/run_pipeline.sh bat /path/to/ES03 /path/to/ES12
```

The pipeline prints progress as it runs:

```
[run_pipeline] SCRNA_SPECIES=human
[run_pipeline] Samples: SampleA SampleB
[run_pipeline] Results directory: /path/to/scRNA/Results/results_SampleA-SampleB_filtered
--- Step 01: Load & QC ---
  SampleA: 4,231 cells after QC (removed 312)
  SampleB: 3,889 cells after QC (removed 201)
--- Step 02: Doublets ---
--- Step 03: Per-sample normalisation & clustering ---
...
--- Step 07: Finalise reports ---
  Overall_report.pdf
DONE Step 07 finished in 847s
```

You have results as soon as Step 03 finishes (~5–10 minutes) — the per-sample UMAPs are already written and readable.

---

## Step 4: Open the QC report

```bash
xdg-open Results/results_SampleA-SampleB_filtered/reports/01-QC_report.pdf
```

Check the violin plots on the first two pages. You are looking for:

- **Gene count (nFeature_RNA):** A unimodal peak around 1,000–3,000 genes. A bimodal distribution suggests empty droplets weren't fully filtered.
- **UMI count (nCount_RNA):** Mirrors gene count. Should peak around 2,000–8,000.
- **Mitochondrial %:** Most cells below 10–15%. A high-MT tail (> 20%) means dying cells passed the filter — lower `QC$max_percent_mt` in `config.R` and re-run step 01.

If QC looks good, proceed.

---

## Step 5: Review annotation and fill the cluster map

```bash
# See what SingleR suggested
grep -A 30 "CLUSTER_CELLTYPE_MAP" \
  Results/results_SampleA-SampleB_filtered/logs/05_annotate.log
```

Open the dot plot to verify the suggestions:

```bash
xdg-open Results/results_SampleA-SampleB_filtered/annotation/canonical_markers_dotplot.pdf
```

For each cluster, the dot plot shows which canonical marker genes are expressed. Cross-check against the SingleR suggestions. A cluster SingleR calls `"NK"` but with high CD3E expression is a CD8 T cell — override it.

Edit `pipeline/config.R`:

```r
CLUSTER_CELLTYPE_MAP <- c(
  "0"  = "CD4 T (naive)",
  "1"  = "NK",
  "2"  = "CD14+ Mono",
  # ... paste the full suggested map here, correcting wrong entries
)
```

Save `config.R`, then re-run annotation and reporting:

```bash
bash pipeline/run_pipeline.sh /path/to/SampleA /path/to/SampleB 05 06 07
```

This takes **2–5 minutes** (no re-processing of QC or clustering).

---

## Step 6: Open the final report

```bash
xdg-open Results/results_SampleA-SampleB_filtered/reports/Overall_report.pdf
```

This is a curated A4 summary of every stage. Each figure has a bold title banner and a `Good: … | Bad: …` interpretation caption at the bottom.

Key pages to check:

| Page title | What to look for |
|-----------|-----------------|
| UMAP — Cell Types | Named clusters, no large "Unknown" blob |
| Dot Plot — Canonical Markers × Cell Type | Each cell type expressing its expected markers |
| Cell Type Proportions | Biologically sensible composition per sample |
| SingleR Delta Score | Most cells above 0.1 delta (confident annotation) |

---

## What you built

You now have in `Results/results_SampleA-SampleB_filtered/`:

- `annotation/integrated_annotated.rds` — the annotated Seurat object; load this in R for any downstream analysis
- `reports/Overall_report.pdf` — shareable summary report
- `reports/01-05_report.pdf` — detailed per-stage reports
- `integrated/celltype_proportions_bar.pdf` — cell type composition per sample

To continue:

- Run additional steps 08–10 for comparison reports, bootstrap proportions, and rarefaction analysis
- Adjust QC thresholds and re-run from step 01 if QC plots showed issues
- Add more samples: see [How to Add a Sample](howto-add-sample.md)
- Correct any remaining mislabelled clusters: see [How to Override Annotations](howto-override-annotations.md)

---

## Troubleshooting

**`conda activate scrna_seurat` fails: no such environment**
Run `bash pipeline/setup_env.sh` first. If mamba is not installed, the script uses conda — installation will take longer but produce the same result.

**Step 01 fails: `No matrix found in /path/to/SampleA`**
The pipeline looks for `filter_matrix/` or `filtered_feature_bc_matrix/` inside the sample folder. Check that the path you passed points to the sample root (not to the matrix subfolder itself). Correct: `/data/ES03`; Wrong: `/data/ES03/filter_matrix`.

**All cells labelled "Unknown" after annotation**
`CLUSTER_CELLTYPE_MAP` contains a cell type name that does not match a key in `CELLTYPE_COLORS`. Check spelling exactly against the `CELLTYPE_COLORS` keys in `config.R`.

**Run takes > 60 minutes for < 5,000 cells**
SingleR is the slowest step. Check that `WORKERS` in `config.R` is greater than 1 (`parallel::detectCores() - 2`). On a single-core VM, the pipeline falls back to sequential execution.
