# How to Override Cell Type Annotations

After step 05 runs, review the suggested cluster labels and correct any that are wrong. This is the most impactful manual step in the pipeline.

## Prerequisites

- Step 05 has completed successfully
- `Results/<run>/logs/05_annotate.log` exists
- `Results/<run>/annotation/canonical_markers_dotplot.pdf` exists

---

## Steps

### 1. Read the suggested map from the log

```bash
grep -A 30 "CLUSTER_CELLTYPE_MAP" Results/results_ES03-ES12_filtered/logs/05_annotate.log
```

You will see output like:

```r
CLUSTER_CELLTYPE_MAP <- c(
  "0"  = "CD4 T (naive)",
  "1"  = "NK",
  "2"  = "CD14+ Mono",
  "3"  = "B cell (memory)",
  "4"  = "NK",           # <-- review this one
  "5"  = "CD8 T",
  ...
)
```

### 2. Open the canonical markers dot plot

```bash
xdg-open Results/results_ES03-ES12_filtered/annotation/canonical_markers_dotplot.pdf
# or on macOS:
open Results/results_ES03-ES12_filtered/annotation/canonical_markers_dotplot.pdf
```

For each cluster flagged by SingleR, check which genes are expressed. Use these as ground truth:

| Marker genes high | Correct label |
|-------------------|---------------|
| CD3D, CD3E + CD8A, CD8B | `"CD8 T"` |
| CD3D, CD3E + CD4, IL7R | `"CD4 T"` |
| GNLY, NKG7, KLRD1 + **no CD3** | `"NK"` |
| CD3D high + GNLY, NKG7 high | `"NKT"` |
| S100A8, S100A9, LYZ, CD14 | `"CD14+ Mono"` |
| CDKN1C, FCGR3A, MS4A7 | `"FCGR3A+ Mono"` |
| MS4A1, CD79A, TCL1A | `"B cell"` |
| PPBP, PF4, ITGB3 | `"Platelet"` |
| S100A12, BST1, G0S2 | `"Neutrophil"` |
| HBB, HBA1, HBA2 | `"RBC"` |
| FCER1A, CLEC10A, FLT3 | `"cDC2"` |

A common mislabelling: SingleR calls a cluster `"NK"` but the dot plot shows CD3E > 80%. That is a cytotoxic T cell (`"CD8 T"`), not NK.

### 3. Edit `pipeline/config.R`

Open `pipeline/config.R`. Paste the suggested map from the log and correct the entries you identified:

```r
CLUSTER_CELLTYPE_MAP <- c(
  "0"  = "CD4 T (naive)",
  "1"  = "NK",
  "2"  = "CD14+ Mono",
  "3"  = "B cell (memory)",
  "4"  = "CD8 T",          # corrected from NK
  "5"  = "CD8 T"
  # clusters not listed fall back to SingleR automatically
)
```

You do not need to list every cluster. Omit clusters where SingleR was correct — they fall back to SingleR automatically.

### 4. Re-run steps 05–07

```bash
bash pipeline/run_pipeline.sh bat /path/to/ES03 /path/to/ES12 05 06 07
```

### 5. Verify the correction

Open `Results/<run>/annotation/annotation_umap.pdf`. The corrected cluster should now show the right label. Cross-check with `canonical_markers_dotplot.pdf` again to confirm marker-label alignment.

---

## Verification

```bash
# Count cells per final cell type
conda run -n scrna_seurat Rscript -e '
  source("pipeline/config.R")
  seu <- readRDS(file.path(DIRS$annotation, "integrated_annotated.rds"))
  print(sort(table(seu$final_cell_type), decreasing=TRUE))
'
```

Expected output shows named cell type counts with no obvious labelling errors.

---

## Troubleshooting

**Cluster numbers changed after re-running step 03**
Cluster numbers are not stable across runs. Never copy a `CLUSTER_CELLTYPE_MAP` from a different dataset or a previous run at a different resolution. Always start from the log of the current run.

**`"Unknown"` cells appearing after override**
A cell type name in `CLUSTER_CELLTYPE_MAP` does not match a key in `CELLTYPE_COLORS`. Check spelling exactly — `"CD8 T"` not `"CD8T"` or `"CD8T cell"`. The canonical names are the keys in `CELLTYPE_COLORS` in `config.R`.

**All B cells disappeared**
`"B cell"` is a contamination type in `CONTAMINATION_TYPES` for some species overrides. If you want B cells labelled by cluster majority (not per-cell SingleR), remove `"B cell"` from `CONTAMINATION_TYPES` in `config.R`.

---

## Related

- [Annotation Strategy](explanation-annotation.md) — why two routes exist and how they interact
- [How to Re-run from a Specific Step](howto-rerun-steps.md) — re-running steps 05–07
- [Configuration Reference](reference-config.md) — `CLUSTER_CELLTYPE_MAP`, `CONTAMINATION_TYPES`
