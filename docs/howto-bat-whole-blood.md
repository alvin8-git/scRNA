# How to Run on Bat Whole Blood

This guide covers running the pipeline on *Eonycteris spelaea* whole-blood samples. The `bat` species keyword activates a set of pre-tuned overrides that address the biology of bat whole blood.

## Prerequisites

- Conda environment active: `conda activate scrna_seurat`
- 10x Genomics matrices aligned to the *E. spelaea* genome (ESpe GTF; HGNC-style gene names via orthology mapping)
- Each sample folder contains a `filter_matrix/` or `filtered_feature_bc_matrix/` subfolder

---

## What the `bat` keyword changes

Passing `bat` as the first argument to `run_pipeline.sh` sets `SCRNA_SPECIES=bat`, which applies these overrides in `config.R`:

| Parameter | Human default | Bat override |
|-----------|--------------|--------------|
| `SINGLER_REF` | `"HumanPrimaryCellAtlas"` | `"MonacoImmune"` |
| `CLUSTER$default_res` | 0.5 | 1.0 |
| `MARKERS` | Human PBMC markers | Bat-validated orthologues |
| `REFINEMENT_MARKERS` | Human CD4/CD8 sub-types | Bat-validated T cell sub-types |
| `CONTAMINATION_TYPES` | Includes Neutrophil | Neutrophil stays (abundant in whole blood) |
| γδ T markers | Not included | TRGC2, TRGC1 added |

`MonacoImmune` is used because it resolves 29 blood cell subtypes including γδ T cells, monocyte subtypes, and pDC/mDC — critical for whole-blood samples where these populations are present.

Resolution 1.0 is used because whole blood contains granulocytes, RBCs, and platelets that form tight clusters requiring higher resolution to split from lymphoid populations.

---

## Steps

### 1. Verify your matrices use HGNC gene names

```bash
zcat /path/to/ES03/filter_matrix/features.tsv.gz | head -5
```

Expected: `ENSG00000... CD3D Gene Expression` — gene names (column 2) should be HGNC symbols, not *E. spelaea* gene IDs. If they are Ensembl IDs, the canonical marker lookup will silently fail.

### 2. Run the pipeline with the `bat` keyword

```bash
conda activate scrna_seurat
bash pipeline/run_pipeline.sh bat /path/to/ES03 /path/to/ES12
```

For a single sample:
```bash
bash pipeline/run_pipeline.sh bat /path/to/ES03
```

For more than two samples:
```bash
bash pipeline/run_pipeline.sh bat /path/to/ES03 /path/to/ES14 /path/to/ES332
```

### 3. First run: leave `CLUSTER_CELLTYPE_MAP = NULL`

Do not set `CLUSTER_CELLTYPE_MAP` on the first run. Let the pipeline auto-annotate via SingleR (MonacoImmune). After step 05, review the log:

```bash
grep -A 30 "CLUSTER_CELLTYPE_MAP" Results/results_ES03-ES12_filtered/logs/05_annotate.log
```

### 4. Review annotation with bat biology in mind

Open `annotation/canonical_markers_dotplot.pdf` and `annotation/singler_delta_umap.pdf`.

Expected bat whole-blood cell types and their approximate proportions:

| Cell type | Approximate % | Key markers |
|-----------|--------------|-------------|
| Neutrophil | 20–60% | S100A12, S100A9, BST1 |
| CD14+ Mono | 15–35% | S100A8, LYZ, CD14 |
| CD4 T | 5–25% | CD3D, CD4, IL7R |
| NK | 5–15% | GNLY, NKG7, KLRD1 |
| RBC | 5–15% | HBB, HBA1 |
| Platelet | 1–10% | PPBP, PF4 |
| B cell | 0.1–2% | MS4A1, CD79A |

B cells are rare in *E. spelaea* whole blood — proportions below 2% are expected (not a pipeline error). For meaningful B cell analysis, PBMC isolation before library prep is required.

### 5. Common annotation corrections for bat

SingleR on bat whole blood frequently makes these mistakes:

| SingleR label | Actual label | Evidence |
|---------------|-------------|---------|
| `"NK"` | `"CD8 T"` | CD3E > 80%, NCAM1 < 1% |
| `"NK"` | `"γδ T"` | TRGC2 or TRGC1 > 20% |
| `"Monocyte"` | `"DC"` | FCER1A > 10%, FLT3 > 5% |
| `"Pre-B cell"` | `"Neutrophil"` | S100A12, G0S2 high |

Correct these in `CLUSTER_CELLTYPE_MAP` and re-run steps 05–07:

```bash
bash pipeline/run_pipeline.sh bat /path/to/ES03 /path/to/ES12 05 06 07
```

### 6. Interpret the composition plots

Open `integrated/celltype_proportions_bar.pdf`. For whole-blood samples:

- High Neutrophil % (30–60%) is normal — this is not a QC failure
- Samples with low Neutrophil % (< 10%) likely have more lymphoid cells and better represent adaptive immunity
- The `_500umi` sample variants (lower minUMI threshold) will show higher Neutrophil % due to admission of low-RNA granulocytes

---

## Verification

After step 07 completes:

```bash
ls -lh Results/results_ES03-ES12_filtered/reports/Overall_report.pdf
```

Open the report and check:
- Page 1 (QC): cells kept per sample, %MT below 20
- Annotation pages: Neutrophil, CD14+ Mono, CD4 T clusters visible on UMAP
- Composition bar: Neutrophil dominates as expected

---

## Troubleshooting

**Most cells labelled `"Unassigned"`**
The MonacoImmune reference could not match bat transcriptomes. Check that gene names in your matrices are HGNC symbols (step 1). If they are bat-specific IDs, the orthology mapping step in alignment was skipped.

**B cells absent from UMAP**
B cells are genuinely rare in bat whole blood. This is expected. If you need B cells, run PBMC isolation (density gradient centrifugation) before library prep.

**Error: `object 'MARKERS' not found`**
The `bat` species override was not applied. Verify `SCRNA_SPECIES=bat` by checking the first lines of `logs/05_annotate.log` for the species confirmation message.

---

## Related

- [How to Override Annotations](howto-override-annotations.md) — correcting bat-specific mislabellings
- [Configuration Reference](reference-config.md) — all bat-override parameters
- [Annotation Strategy](explanation-annotation.md) — understanding contamination-type overrides
