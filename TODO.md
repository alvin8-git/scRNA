# TODO

## Done

- [x] **Samba share** — installed samba 4.15.13; map from Windows: `\\192.168.1.168\alvin`
- [x] **Pipeline: dynamic sample folders** — `run_pipeline.sh` accepts any number of folder paths; names auto-derived; `SINGLE_SAMPLE` mode skips Harmony
- [x] **Pipeline: raw_matrix support** — `raw_matrix` / `raw_feature_bc_matrix` recognised alongside `filter_matrix`; results folder appends `_filtered` or `_raw` to differentiate runs
- [x] **Pipeline: combined PDF reports** — 5 numbered reports (`01-QC_report.pdf` … `05-Integrated_report.pdf`) + `Overall_report.pdf` (curated cross-stage summary); each figure has bold title banner + Good/Bad caption
- [x] **Pipeline: report naming** — reports renamed `01-QC_report.pdf` … `05-Integrated_report.pdf`; `Overall_report.pdf` combines QC + doublets + SingleR heatmap + individual HVG/markers + all integrated plots
- [x] **Pipeline: remove hardcoded H1/H2** — all scripts loop over `SAMPLE_NAMES` dynamically
- [x] **Step 07** — `07_finalize_reports.R` merges partial PDFs into category reports + `Overall_report.pdf`
- [x] **conda auto-activation** — `run_pipeline.sh` activates `scrna_seurat` automatically if not already active
- [x] **Relative path support** — `run_pipeline.sh` accepts relative paths (e.g. `H1/raw_matrix/`) as sample inputs
- [x] **Per-plot captions** — all 5 report PDFs include a description + Good/Bad interpretation guide under each figure
- [x] **ReportGuide.md** — detailed guide with plot descriptions and embedded images from H1+H2 run
- [x] **Expanded README** — QC threshold table, cell type detection explanation, marker contamination guide, annotation workflow, species adaptation table
- [x] **Extended marker list** — added Neutrophil (FCGR3B, CSF3R, CXCR2, CEACAM8), Treg (FOXP3, IL2RA, CTLA4), Plasma (MZB1, JCHAIN, SDC1), RBC (HBB, HBA1, HBA2, GYPA), HSPC (CD34, GATA2, AVP)
- [x] **Feature plot pagination** — all feature plots paginated 3 cols × 2 rows per page; dedicated Neutrophil and RBC plots added; "All Canonical Markers" replaced by "Remaining Markers (Treg/Plasma/HSPC)" to avoid duplicating genes already shown elsewhere
- [x] **Unassigned → grey** — Unassigned/Unknown forced to `#AAAAAA` in SingleR UMAP, cell type UMAP, and heatmap group bar
- [x] **Per-cell SingleR labels** — when no `CLUSTER_CELLTYPE_MAP` is set, per-cell labels used (not majority-vote per cluster) so rare types (Neutrophils, Unassigned) appear in the cell type UMAP
- [x] **Colour fixes** — Neutrophil → `#E377C2` (pink, distinct from NK teal); HSPC → `#8C564B` (brown, distinct from Monocyte salmon)
- [x] **Dot plot rotation** — gene labels rotated 90° CCW (`angle=90, hjust=1, vjust=0.5`) in integrated dot plot and individual Top 5 Markers plot (size 7)
- [x] **UMAP triptych labels** — cell type panel uses `label.size=2.5` to reduce clutter
- [x] **Cell composition plot** — stacked % bar (one bar per sample) with integer % labels; cell type name shown in large segments (≥7%); compact horizontal legend strip between panels; overall title removed from combined figure
- [x] **QC violin x-axis** — labels at 0° (horizontal)
- [x] **Doublet UMAP title** — includes count and % e.g. `H1 - Doublet UMAP (31 in 400 cells (7.8%))`
- [x] **Violin key lineage markers** — each sub-plot titled with cell type association e.g. `PPBP (Platelet)`, `CD3D (T cells)`
- [x] **Cell fate tracking** — `cell_fate.csv` saved to `results/qc/`; accounts for every GEM barcode: load losses (min.features < 200) → QC losses → doublet removal → final retained cells + %

---

## Immediate — Manual Cell Type Annotation

The pipeline ran with SingleR automated labels as a fallback. T cells (~60%) are not yet split into CD4+/CD8+ subtypes. Complete manual annotation:

- [ ] Open `results_*/annotation/canonical_markers_dotplot.pdf`
- [ ] Open `results_*/annotation/singler_labels_umap.pdf` (reference)
- [ ] Open `results_*/annotation/cluster_annotation_table.csv`
- [ ] Map each cluster to a specific PBMC cell type
- [ ] Fill `CLUSTER_CELLTYPE_MAP` in `pipeline/config.R`:
  ```r
  CLUSTER_CELLTYPE_MAP <- c(
    "0" = "CD4 T",
    "1" = "CD8 T",
    "2" = "NK",
    # ... complete for all clusters
  )
  ```
- [ ] Re-run: `bash pipeline/run_pipeline.sh H1 H2 05 06 07`
- [ ] Verify `results_*/integrated/integrated_umap_celltype.pdf` looks biologically correct

---

## Analysis Improvements

- [ ] **Split CD4/CD8 T cells** — sub-cluster T cell cluster at higher resolution to separate CD4 naive, CD4 memory, CD8 effector, Tregs
- [ ] **Monocyte sub-typing** — confirm both CD14+ and FCGR3A+ monocytes are resolved; may need `res = 0.6` or `0.8`
- [ ] **Dendritic cells** — low abundance; check FCER1A/CLEC9A expression; may be absent in small samples
- [ ] **Tune Harmony integration** — inspect `harmony_before_after.pdf`; if H1/H2 still separate, increase `HARMONY$theta` from 2 → 3–5
- [ ] **Adjust clustering resolution** — `res = 0.5` gives ~7 clusters; try `res = 0.6` or `0.8` to better separate T cell subtypes

---

## Differential Expression

- [ ] Compare H1 vs H2 within each cell type (donor effect vs biology)
- [ ] Run `FindMarkers(group.by = "sample")` per cell type to find donor-specific genes
- [ ] Check for inflammatory signature differences between H1 and H2

---

## QC Review

- [ ] Review `cell_fate.csv` — confirm removal counts at each stage per sample
- [ ] H1 has only ~360 post-QC cells; confirm QC thresholds are not too aggressive
- [ ] Check doublet score histograms — should be bimodal (clear singlet peak near 0)

---

## Figures for Publication / Report

- [ ] Finalise cell type labels (depends on manual annotation above)
- [ ] Export `umap_triptych.pdf` as main figure
- [ ] Export `integrated_dotplot.pdf` as supplementary marker validation figure
- [ ] Export `celltype_composition_combined.pdf` for sample comparison
- [ ] Consider PNG export for PowerPoint:
  ```bash
  conda activate scrna_seurat
  pdftoppm -r 150 results_*/integrated/umap_triptych.pdf umap_triptych
  ```

---

## Optional Extensions

- [ ] **Trajectory analysis** — run PAGA or Monocle3 on monocyte/DC lineage
- [ ] **Gene set scoring** — score cells for interferon response, exhaustion, activation signatures
- [ ] **Cell-cell communication** — run CellChat or NicheNet on integrated object
- [ ] **Reference mapping** — map to Azimuth PBMC reference for higher-resolution annotation
- [ ] **Batch effect validation** — compute LISI scores to quantify integration quality

---

## Environment

- [ ] Install `qpdf` R package for PDF merging: `install.packages("qpdf")`
- [ ] Add `scDblFinder` to `setup_env.sh` for fresh environment installs:
  ```bash
  mamba install -n scrna_seurat -c conda-forge -c bioconductor bioconductor-scdblfinder -y
  ```
- [ ] Remove DoubletFinder GitHub install from `setup_env.sh` (replaced by scDblFinder)
