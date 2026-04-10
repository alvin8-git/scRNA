# TODO

## Done

- [x] **Samba share** — installed samba 4.15.13, config written; pending 3 activation commands (append config, `smbpasswd -a alvin`, restart smbd/nmbd). Map from Windows: `\\192.168.1.168\alvin`
- [x] **Pipeline: dynamic sample folders** — `run_pipeline.sh` now accepts 1 or 2 folder paths as args; sample names auto-derived from directory names; `SINGLE_SAMPLE` mode skips Harmony
- [x] **Pipeline: combined PDF reports** — 4 final reports in `results/`: `report_doublets.pdf`, `report_individual.pdf`, `report_annotation.pdf`, `report_integrated.pdf`; each figure has a bold title banner; small plots paired 2-per-portrait-page
- [x] **Pipeline: remove hardcoded H1/H2** — all scripts loop over `SAMPLE_NAMES` dynamically; new samples just need folder paths passed to `run_pipeline.sh`
- [x] **New step 07** — `07_finalize_reports.R` merges partial PDFs into the 4 category reports using qpdf/pdfunite/gs

---

## Immediate — Activate Samba Share

```bash
sudo bash -c 'cat /tmp/samba_alvin_share.conf >> /etc/samba/smb.conf'
sudo smbpasswd -a alvin
sudo systemctl restart smbd nmbd && sudo systemctl enable smbd nmbd
```
Then map `\\192.168.1.168\alvin` from Windows.

---

## Immediate — Manual Cell Type Annotation

The pipeline ran with SingleR automated labels as a fallback. T cells (59.7%) are not yet split into CD4+/CD8+ subtypes. Complete manual annotation:

- [ ] Open `results/annotation/canonical_markers_dotplot.pdf`
- [ ] Open `results/annotation/singler_labels_umap.pdf` (reference)
- [ ] Open `results/annotation/cluster_annotation_table.csv`
- [ ] Map each of the 7 clusters to a specific PBMC cell type
- [ ] Fill `CLUSTER_CELLTYPE_MAP` in `pipeline/config.R`:
  ```r
  CLUSTER_CELLTYPE_MAP <- c(
    "0" = "CD4 T",
    "1" = "CD8 T",
    "2" = "NK",
    # ... complete for all 7 clusters
  )
  ```
- [ ] Re-run: `bash /data/alvin/scRNA/pipeline/run_pipeline.sh 05 06`
- [ ] Verify `results/integrated/integrated_umap_celltype.pdf` looks biologically correct

---

## Analysis Improvements

- [ ] **Split CD4/CD8 T cells** — after manual annotation, consider sub-clustering the T cell cluster at higher resolution to separate CD4 naive, CD4 memory, CD8 effector, regulatory T cells
- [ ] **Monocyte sub-typing** — check if both CD14+ and FCGR3A+ monocyte populations are present (currently merged as "Monocyte" by SingleR)
- [ ] **Dendritic cells** — low abundance; may be absent in H1 (only 359 cells); check FCER1A/CLEC9A expression
- [ ] **Tune Harmony integration** — inspect `harmony_before_after.pdf`; if H1/H2 still separate, increase `HARMONY$theta` from 2 to 3
- [ ] **Adjust clustering resolution** — currently `res = 0.5` gives 7 clusters; consider `res = 0.6` or `res = 0.8` to better separate T cell subtypes

---

## Differential Expression

- [ ] Compare H1 vs H2 within each cell type (donor effect vs biology)
- [ ] Run `FindMarkers(group.by = "sample")` per cell type to find donor-specific genes
- [ ] Check for inflammatory signature differences between H1 and H2

---

## QC Review

- [ ] Review `results/qc/H1_violin_qc.pdf` — H1 has only 359 post-QC cells; confirm thresholds are not too aggressive
- [ ] Check `results/doublets/H1_doublet_score_hist.pdf` — score distribution should be bimodal (clear singlet peak)

---

## Figures for Publication / Report

- [ ] Finalise cell type labels (depends on manual annotation above)
- [ ] Export `umap_triptych.pdf` as main figure
- [ ] Export `integrated_dotplot.pdf` as supplementary marker validation figure
- [ ] Export `celltype_proportions_bar.pdf` for sample comparison
- [ ] Consider converting PDFs to PNG for PowerPoint compatibility:
  ```bash
  conda activate scrna_seurat
  Rscript -e "
  library(Seurat); library(ggplot2)
  seu <- readRDS('results/integrated/integrated_annotated.rds')
  p <- DimPlot(seu, group.by='cell_type')
  ggsave('results/integrated/umap_celltype.png', p, width=8, height=7, dpi=300)
  "
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

- [ ] Install `qpdf` R package for PDF merging in step 07: `install.packages("qpdf")`
- [ ] Add `scDblFinder` to `setup_env.sh` so it installs via conda on fresh environments:
  ```bash
  mamba install -n scrna_seurat -c conda-forge -c bioconda bioconductor-scdblfinder -y
  ```
- [ ] Remove DoubletFinder GitHub install from `setup_env.sh` (replaced by scDblFinder)
