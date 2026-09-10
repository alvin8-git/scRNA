# =============================================================================
# config_species_cm.R — overlay for human iPSC/ESC -> cardiomyocyte differentiation
# Sourced by config.R when SCRNA_SPECIES=cm. Rewrites MARKERS, ALL_MARKERS,
# SINGLER_REF, CONTAMINATION_TYPES, CLUSTER, SUBTYPE_MARKERS and QC caps.
# CELLTYPE_COLORS is NOT touched here — CM labels live in the base palette in
# config.R, matching the convention used by the bat/bat_wing overlays.
#
# Calibrated against the 7-sample H1 time course in Samples/Cardio (D0/D11/D20/D30,
# 107,658 cells, MGI DNBelab C4, GRCh38 37,488 genes) on 2026-09-10. Every gene named
# below was checked present in that annotation; percentages quoted in the comments are
# per-cell detection rates measured directly from the count matrices.
# =============================================================================

if (.species == "cm") {

  message("[Species] cm (human iPSC/ESC -> cardiomyocyte) — applying differentiation overrides")

  # ---- QC: iPSC-CMs are large, RNA-rich cells --------------------------------
  # Measured on the D0-D30 cohort: p99 nFeature reaches 8,071 and p99 nCount 54,422.
  # The human defaults (5000 / 25000) discard 12.7% of all cells and up to 21.2% of
  # H1D0_2 / H1D30_2 — real cells, not multiplets. 9000/70000 sits just above p99 for
  # every sample and drops the loss to 0.2%; scDblFinder (step 02) remains the actual
  # doublet filter, so these caps are only a backstop against extreme outliers (one
  # H1D30_2 cell reaches 435,190 UMI).
  QC$max_features <- 9000
  QC$max_counts   <- 70000
  # max_percent_mt deliberately LEFT at the base value. Cardiac tissue is
  # mitochondria-rich and a raised ceiling was expected, but this platform runs a
  # ~0.1-1.0% mito baseline (p99 = 4.6%): only 13 of 107,658 cells exceed 20%, so the
  # gate is already inert. Same DNBelab C4 property documented for the bat data.

  # ---- Merged-phase memory: this annotation is gene-rich ---------------------
  # Step 04 died at the shared 24 GiB ceiling on the first run: NormalizeData tried to
  # export 24.18 GiB, missing by 0.18. Cell count is not the driver — 94,746 cells here
  # vs 92,864 in the 8-sample bat run, but 37,488 genes vs 22,140, so each worker
  # carries ~1.7x the matrix. The bulk is Seurat's normalisation closure (FUN, 21.54
  # GiB) capturing the expression data; the Seurat object itself is only 2.64 GiB.
  # Raised here rather than in config.R so bat runs keep 3 merged workers at 24 GiB.
  # regovern() recomputes the worker count from MemAvailable so the raised ceiling does
  # not silently over-subscribe RAM (2 workers x 40 GiB = 80 of ~116 GiB available).
  PARALLEL$merge_mem_gb  <- 40L
  PARALLEL$merge_workers <- PARALLEL$regovern(40L)

  # ---- SingleR: HumanPrimaryCellAtlas, with a known blind spot ---------------
  # HPCA is the best available off-the-shelf choice here: of its 36 label.main values
  # it supplies Embryonic_stem_cells and iPS_cells (D0), Fibroblasts, Endothelial_cells,
  # Epithelial_cells, Smooth_muscle_cells, Chondrocytes, MSC, Tissue_stem_cells and
  # Hepatocytes — which covers the D0 end and the stromal/off-target populations.
  #
  # IMPORTANT: HPCA has NO cardiomyocyte label (verified 2026-09-10). The main
  # population of this entire experiment therefore CANNOT be called by SingleR; it will
  # land on Smooth_muscle_cells, Fibroblasts or MSC. This is the same class of failure
  # as bat neutrophils under MonacoImmune. Cardiomyocyte identity has to come from the
  # marker route (scType over MARKERS below, then a guarded CLUSTER_CELLTYPE_MAP built
  # from annotation/canonical_markers_dotplot.pdf). Expect a high Unassigned rate and
  # do NOT read the raw SingleR labels as the answer.
  SINGLER_REF <- "HumanPrimaryCellAtlas"

  # ---- Marker panels ---------------------------------------------------------
  # Supersedes the PBMC panel entirely — no immune cell types are expected.
  # Detection rates below are (D0 -> D11 -> D20 -> D30).
  MARKERS <- list(
    # Core sarcomere. MYL7 (0.9 -> 89.8 -> 56.1 -> 89.8) and TNNI1 (0.5 -> 73.5 ->
    # 56.2 -> 72.5) are the strongest CM discriminators in this data — both stronger
    # than MYL4 (54.5 at D30). TTN and TNNT2 confirm but are less sensitive.
    Cardiomyocyte = c("TNNT2", "TTN", "ACTC1", "ACTN2", "MYL7", "MYL4", "TNNI1",
                      "MYH6", "DES", "MYBPC3"),
    # Ventricular. Near-absent in this cohort (MYL2 1.4%, IRX4 5.0% at D30) — kept
    # precisely so the report shows that absence rather than hiding it.
    `Cardiomyocyte (ventricular)` = c("MYL2", "IRX4", "MYH7", "HEY2"),
    # Atrial / immature. NR2F2 (23.6) and NPPA (11.4) rise through D30; SLN peaks at
    # D20 (34.9) then falls to 4.7, which is worth a look.
    `Cardiomyocyte (atrial)` = c("NPPA", "NR2F2", "SLN", "GJA5", "KCNJ3"),
    # Cardiac progenitor / mesoderm — peaks sharply at D11 (NKX2-5 18.5, TBX5 19.7,
    # ISL1 15.5, GATA4 71.2) and is largely gone by D20. The D11 sample is the only
    # window on this state, and it has no replicate.
    `Cardiac progenitor` = c("NKX2-5", "ISL1", "TBX5", "GATA4", "MESP1", "KDR"),
    # Pluripotent. POU5F1 (89.2% at D0) is far more sensitive than NANOG (27.4%);
    # the source figures used only SOX2/NANOG. LIN28A persists into D11 (66.6%).
    Pluripotent = c("POU5F1", "SOX2", "NANOG", "LIN28A", "DNMT3B"),
    # Fibroblast. CCDC80 was supplied as a cardiac marker but tracks the fibroblast
    # compartment here (4.9 -> 61.5 -> 28.9 -> 68.3), so it is grouped accordingly.
    Fibroblast = c("COL1A1", "COL3A1", "DCN", "LUM", "POSTN", "TCF21", "CCDC80"),
    Myofibroblast = c("ACTA2", "TAGLN", "FN1", "POSTN"),
    `Smooth Muscle` = c("MYH11", "CNN1", "TAGLN", "ACTA2"),
    Pericyte = c("PDGFRB", "RGS5", "KCNJ8", "NOTCH3", "ANPEP"),
    # Endothelial is minimal in this cohort (PECAM1 1.1%, CDH5 1.4% at D30) despite
    # being a clear cluster in the Day-21 reference figure; APLNR peaks at D11 (15.5%).
    Endothelial = c("PECAM1", "CDH5", "APLNR", "ENG"),
    Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19"),
    # Epicardial — absent from the source figures but rising here (WT1 0.1 -> 18.4,
    # ALDH1A2 0.9 -> 9.7). A genuine expected population in cardiac differentiation.
    Epicardial = c("WT1", "TBX18", "ALDH1A2", "UPK3B"),
    # Hepatic / yolk-sac endoderm off-target. AFP goes 6.4 -> 12.8 -> 95.6 -> 99.4 and
    # FOXA2 peaks at D11 (23.4%). Detection that high is either a large endoderm
    # population or heavy ambient RNA from one — see the caveat at the foot of this file.
    `Hepatic/Endoderm` = c("AFP", "FOXA2", "SOX17", "TTR", "APOA1"),
    Proliferating = c("MKI67", "TOP2A", "CCNB1")
  )
  # Gate the integrated FindAllMarkers pass the same way the other overlays do.
  MARKERS$compute_integrated <- TRUE

  ALL_MARKERS <- unique(unlist(MARKERS[setdiff(names(MARKERS), "compute_integrated")]))

  # ---- Contamination / rare-type per-cell override ---------------------------
  # These are real but small enough that cluster majority vote would swallow them.
  # Endothelial and Pericyte are both ~1-4%; residual Pluripotent cells at D11+ are
  # the ones worth catching individually (a leftover undifferentiated pocket is a
  # quality signal, not noise). Cardiomyocyte is deliberately NOT here — it is the
  # majority population, not a contaminant.
  # Pluripotent and Proliferating were here on the first pass, but both are now assigned
  # by cluster in CLUSTER_CELLTYPE_MAP (they are 25% and 5% of the run — populations, not
  # contaminants). A per-cell override on a quarter of the data fights the curated map and
  # scatters labels; leave only the genuinely rare types the map could miss.
  CONTAMINATION_TYPES <- c("Endothelial", "Pericyte")

  # ---- Clustering ------------------------------------------------------------
  # The source figures used RNA_snn_res.0.1, which yields only 5-6 clusters. That is
  # too coarse to separate atrial/ventricular CM, epicardium and the endoderm
  # off-target across 107k cells, so the sweep starts higher. default_res 0.5 is a
  # compromise: high enough to split the stromal compartment, low enough not to
  # shatter the dominant CM population.
  CLUSTER$resolutions <- c(0.1, 0.3, 0.5, 0.8)
  CLUSTER$default_res <- 0.5
  CLUSTER$compare_res <- c(0.1, 0.3, 0.5, 0.8)

  # ---- Sub-type refinement ---------------------------------------------------
  # Splits the coarse Cardiomyocyte call by chamber identity. On this cohort the
  # expected outcome is "immature" for most clusters: MYL2 is essentially absent, so
  # a ventricular call should be treated with suspicion.
  # NOTE the naming convention: the inner names are the FULL final labels, not bare
  # suffixes — 05_annotate.R assigns `label_map[cl] <- best` verbatim, so a bare
  # "resting" would land in the metadata as the cell type (it did, on the first run)
  # and would have no CELLTYPE_COLORS entry. Match the bat overlay, which spells out
  # "CD4 T (naive)" etc.
  SUBTYPE_MARKERS <- list(
    Cardiomyocyte = list(
      "Cardiomyocyte (ventricular)" = c("MYL2", "IRX4", "MYH7", "HEY2"),
      "Cardiomyocyte (atrial)"      = c("NPPA", "NR2F2", "SLN", "GJA5"),
      "Cardiomyocyte (immature)"    = c("MYL7", "TNNI1", "MYL4", "ACTC1")
    ),
    Fibroblast = list(
      "Fibroblast (activated)" = c("POSTN", "ACTA2", "TAGLN", "FN1"),
      "Fibroblast (resting)"   = c("DCN", "LUM", "TCF21", "CCDC80")
    )
  )

  # ---- Wound modules are meaningless here ------------------------------------
  if (exists("WOUND_MODULES")) WOUND_MODULES <- NULL

  # ---- CAVEAT recorded in config so it travels with the analysis -------------
  # COL1A1 (100.0%), COL3A1 (100.0%), LUM (99.8%), DCN (99.1%) and AFP (99.4%) are
  # detected in essentially EVERY cell at D30. Near-universal detection of secreted /
  # highly abundant transcripts is the ambient-RNA signature, not biology — the same
  # pattern that made bat neutrophils unrecoverable from raw detection rates. Judge
  # these populations on soup-corrected signal (per-cluster mean minus the mean in a
  # compartment that cannot express the gene), never on detection rate alone.
  #
  # Also note H1D20_1/H1D20_2 are ~2.5-3x shallower than the other samples (median
  # 1,300 / 1,019 genes vs 3,208 at D11 and 3,508 at D30). The apparent D20 dip in
  # every CM marker is at least partly library depth, not loss of cardiomyocytes.
  # Use step 09 (depth-normalised bootstrap proportions) and step 10 (rarefaction)
  # before making any claim about D20 composition.
}
