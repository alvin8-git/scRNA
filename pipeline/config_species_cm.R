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
  # min_features 200 -> 1000 (2026-09-10). The 200 floor let a large ambient-dominated
  # population through: cluster 0 of the first run held 23,498 cells (24.8%) at a median
  # of 1,141 genes, and its top FindAllMarkers genes were all abundant secreted plasma
  # proteins (RBP4/TTR/FGB/AHSG/APOC3/APOA1/APOA2/AFP) plus MALAT1 — the ambient
  # signature, not a lineage. Restricting to >2000-gene cells retained only 4.9% of that
  # cluster while the genuine ALB+/HNF4A+ hepatic cluster retained 70.1%, and D20
  # composition changed from 51% "endoderm" to Epicardial 43.4% / Fibroblast 39.1% /
  # endoderm 2.7%. SoupX would be the principled fix but needs the raw unfiltered
  # matrices, and Samples/Cardio ships only filter_matrix.
  #
  # Cost, stated because it is a real bias: the loss is not uniform across the time
  # course — 20.6% of D0, 5.6% of D11, 25.7% of D20 and 2.0% of D30. Comparisons across
  # timepoints are therefore made on differently-filtered populations. 1000 was chosen
  # over 1500 (which removes 68.7% of D20) and over 800 (too lenient: D20's Q2/Q3 sit at
  # 1,123-1,337 genes and remain ambient-dominated).
  QC$min_features <- 1000
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

  # --- H1 cardiomyocyte differentiation (SCRNA_SPECIES=cm), res 0.5, 21 clusters -------
  # Guarded on the exact sample set: cluster numbering is not stable across runs, so this
  # must never leak into a different cohort. Curated 2026-09-10 from
  # integrated/integrated_cluster_markers.csv (FindAllMarkers) — the marker-panel scores
  # alone were unusable because ambient collagen (COL1A1/DCN ~2-3 in every cluster) made
  # every cluster look fibroblast-like. EVERY cluster is mapped: a partial map falls back
  # to per-cell SingleR labels for the rest, which would shatter them.
  #
  # SingleR (HumanPrimaryCellAtlas) called clusters 2/6/7/10/11/12/14/15 "Neurons" — that
  # is a reference artefact, not biology. Verified TUBB3 0.00, MAP2 0.00-0.02, ELAVL3 0.01,
  # SOX2 0.00, PAX6 0.00, SOX10 0.00 across all of them: there are no neurons in this
  # culture. HPCA has no cardiomyocyte label, so it files excitable cells under Neurons.
  # The QC$min_features == 200 term is NOT redundant with the sample-name guard. Cluster
  # numbering depends on which cells survive QC, so raising min_features renumbers every
  # cluster while the sample set stays identical — the name guard alone would silently
  # apply these stale numbers to the new run. Curated under min_features 200; when a run
  # uses a different QC floor this map correctly withholds itself and step 05 falls back
  # to auto-annotation, which prints a fresh paste-ready map to re-curate from.
  if (length(SAMPLE_NAMES) == 7 && isTRUE(QC$min_features == 200) &&
        setequal(SAMPLE_NAMES, c("H1D0_1", "H1D0_2", "H1D11_2", "H1D20_1",
                               "H1D20_2", "H1D30_1", "H1D30_2"))) {
    CLUSTER_CELLTYPE_MAP <- c(
      # --- cardiac lineage -----------------------------------------------------------
      "10" = "Cardiomyocyte",      # MYH6/TTN/ACTC1/ACTN2/MYL7/MYOCD/SLC8A1/LDB3/CCDC141;
                                   # TTN 2.55 vs <=0.84 elsewhere, NKX2-5 0.63, MEF2C 0.57
      "19" = "Cardiomyocyte",      # PLN/MYL3/HSPB6/CRYAB/SMIM3 — the most mature CM here
                                   # (98% D30); subtype refinement should call it
      "12" = "Cardiac progenitor", # GATA4 1.51 (highest), TBX5 0.37, TECRL (cardiac-
                                   # specific), ITGA8, CCBE1, LIX1; 46% D11
      # --- epicardium / mesothelium --------------------------------------------------
      "18" = "Epicardial",         # ITLN1 0.87 (unique), TBX18, ALDH1A2, NPR3, SFRP5, UPK3B
      "3"  = "Epicardial",         # UPK3B/SFRP2/PTGDS/SLPI/NPY mesothelial signature
      "14" = "Epicardial",         # same signature as 3 (SPRR2F/UPK3B/SLC7A7) + high MT
      # --- stromal -------------------------------------------------------------------
      "4"  = "Fibroblast",         # FMOD/COL6A3/FBN1/DLK1/LOX/SERPINE2 — definitive
      "2"  = "Fibroblast",         # CNTN5/TENM2/SOX6/PDE3A/ZFPM2; sarcomere-negative
                                   # (TTN 0.71, TNNT2 0.20) and neural-negative
      "7"  = "Fibroblast",         # same programme as 2, 49% D20
      # --- off-target endoderm (the largest single lineage) --------------------------
      "0"  = "Hepatic/Endoderm",   # RBP4/TTR/FGB/AHSG/APOC3/APOA1/APOA2/AFP — visceral
                                   # /yolk-sac endoderm; 23,498 cells, 70% D20
      "13" = "Hepatic/Endoderm",   # ALB/APOB/MTTP/CEBPA/F2/AMN — hepatocyte-like
      "11" = "Hepatic/Endoderm",   # HNF4A/HNF1A-AS1/ONECUT1/HHEX/FOXA2/NR5A2
      "20" = "Hepatic/Endoderm",   # FOXA2/HHEX/ONECUT1/FOXA1 + cell cycle; 98% D11, n=64
      # --- pluripotent ---------------------------------------------------------------
      "5"  = "Pluripotent",        # UTF1/NANOG/SOX2/ALPL/POU5F1/LNCPRESS1; 91% D0
      "8"  = "Pluripotent",        # DPPA4/L1TD1/MIR302CHG/ESRG/XACT; 96% D0
      "9"  = "Pluripotent",        # XACT/CADM2/GRID2/RMST; 94% D0
      "1"  = "Pluripotent",        # POU5F1/ESRG/DPPA4/MIR302CHG/CRABP1
      # --- other ---------------------------------------------------------------------
      "6"  = "Proliferating",      # KIF20A/PBK/MKI67/NEK2/TOP2A/CDCA3/TPX2 — pure cycle,
                                   # no lineage genes in its top markers
      "15" = "Epithelial",         # GABRP/CLDN4/CLDN7/GRHL2/PRSS8/RAB25/WFDC2
      "16" = "Endothelial",        # CDH5/ICAM2/TIE1/ESAM/SOX7/ECSCR/CD34/GJA4; 65% D11
      "17" = "Unknown"             # top markers are ALL MT- genes — mito-high/dying,
                                   # 54% D20 (the shallow libraries). Do not interpret.
    )
  }

  # --- Same cohort at min_features 1000, res 0.5, 20 clusters ------------------
  # Curated 2026-09-10 from the re-run's integrated_cluster_markers.csv. Guarded on the
  # QC floor as well as the sample names, so this and the min_features 200 map above are
  # mutually exclusive and neither can apply to the other's clustering.
  #
  # Raising min_features to 1000 did NOT dissolve the ambient population: cluster 0 still
  # holds 17,581 cells (21.3%) at a median of 1,216 genes, i.e. sitting just above the
  # new floor. It is labelled Unknown, not Hepatic/Endoderm, on direct evidence — it
  # carries the secreted CARGO without the lineage IDENTITY:
  #
  #            FOXA2 HNF4A ONECUT1 HHEX CEBPA | ALB  AFP  TTR APOA1 | COL1A1 DCN TAGLN
  #   cl0       0.02  0.02   0.01  0.01  0.03 | 1.06 2.75 1.82 2.18 |  2.60 2.57  1.37
  #   cl9  real 1.04  0.40   0.48  0.63  0.12 | 0.49 1.06 1.77 1.23 |    -    -     -
  #   cl12 real 0.82  0.40   0.21  0.29  0.46 | 1.61 2.89 3.19 3.84 |    -    -     -
  #
  # Ambient RNA carries abundant secreted transcripts, not transcription factors, and
  # cl0 is simultaneously collagen-high, plasma-protein-high AND ACTA2/TAGLN-positive —
  # the average of the whole culture rather than any one lineage. Calling it Unknown
  # costs 21% of the run but does not invent a population. SoupX on the raw matrices is
  # the real fix; Samples/Cardio ships only filter_matrix (see docs).
  if (length(SAMPLE_NAMES) == 7 && isTRUE(QC$min_features == 1000) &&
      setequal(SAMPLE_NAMES, c("H1D0_1", "H1D0_2", "H1D11_2", "H1D20_1",
                               "H1D20_2", "H1D30_1", "H1D30_2"))) {
    CLUSTER_CELLTYPE_MAP <- c(
      # --- cardiac lineage ---------------------------------------------------------
      "7"  = "Cardiomyocyte",      # NKX2-5 0.64, MEF2C 0.62, GATA4 1.11, TTN 2.72,
                                   # ACTC1 2.06, TNNT2 1.12 + MYH6/MYOCD/SLC8A1/LDB3/CMYA5
      "10" = "Proepicardial",      # was "Cardiac progenitor". % cells: TBX18 17, WT1 25, TCF21 27,
                                   # TBX5 38, GATA4 88, PDGFRA 49 + LHX2/SFRP5/C7/HGF/COLEC11
                                   # top markers; NKX2-5 only 13 (vs 38 in CM cl7). Trajectory
                                   # places it at the tip of the stromal branch, not before CMs.
      # --- epicardium / mesothelium ------------------------------------------------
      "3"  = "Epicardial",         # UPK3B/SFRP2/PTGDS/SLPI/NPY/SLC34A2
      "11" = "Epicardial",         # same programme (SPRR2F/SLPI/UPK3B/SLC7A7)
      "17" = "Myofibroblast",      # TAGLN 2.47 / ACTA2 1.80 with WT1 0.23, TBX18 0.26,
                                   # ALDH1A2 0.28 — epicardium-derived, plus ANKRD1/CCN2/CCN1
      # --- stromal -----------------------------------------------------------------
      "5"  = "Fibroblast",         # FMOD/COL6A3/FBN1/LOX/DLK1; COL1A1 3.74
      "1"  = "Fibroblast",         # placeholder only — cluster 1 is split by
                                   # CLUSTER_SUBCLUSTER_MAP below; this label is never
                                   # the final one for any of its cells
      "6"  = "Fibroblast",         # COL1A1 2.63/COL3A1 3.23/DCN 2.35/LUM 2.55
      # --- hepatic endoderm (the genuine fraction) ---------------------------------
      "9"  = "Hepatic/Endoderm",   # FOXA2 1.04/HHEX 0.63/NR5A2 0.61/ONECUT1 0.48/HNF1A-AS1
      "12" = "Hepatic/Endoderm",   # APOB/MTTP/ALB/CEBPA 0.46/PLG/AMN/CIDEB
      "19" = "Hepatic/Endoderm",   # ALB 3.31/A2M/FGA/FGG/FGB/FABP1/SERPINA1; n=58, 100% D30
      # --- pluripotent -------------------------------------------------------------
      "2"  = "Pluripotent",        # UTF1/SOX2/GAL/TDGF1/FOXD3-AS1/LNCPRESS1; 88% D0
      "4"  = "Pluripotent",        # POU5F1/MIR302CHG/ESRG/CRABP1
      "16" = "Pluripotent",        # POU5F1 1.75/SOX2 0.79/XACT; 96% D0
      # --- other -------------------------------------------------------------------
      "8"  = "Proliferating",      # KIF20A/NEK2/MKI67/TOP2A/CDCA8/ASPM — pure cycle
      "14" = "Epithelial",         # GABRP/CLDN4/CLDN7/GRHL2/PRSS8/RAB25/WFDC2
      "15" = "Endothelial",        # CDH5/ICAM2/TIE1/ESAM/SOX7/CD34/GJA4; 67% D11
      # --- not interpretable -------------------------------------------------------
      "0"  = "Unknown",            # ambient-dominated, see the note above (21.3%)
      "13" = "Unknown",            # MT-genes dominate the markers, 1,217 genes, MT 2.1
      "18" = "Unknown"             # MT 6.7% and 1,211 genes — dying
    )
    # Cluster 1 (9,555 cells) mixed two populations: labelled wholesale as Fibroblast it
    # put 22.9% of D0 — undifferentiated hESC — into "Fibroblast". FindSubCluster at
    # res 0.2 on RNA_snn separates it cleanly (verified 2026-09-14):
    #   1_0  n=3,530  87% D0  POU5F1 1.14 DNMT3B 1.96 L1TD1 1.41 LIN28A 0.85, COL1A1 0.72
    #   1_1  n=2,534  D11-D30 COL1A1 2.43 COL3A1 2.79 DCN 1.67 LUM 2.10, POU5F1 0.07
    #   1_2  n=2,234  D20-D30 COL1A1 2.42 COL3A1 2.99 DCN 1.94 POSTN 0.88, POU5F1 0.05
    #   1_3  n=  882  63% D30 COL1A1 3.06 COL3A1 3.31 DCN 1.97, POU5F1 0.03
    #   1_4  n=  375  85% D11 COL1A1 2.03 COL3A1 1.70, DCN 0.34 — early mesenchyme
    # Res 0.4 gives the same pluripotent/fibroblast split in 7 pieces; 0.2 is preferred.
    CLUSTER_SUBCLUSTER_MAP <- list(
      cluster    = "1",
      resolution = 0.2,
      labels     = c("0" = "Pluripotent", "1" = "Fibroblast", "2" = "Fibroblast",
                     "3" = "Fibroblast",  "4" = "Fibroblast")
    )
  }

  # --- SoupX-corrected cohort (*_sx), curated 2026-09-23 -----------------------
  # Same 7 libraries after SoupX adjustCounts (rho 0.076-0.673; H1D20_1 lost 67.3% of
  # counts and 76.9% of its cells, confirming it was ambient-dominated). Correction
  # re-clusters the data, so the uncorrected map above MUST NOT be reused: different
  # sample names keep this block and that one mutually exclusive.
  # Curated from per-cluster detection rates (22 clusters, 61,697 cells).
  if (length(SAMPLE_NAMES) == 7 && isTRUE(QC$min_features == 1000) &&
      setequal(SAMPLE_NAMES, c("H1D0_1_sx", "H1D0_2_sx", "H1D11_2_sx", "H1D20_1_sx",
                               "H1D20_2_sx", "H1D30_1_sx", "H1D30_2_sx"))) {
    CLUSTER_CELLTYPE_MAP <- c(
      # --- pluripotent: POU5F1/SOX2/LIN28A/L1TD1 high, D0-dominant ---------------
      "6"  = "Pluripotent",      # POU5F1 82 SOX2 68 L1TD1 77; 94% D0
      "8"  = "Pluripotent",      # POU5F1 86 SOX2 81 LIN28A 88; 87% D0
      "7"  = "Pluripotent",      # POU5F1 78 SOX2 54; 79% D0
      "4"  = "Pluripotent",      # POU5F1 31 SOX2 26; 76% D0, lower content (1,973 genes)
      # --- cardiomyocyte ----------------------------------------------------------
      "9"  = "Cardiomyocyte",    # TTN 93 MYL7 83 MYL4 72 MYH6 68 ACTC1 68 TNNT2 61,
                                 # NKX2-5 41 MEF2C 44 TBX5 37 GATA4 66; 47% D11
      # --- proepicardial: TBX18/WT1/TBX5/SFRP5/LHX2 with GATA4, low NKX2-5 -------
      "16" = "Proepicardial",    # GATA4 91 TBX5 45 WT1 27 SFRP5 23 LHX2 20 TBX18 17
      "14" = "Proepicardial",    # GATA4 69 TBX5 29 TBX18 27 WT1 26 SFRP5 43 UPK3B 30
      # --- epicardium / mesothelium ----------------------------------------------
      "0"  = "Epicardial",       # UPK3B 62 with COL1A1 93 DCN 91 POSTN 60 KDR 62
      # --- stromal ----------------------------------------------------------------
      "5"  = "Myofibroblast",    # ACTA2 76 TAGLN 78 on COL1A1 99 LUM 97 PDGFRB 54
      "10" = "Fibroblast",       # PDGFRA 60 PDGFRB 53 TCF21 31 LUM 92 DCN 73
      "17" = "Fibroblast",       # PDGFRA 47 PDGFRB 50 TCF21 34 GATA4 77 COL1A1 98
      "3"  = "Fibroblast",       # COL1A1 88 COL3A1 92 LUM 82 DCN 74 POSTN 48
      "11" = "Fibroblast",       # COL1A1 76 DCN 61 LUM 62 PDGFRB 35
      # --- hepatic endoderm: definitive TFs, not cargo alone ----------------------
      "13" = "Hepatic/Endoderm", # FOXA2 75 HNF4A 50 HHEX 35 ONECUT1 25, TTR 96 APOA1 93
      "12" = "Hepatic/Endoderm", # FOXA2 64 HHEX 47 ONECUT1 42 HNF4A 41, EPCAM 65
      "21" = "Hepatic/Endoderm", # FOXA2 94 HHEX 69 ONECUT1 50; n=54, cycling, 100% D11
      # --- other ------------------------------------------------------------------
      "18" = "Endothelial",      # CDH5 75 KDR 85 PECAM1 58
      "15" = "Epithelial",       # EPCAM 75 GABRP 54 CLDN4 55 KRT8 90
      "2"  = "Proliferating",    # TOP2A 90 MKI67 67 (POU5F1 30 — cycling pluripotent)
      # --- unresolved: low content + secreted cargo, no identity TFs --------------
      "1"  = "Unknown",          # 1,151 genes / 2,155 UMI (lowest); AFP 88 ALB 69 with
                                 # COL1A1 95 LUM 96 MYL7 75 — cargo without identity
      "20" = "Unknown",          # 1,207 genes / 2,393 UMI; AFP 78 ALB 65; n=121
      "19" = "Unknown"           # n=212, 100% D30; MYL7 99 ACTC1 81 AND AFP 100 ALB 99
                                 # — CM markers and hepatic cargo together; likely doublets
    )
  }

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
