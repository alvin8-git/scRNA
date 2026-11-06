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

  # ---- Customer report: target-population panel (08b --customer) -------------
  # The read-outs a cardiac differentiation is judged by. `label` is a regex on the cell
  # label; the marker is called at >= min_umi so one ambient read does not count.
  TARGET_POPULATION <- list(
    title   = "Differentiation efficiency",
    label   = "^Cardiomyocyte", name = "cardiomyocytes", short = "CMs",
    marker  = "TNNT2", min_umi = 2,
    subtypes = c(Ventricular = "ventricular", Atrial = "atrial"),
    ratios  = list(`MYH7/MYH6` = c("MYH7", "MYH6"), `MYL2/MYL7` = c("MYL2", "MYL7")),
    desc    = paste0(
      "The read-outs a cardiac differentiation is usually judged by, per sample. ",
      "<strong>TNNT2+ cells</strong> is the single-cell equivalent of a cardiac troponin T flow ",
      "cytometry stain: a cell counts when it has at least 2 TNNT2 reads, so one stray ambient read ",
      "does not make a cardiomyocyte. It does not depend on the cell-type labels; ",
      "<strong>Labelled cardiomyocytes</strong> does, and the two should agree. Among cardiomyocytes, ",
      "<strong>MYH7/MYH6</strong> and <strong>MYL2/MYL7</strong> (summed reads) rise as cells mature ",
      "and take on a ventricular identity."))

  # ---- QC: iPSC-CMs are large, RNA-rich cells --------------------------------
  # Measured on the D0-D30 cohort: p99 nFeature reaches 8,071 and p99 nCount 54,422.
  # The human defaults (5000 / 25000) discard 12.7% of all cells and up to 21.2% of
  # H1D0_2 / H1D30_2 — real cells, not multiplets. 9000/70000 sits just above p99 for
  # every sample and drops the loss to 0.2%; scDblFinder (step 02) remains the actual
  # doublet filter, so these caps are only a backstop against extreme outliers (one
  # H1D30_2 cell reaches 435,190 UMI).
  QC$max_features <- 9000
  QC$max_counts   <- 70000
  # Do not regress %mt: cardiomyocytes are mitochondria-dense by function, so mito
  # load tracks the CM/progenitor axis here rather than damage (p99 is only 4.6%).
  NORM$vars_to_regress <- NULL
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
    # Epicardium vs proepicardium split 2026-10-24. These were one panel
    # (WT1/TBX18/ALDH1A2/UPK3B), which is really the PROEPICARDIAL transcription
    # factor set — so "Epicardial" was defined by progenitor markers and the two
    # states could not be separated. build_reference.R dropped Proepicardial
    # entirely for want of a panel, and the resulting frozen reference forced
    # those cells into Myofibroblast (28 pp drift on the held-out anchor).
    # Split on the biology: proepicardium is the TF-defined progenitor,
    # epicardium the differentiated mesothelial sheet. All genes verified present
    # in the cardiac features files.
    Epicardial    = c("UPK3B", "MSLN", "KRT19", "BNC1", "CCDC80"),
    Proepicardial = c("TBX18", "WT1", "SFRP5", "LHX2", "SEMA3D", "ALDH1A2"),
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
  # Cardiac-lineage and mesoderm-derived states of the differentiation. Off-target arms
  # (Hepatic/Endoderm, Amnion-like epithelium, Epithelial) are what the flag is for.
  EXPECTED_TYPES <- c("Pluripotent", "Cardiac progenitor", "Cardiomyocyte", "Epicardial",
                      "Proepicardial", "Fibroblast", "Myofibroblast", "Smooth Muscle",
                      "Pericyte", "Endothelial", "Proliferating",
                      "Mesenchyme")   # lateral plate mesoderm, the route to cardiac

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
  if (length(SAMPLE_NAMES) == 7 && isTRUE(QC$min_features == 200) && isTRUE(CLUSTER$algorithm == 1) &&
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
  if (length(SAMPLE_NAMES) == 7 && isTRUE(QC$min_features == 1000) && isTRUE(CLUSTER$algorithm == 1) &&
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

  # --- RETIRED 2026-10-24: SoupX-corrected cohort (*_sx), curated 2026-09-23 ---
  # Kept for the curation record only; the `FALSE &&` below makes it unreachable.
  # Retired because the cohort it describes can no longer be reproduced:
  #   * H1D30_1's source matrix was replaced by the 2026-10-22 re-run (rho 0.250 ->
  #     0.187), so H1D30_1_sx cannot be regenerated identically.
  #   * H1D0_2's rho here was a manual fallback (0.156) after autoEstCont failed;
  #     it estimates 0.266 automatically, so this cohort was under-corrected.
  # Superseded by the 2026-10-24 block below. Original note follows.
  #
  # Same 7 libraries after SoupX adjustCounts (rho 0.076-0.673; H1D20_1 lost 67.3% of
  # counts and 76.9% of its cells, confirming it was ambient-dominated). Correction
  # re-clusters the data, so the uncorrected map above MUST NOT be reused: different
  # sample names keep this block and that one mutually exclusive.
  # Curated from per-cluster detection rates (22 clusters, 61,697 cells).
  if (FALSE && length(SAMPLE_NAMES) == 7 && isTRUE(QC$min_features == 1000) &&
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

  # --- SoupX-corrected cohort with H1D11_1, H1D20_1 dropped, curated 2026-10-24 -
  # 7 libraries: H1D11_1 added (2026-10-22 re-run), H1D30_1 re-run matrix, and
  # H1D20_1 EXCLUDED (rho 0.795; lost 88.8% of cells to QC after correction, so
  # D20 is represented by its clean replicate H1D20_2_sx, rho 0.086, alone).
  # Curated from per-cluster detection rates (22 clusters, 76,314 cells, res 0.5).
  # Identity TFs decide each call; secreted cargo (AFP/ALB/TTR/APOA1, collagens)
  # does NOT, since it survives correction where real endoderm exists.
  # Unassignable fell from 26% (raw cohort) to 9.1% here.
  if (length(SAMPLE_NAMES) == 7 && isTRUE(QC$min_features == 1000) && isTRUE(CLUSTER$algorithm == 1) &&
      setequal(SAMPLE_NAMES, c("H1D0_1_sx", "H1D0_2_sx", "H1D11_1_sx", "H1D11_2_sx",
                               "H1D20_2_sx", "H1D30_1_sx", "H1D30_2_sx"))) {
    CLUSTER_CELLTYPE_MAP <- c(
      # --- pluripotent: POU5F1/SOX2/L1TD1, D0-dominant -------------------------
      "3"  = "Pluripotent",       # POU5F1 92 SOX2 83 LIN28A 87 L1TD1 86 NANOG 35; 98% D0
      "9"  = "Pluripotent",       # POU5F1 90 SOX2 63 L1TD1 67 LIN28A 81; 91% D0
      "11" = "Pluripotent",       # POU5F1 36 SOX2 32 L1TD1 71; 97% D0, low content (1,928 genes)
      # --- cardiomyocyte: the only cluster with sarcomere + NKX2-5 -------------
      "12" = "Cardiomyocyte",     # TTN 96 MYL7 93 MYL4 86 MYH6 81 ACTC1 84 TNNT2 69,
                                  # NKX2-5 56 MEF2C 49 TBX5 38 GATA4 66; 75% D11 20% D30
      # --- lateral-plate mesenchyme (renamed 2026-10-24, was "Cardiac progenitor"):
      # --- HAND1 89-92 MSX1 65 TMEM88 89 BMP4 70 with NKX2-5 3 ISL1 6 TBX5 5 — the
      # --- paper's "mesenchymal" population, precursor of the D20/D30 stroma -----
      "0"  = "Mesenchyme (lateral plate)", # HAND1 89 MSX1 65 TMEM88 89 LUM 85; NKX2-5 3; 91% D11
      "4"  = "Mesenchyme (lateral plate)", # HAND1 92 MSX1 67 TMEM88 93; NKX2-5 3; 91% D11
      "14" = "Mesenchyme (lateral plate)", # HAND1 46 GATA4 56 CDH11 58; low depth; 64% D11
      # --- cardiac progenitor (renamed 2026-10-24, was "Proepicardial"): TBX5+ GATA4+
      # --- first-heart-field progenitor = the lab's day-11 CCP. WT1 16 TBX18 11
      # --- ALDH1A2 2 rule out proepicardium; nearest clusters 0 (0.963) and 12 (0.941)
      "10" = "Cardiac progenitor", # GATA4 83 TBX5 37 NKX2-5 13 MYL7 82 SLC8A1 41 CCDC80 79; 90% D11
      # --- proepicardial: TBX18/WT1/ALDH1A2 — the only genuine epicardium --------
      "18" = "Proepicardial",     # TBX18 71 WT1 60 ALDH1A2 82 UPK3B 73 GATA4 84, NKX2-5 0;
                                  # 99% D30, 3,495 genes, n=215
      # --- epicardial-like mesothelium (renamed 2026-10-24): UPK3B/BNC1/KRT19/CDH1
      # --- but WT1 14-24 TBX18 1-2, HAND1+ TMEM88+ — epicardium-derived, not epicardium
      "1"  = "Epicardial-like",   # UPK3B 69 BNC1 63 KRT19 95 CDH1 64 WT1 14 TBX18 1; 68% D30
      "7"  = "Epicardial-like",   # UPK3B 41 BNC1 66 KRT19 80 WT1 24 TBX18 2; 68% D30
      # --- stromal -------------------------------------------------------------
      "17" = "Myofibroblast",     # ACTA2 89 TAGLN 91 PDGFRB 41 PDGFRA 30 TCF21 15
      "5"  = "Myofibroblast",     # ACTA2 70 TAGLN 72 PDGFRB 65 PDGFRA 59 TCF21 31; 86% D30
      # --- hepatic endoderm: definitive TFs, not cargo alone -------------------
      "6"  = "Hepatic/Endoderm",  # FOXA2 75 HHEX 55 ONECUT1 34 HNF4A 25 EPCAM 67; 93% D11
      "13" = "Hepatic/Endoderm",  # FOXA2 65 HNF4A 54 ONECUT1 16, TTR 93 APOA1 98 EPCAM 80
      "19" = "Hepatic/Endoderm",  # FOXA2 90 HHEX 64 ONECUT1 44; TOP2A 95 (cycling); 100% D11
      # --- other ---------------------------------------------------------------
      "15" = "Endothelial",       # CDH5 65 PECAM1 45 KDR 79 MEF2C 42; 89% D11
      # --- amnion-like epithelium (renamed 2026-10-24, was "Epithelial"): the
      # --- BMP4-driven off-target — GABRP/VTCN1/TFAP2A/GATA3/ISL1 ----------------
      "16" = "Amnion-like epithelium", # GABRP 54 VTCN1 29 TFAP2A 40 GATA3 35 ISL1 37 EPCAM 74
      "8"  = "Proliferating",     # TOP2A 93 MKI67 70; 59% D11
      # --- low-depth fibroblast lineage (renamed 2026-10-24, was "Unknown"):
      # --- DCN 93 LUM 96 COL3A1 99, nearest cluster 5 (Spearman 0.94). The hepatic
      # --- genes detected in 70-97% of cells carry a uniform 1.1% of UMIs per cell
      # --- (real endoderm, cluster 13: 3.4%) — ambient reads, not identity ------
      "2"  = "Fibroblast (low-depth)", # 1,117 genes (LOWEST); 81% D30. n=6,668
      "20" = "Unknown",           # n=184, 99% D30; MYL7 97 ACTC1 82 AND AFP 100 ALB 99
                                  # APOA1 99 DCN 99 — sarcomere and hepatic cargo together,
                                  # likely doublets
      "21" = "Unknown"            # n=83, 1,209 genes; AFP 83 ALB 73, no identity TF
    )

    # --- subtype curation, 2026-10-24 -------------------------------------
    # Three clusters hid substructure the lab needs: the cardiomyocytes, and two
    # clusters both labelled Myofibroblast. Curated from per-subcluster marker
    # ENRICHMENT against the run-wide baseline at res 0.4 on the saved RNA_snn
    # graph (raw detection is useless here — ACTA2 and COL1A1 are near-universal
    # in this cohort).
    #
    # NO ventricular/atrial split is made, deliberately. All six characterisable
    # subclusters of cluster 12 score CM immature 38-65% on MYH6 (+53 to +88
    # enrichment); ventricular appears in one subcluster on FHL2 alone. The
    # apparent "atrial 16-28%" is driven entirely by MYL4, which is an atrial
    # marker in adult heart but also a fetal/immature CM marker — at D11-D30 it
    # reports immaturity, not chamber identity. Chamber specification needs
    # retinoic acid patterning this protocol does not apply, and the source
    # paper (npj Regen Med s41536-023-00302-6) treats day-11 output as committed
    # progenitors that mature only after transplantation.
    CLUSTER_SUBCLUSTER_MAP <- list(
      # cardiomyocytes: immature throughout, with one cycling population
      list(cluster = "12", resolution = 0.4, labels = c(
        "0" = "Cardiomyocyte (immature)",   # n=935  CM pan 57 MYH6 +53; 67% D11
        "1" = "Cardiomyocyte (immature)",   # n=838  CM pan 59 MYH6 +61
        "2" = "Cardiomyocyte (immature)",   # n=657  CM pan 83 MYH6 +88; 92% D11
        "3" = "Cardiomyocyte (immature)",   # n=620  CM pan 71 MYH6 +82
        "4" = "Cardiomyocyte (immature)",   # n=515  CM pan 83 MYH6 +88
        "5" = "Proliferating",              # n=116  TOP2A +56, Proliferating 61
        "6" = "Cardiomyocyte (immature)")), # n=22   too small to characterise
      # cluster 5: a genuine fibroblast<->myofibroblast continuum, 82-90% D30
      list(cluster = "5", resolution = 0.4, labels = c(
        "0" = "Myofibroblast",              # Fibroblast 81 / Myofib 55, DCN +52
        "1" = "Myofibroblast",              # Fibroblast 75 / Myofib 61, TAGLN +31
        "2" = "Myofibroblast",              # Fibroblast 82 / Myofib 71, PDGFRA +52
        "3" = "Myofibroblast",              # Fibroblast 82 / Myofib 77, TAGLN +50
        "4" = "Fibroblast (resting)",       # n=676 TCF21 +37, Myofib only 43 — least activated
        "5" = "Fibroblast (activated)")),   # n=398 Myofib 86, TAGLN +58 — most activated
      # cluster 17: NOT one population. Smooth muscle and a cycling stromal
      # population were both collapsed into Myofibroblast at cluster resolution.
      list(cluster = "17", resolution = 0.4, labels = c(
        "0" = "Myofibroblast",              # n=414 Myofib 88 TAGLN +58; 85% D11
        "1" = "Proliferating",              # n=380 MKI67 +74, Proliferating 87; 88% D30
        "2" = "Smooth Muscle",              # n=281 SMC 64 on CNN1 +67; 81% D11
        "3" = "Smooth Muscle"))             # n=31  CNN1 +78; small, same identity as 17_2
    )
  }

  # --- Raw (uncorrected) 8-library cohort, curated 2026-10-24 -------------------
  # 20 clusters, res 0.5, min_features 1000. Until this date the run carried raw
  # SingleR labels (Neurons 72-93% of every differentiated sample) because every
  # earlier map was guarded on a 7-name set. Same evidence base as the _sx map
  # above (docs/cm_fresh_look_2026-10-24.md); cluster numbers differ, identities
  # correspond one-to-one. Cluster 0 (19,340 cells, 32% of the cohort) is where
  # uncorrected ambient RNA and H1D20_1 (rho 0.795) collect.
  if (length(SAMPLE_NAMES) == 8 && isTRUE(QC$min_features == 1000) && isTRUE(CLUSTER$algorithm == 1) &&
      setequal(SAMPLE_NAMES, c("H1D0_1", "H1D0_2", "H1D11_1", "H1D11_2",
                               "H1D20_1", "H1D20_2", "H1D30_1", "H1D30_2"))) {
    CLUSTER_CELLTYPE_MAP <- c(
      "4"  = "Pluripotent",                # POU5F1 97 SOX2 93 ESRG 97; 97% D0
      "8"  = "Pluripotent",                # POU5F1 74 XACT 90; low depth; 97% D0
      "11" = "Pluripotent",                # POU5F1 95 SOX2 74; 97% D0
      "1"  = "Mesenchyme (lateral plate)", # HAND1 87 MSX1 62 TMEM88 86; NKX2-5 5; 86% D11
      "3"  = "Mesenchyme (lateral plate)", # HAND1 74 MSX1 39 TMEM88 78; NKX2-5 6; low depth
      "10" = "Cardiac progenitor",         # GATA4 79 TBX5 31 NKX2-5 11 ZFPM2 69; WT1 17 TBX18 9; 81% D11
      "9"  = "Cardiomyocyte",              # MYH6 73 SLC8A1 76 TTN 92 TNNT2 64 NKX2-5 48 MYOCD 54
      "12" = "Proliferating",              # TOP2A 95 MKI67 74 with HAND1 87 BMP4 84
      "7"  = "Hepatic/Endoderm",           # FOXA2 75 HHEX 55 TBX3 55 HNF4A 25; 91% D11
      "13" = "Hepatic/Endoderm",           # APOB 85 MTTP 79 RBP4 88 FGB 88 HNF4A 53 FOXA2 63
      "14" = "Endothelial",                # CDH5 64 PECAM1 44 CD34 51 KDR 77; 85% D11
      "15" = "Amnion-like epithelium",     # GABRP 55 VTCN1 30 TFAP2A 40 GATA3 36 ISL1 38 EPCAM 74
      "2"  = "Epicardial-like",            # UPK3B 71 BNC1 66 CDH1 69 WT1 15 TBX18 2; D20/D30
      "5"  = "Epicardial-like",            # UPK3B 41 BNC1 61 WT1 20 TBX18 2; D20/D30
      "18" = "Proepicardial",              # WT1 59 TBX18 66 ALDH1A2 75 UPK3B 77 NPR3 69; 94% D30
      "6"  = "Myofibroblast",              # CDH11 88 FBN1 71 PDGFRB 67 ACTA2 85 TGFBI 92; 76% D30
      "16" = "Myofibroblast",              # ACTA2 95 TAGLN 97 CNN1 67 NPPB 43 — SMC-like/stressed
      "0"  = "Fibroblast (low-depth)",     # 1,220 genes; DCN 91 LUM 99; hepatic 0.9% of UMIs
                                           # per cell (ambient); 56% D20 incl. H1D20_1
      "17" = "Unknown",                    # n=372, percent.mt 6.6 (10x the cohort) — dying
      "19" = "Unknown"                     # n=207, 96% D30; ALB 97 AHSG 93 with ACTA2 97 — doublets
    )
  }

  # =============================================================================
  # Leiden maps (CLUSTER$algorithm == 4, pipeline of 2026-10-26: Leiden on the
  # igraph path, per-sample HVGs, no %mt regression). Cluster ids are 1-based.
  # Same Yap et al. 2023 vocabulary and the same marker rules as the Louvain maps
  # above (docs/cm_fresh_look_2026-10-24.md); evidence per cluster in
  # docs/cm_rerun_leiden_2026-10-26.md and /data/alvin/tmp/look_{7sx,8}.txt.
  # Numbers in comments are detection % in the cluster; "D11 88" = % of cells
  # from that timepoint.
  # Note: the old Louvain maps above are now also guarded on algorithm == 1.
  # =============================================================================

  # ---- 7-sample SoupX cohort, Leiden (results_H1D01sx_7samples_leiden) ----------
  if (length(SAMPLE_NAMES) == 7 && isTRUE(QC$min_features == 1000) && isTRUE(CLUSTER$algorithm == 4) &&
      setequal(SAMPLE_NAMES, c("H1D0_1_sx", "H1D0_2_sx", "H1D11_1_sx", "H1D11_2_sx",
                               "H1D20_2_sx", "H1D30_1_sx", "H1D30_2_sx"))) {
    CLUSTER_CELLTYPE_MAP <- c(
      "10" = "Pluripotent",                # ESRG 64 XACT 87 POU5F1 31; low depth (2,029 genes); D0 90
      "11" = "Pluripotent",                # POU5F1 88 ESRG 87 SOX2 61; D0 89
      "18" = "Pluripotent",                # POU5F1 91 SOX2 78 NANOG 35; D0 94
      "19" = "Pluripotent",                # POU5F1 97 ESRG 96 + CDC20/CCNB1 (cycling); D0 98
      "2"  = "Mesenchyme (lateral plate)", # mixed with CCP — split by CLUSTER_SUBCLUSTER_MAP below
      "12" = "Cardiomyocyte (immature)",   # MYH6 67 SLC8A1 82 TTN 91 TNNT2 56 MYOCD 59 NKX2-5 40 TBX5 39
      "3"  = "Proliferating",              # cycling cells of several lineages — split by lineage below
      "5"  = "Proliferating",              # cycling cells of several lineages — split by lineage below
      "17" = "Proliferating",              # TOP2A 73 MKI67 40 with HAND1 76 MSX1 50 (proliferating mesenchyme)
      "9"  = "Hepatic/Endoderm",           # FOXA2 59 HHEX 41 TBX3 47 HNF4A 27 TTR 59; D11 76
      "13" = "Hepatic/Endoderm",           # APOB 79 MTTP 79 RBP4 83 FGB 81 HNF4A 49 (yolk-sac-like); hepatic UMI 6%
      "14" = "Hepatic/Endoderm",           # FOXA2 58 HHEX 44 TTR 65 EPCAM 71; D11 69
      "15" = "Endothelial",                # CDH5 65 PECAM1 45 CD34 52 KDR 79; D11 89
      "16" = "Amnion-like epithelium",     # GABRP 52 VTCN1 28 TFAP2A 40 GATA3 34 ISL1 32 EPCAM 65
      "1"  = "Epicardial-like",            # UPK3B 45 BNC1 44 WT1 11 TBX18 1; D11 part split off below
      "8"  = "Epicardial-like",            # BNC1 55 UPK3B 34 CDH1 47 WT1 19 TBX18 1; D30 69
      "20" = "Proepicardial",              # WT1 50 TBX18 60 ALDH1A2 72 UPK3B 67 NPR3 71; D30 99 (true epicardium)
      "22" = "Proepicardial",              # WT1 72 TBX18 81 ALDH1A2 92 UPK3B 75 NPR3 69; D30 97
      "6"  = "Myofibroblast",              # TGFBI 71 CDH11 68 FBN1 49 PDGFRB 45 ACTA2 69 SNAI2/LOX/MMP2
      "7"  = "Myofibroblast",              # CDH11 75 PDGFRB 47 FBN1 48 ACTA2 55 TCF21 25; culture-activated fibroblast
      "4"  = "Fibroblast (low-depth)",     # 1,123 genes; DCN 81 LUM 85 with AFP 85 ALB 74 (ambient); D30 77
      "21" = "Unknown"                     # n=146; ALB 99 AHSG 91 AFP 100 with ACTA2 96 TAGLN 95 — doublets
    )
    # Cluster 2 (7,640, D11 88%) at res 0.15 on the SNN graph: 2_0 and 2_2 carry the
    # Yap CCP signature (GATA4 84/68, TBX5 36/19, NKX2-5 14/5, ZFPM2 68, SLC8A1 42/28,
    # WT1 19/12 — the paper's day-11 committed progenitor is GATA4 83 TBX5 37 NKX2-5
    # 13 WT1 16); 2_1 and 2_3 are HAND1 77/89 MSX1 54/63 TBX5 3/4 lateral plate.
    # Cell cycle is not a lineage (Yap: D0 is pluripotent, full stop). Clusters 3 and 5
    # are cycling cells drawn from D0 hESC, D11 mesenchyme/CCP and D30 mesothelium, so
    # they are split by lineage; "Proliferating" is reserved for HAND1+ cycling
    # mesenchyme as in the Louvain maps. Cluster 1 carries a D11 HAND1 88 MSX1 57
    # UPK3B 11 subcluster that is mesenchyme, not mesothelium.
    CLUSTER_SUBCLUSTER_MAP <- list(
      list(cluster = "2", resolution = 0.15,
           labels = c("0" = "Cardiac progenitor", "2" = "Cardiac progenitor",
                      "1" = "Mesenchyme (lateral plate)", "3" = "Mesenchyme (lateral plate)")),
      list(cluster = "1", resolution = 0.15,       # 1_0 UPK3B 60 BNC1 56 D30 58; 1_1 D11 73 HAND1 88; 1_2 D20/D30 DCN 91
           labels = c("0" = "Epicardial-like", "2" = "Epicardial-like", "1" = "Mesenchyme (lateral plate)")),
      list(cluster = "3", resolution = 0.30,       # 3_0 D0 76 POU5F1 76; 3_1/3_4 D11 HAND1 86/89 TBX5 1-2;
           labels = c("0" = "Pluripotent",         # 3_2 D11 97 GATA4 60 TBX5 15 NKX2-5 16 MYL7 89; 3_3 D30 66 UPK3B 56 BNC1 54
                      "1" = "Proliferating", "4" = "Proliferating",
                      "2" = "Cardiac progenitor", "3" = "Epicardial-like")),
      list(cluster = "5", resolution = 0.30,       # 5_0/5_3 D0 85/92 POU5F1 84/87; 5_1/5_5 D11 HAND1 86/97 TOP2A 94/98;
           labels = c("0" = "Pluripotent", "3" = "Pluripotent",   # 5_2 D0 34/D11 55 POU5F1 30 HAND1 47 GATA4 55 (exiting pluripotency, cycling);
                      "1" = "Proliferating", "5" = "Proliferating", "2" = "Proliferating",
                      "4" = "Cardiac progenitor"))   # 5_4 D11 91 GATA4 89 TBX5 45 NKX2-5 36 SLC8A1 62
    )
  }

  # ---- 8-sample raw cohort, Leiden (results_H1D01_8samples_leiden) ---------------
  # No SoupX, H1D20_1 (rho 0.795) included; clusters 1 and 4 are where ambient
  # collects (1,216 / 1,232 genes, AFP 78/87, DCN 77/86). No true-epicardium
  # (WT1+ TBX18+) cluster is resolved in this run; the closest cells sit in 14_0.
  if (length(SAMPLE_NAMES) == 8 && isTRUE(QC$min_features == 1000) && isTRUE(CLUSTER$algorithm == 4) &&
      setequal(SAMPLE_NAMES, c("H1D0_1", "H1D0_2", "H1D11_1", "H1D11_2",
                               "H1D20_1", "H1D20_2", "H1D30_1", "H1D30_2"))) {
    CLUSTER_CELLTYPE_MAP <- c(
      "6"  = "Pluripotent",                # cycling cells of several lineages — split by lineage below
      "5"  = "Pluripotent",                # mixed D0/D11 — split by CLUSTER_SUBCLUSTER_MAP below
      "13" = "Cardiac progenitor",         # GATA4 87 TBX5 37 NKX2-5 14 ZFPM2 73 WT1 22 TBX18 13; D11 79 (Yap CCP)
      "15" = "Cardiomyocyte",              # MYH6 63 SLC8A1 77 TTN 87 TNNT2 54 MYOCD 53 NKX2-5 38
      "3"  = "Proliferating",              # cycling cells of several lineages — split by lineage below
      "10" = "Proliferating",              # cycling cells of several lineages — split by lineage below
      "9"  = "Hepatic/Endoderm",           # FOXA2 62 HHEX 43 TBX3 49 TTR 67 HNF4A 29; D11 78
      "12" = "Hepatic/Endoderm",           # APOB 63 MTTP 62 RBP4 66 HNF4A 40 FOXA2 71 (yolk-sac-like); hepatic UMI 4.6%
      "18" = "Hepatic/Endoderm",           # FOXA2 74 HHEX 55 TBX3 47 TTR 84; n=531
      "16" = "Endothelial",                # CDH5 64 PECAM1 45 CD34 51 KDR 78; D11 85
      "17" = "Amnion-like epithelium",     # GABRP 46 TFAP2A 35 VTCN1 25 GATA3 31 ISL1 33 EPCAM 66
      "2"  = "Epicardial-like",            # UPK3B 50 BNC1 49 WT1 12; D0/D11 parts and true epicardium split off below
      "7"  = "Epicardial-like",            # BNC1 45 UPK3B 34 CDH1 43 WT1 15; DCN 89 LUM 94; D20/D30
      "14" = "Epicardial-like",            # mixed — split by CLUSTER_SUBCLUSTER_MAP below
      "8"  = "Myofibroblast",              # CDH11 69 FBN1 47 PDGFRB 47 ACTA2 79 TGFBI 74 SNAI2/LOX/MMP2; D30 48
      "11" = "Myofibroblast",              # CDH11 78 FBN1 50 PDGFRB 48 ACTA2 71 TCF21 30; D30 46
      "1"  = "Fibroblast (low-depth)",     # 1,216 genes; DCN 77 LUM 81 AFP 78 (ambient); shallow D0 hESC split off below
      "4"  = "Fibroblast (low-depth)",     # 1,232 genes, %mt 1.3; DCN 86 LUM 89 AFP 87, MT- genes as markers
      "19" = "Unknown"                     # n=222; ALB 96 AHSG 92 with ACTA2 95 — doublets
    )
    CLUSTER_SUBCLUSTER_MAP <- list(
      # Cluster 5 (7,275): 5_0 and 5_3 are D0 hESC (POU5F1 55/75, XACT 73/91);
      # 5_1 and 5_2 are D11 lateral plate (HAND1 87/90, MSX1 64/65, TBX5 3/6).
      list(cluster = "5", resolution = 0.15,
           labels = c("0" = "Pluripotent", "3" = "Pluripotent",
                      "1" = "Mesenchyme (lateral plate)", "2" = "Mesenchyme (lateral plate)")),
      # Cluster 14 (3,402): 14_0 D20/D30 BNC1 63 UPK3B 21 WT1 16 (epicardial-like);
      # 14_1 D11 81% GATA4 57 ZFPM2 64 TBX5 7 (early CCP, weaker than cl 13; the old
      # curated run called these cells Cardiac progenitor); 14_2 D30 CDH11 89 PDGFRB 60
      # ACTA2 68 (fibroblast); 14_3 GATA4 83 TBX5 38 SLC8A1 83 MYH6 53 TNNT2 36.
      list(cluster = "14", resolution = 0.15,
           labels = c("0" = "Epicardial-like", "1" = "Cardiac progenitor",
                      "2" = "Myofibroblast", "3" = "Cardiomyocyte")),
      # Cluster 2 (13,611): 2_0 D20/D30 UPK3B 62 BNC1 59; 2_1 D11 93 HAND1 96 MSX1 70;
      # 2_2 D0 83 POU5F1 83; 2_3 (293, D30 85) WT1 41 TBX18 38 ALDH1A2 44 UPK3B 66 —
      # the true-epicardium population this run otherwise lacks.
      list(cluster = "2", resolution = 0.15,
           labels = c("0" = "Epicardial-like", "1" = "Mesenchyme (lateral plate)",
                      "2" = "Pluripotent", "3" = "Proepicardial")),
      # Cluster 3 (8,972, cycling, res 0.30): 3_5/3_7 D0 87/77 POU5F1 87/77; 3_0 D0 49
      # D20 26 POU5F1 49 HAND1 39 (ambient-mixed, majority D0); 3_1/3_3 D11 HAND1 92;
      # 3_2 D30 61 UPK3B 64; 3_4/3_6 D11 96-98 GATA4 72/65 TBX5 25/12 NKX2-5 22/19.
      list(cluster = "3", resolution = 0.30,
           labels = c("5" = "Pluripotent", "7" = "Pluripotent", "0" = "Pluripotent",
                      "1" = "Proliferating", "3" = "Proliferating",
                      "2" = "Epicardial-like",
                      "4" = "Cardiac progenitor", "6" = "Cardiac progenitor")),
      # Cluster 10 (4,246, cycling): 10_1 D0 78 POU5F1 74; 10_0 D11 67 HAND1 74 TOP2A 80;
      # 10_2 D11 94 GATA4 93 TBX5 49 NKX2-5 23 SLC8A1 53; 10_3 D30 94 UPK3B 78 BNC1 97
      # WT1 48 but TBX18 3 (mesothelium, not epicardium proper).
      list(cluster = "10", resolution = 0.15,
           labels = c("1" = "Pluripotent", "0" = "Proliferating",
                      "2" = "Cardiac progenitor", "3" = "Epicardial-like")),
      # Cluster 1 (15,119, low depth): 1_1 (3,323) is D0 71% POU5F1 70 ESRG 70 — shallow
      # hESC that co-cluster with the ambient-heavy fibroblast bulk on depth alone.
      list(cluster = "1", resolution = 0.15,
           labels = c("0" = "Fibroblast (low-depth)", "1" = "Pluripotent")),
      # Cluster 6 (6,925, cycling, res 0.30): 6_1/6_2/6_4 D0 85-89 POU5F1 85; 6_0 D0 48
      # POU5F1 48 HAND1 38 (cycling, majority D0); 6_3 D11 67 HAND1 93 TOP2A 98;
      # 6_5 D11 71 GATA4 82 TBX5 31 NKX2-5 24 SLC8A1 48.
      list(cluster = "6", resolution = 0.30,
           labels = c("1" = "Pluripotent", "2" = "Pluripotent", "4" = "Pluripotent", "0" = "Pluripotent",
                      "3" = "Proliferating", "5" = "Cardiac progenitor"))
    )
  }

  # =============================================================================
  # ASF5 and CLMC lines (flowcells C and D, 2026-10-27), Leiden, same Yap rules.
  # Evidence: /data/alvin/tmp/look_{asf5,clmc}.txt, tp_{asf5,clmc}.txt, sub_{asf5,clmc}.txt.
  # Both cohorts carry heavy ambient RNA (ASF5: MYL7/TTR/AFP detected in 76-100 % of
  # cells in nearly every cluster), so labels rest on per-cluster SPECIFIC markers
  # (pct_in - pct_rest) and timepoint mix, not on raw detection. Clusters that span
  # D0-D30 are depth or cell-cycle groupings and are split by lineage below.
  # =============================================================================

  # ---- ASF5 (results_ASF5D01_8samples_filtered, 20 clusters) -------------------
  if (length(SAMPLE_NAMES) == 8 && isTRUE(QC$min_features == 1000) && isTRUE(CLUSTER$algorithm == 4) &&
      setequal(SAMPLE_NAMES, paste0("ASF5_D", rep(c(0, 11, 20, 30), each = 2), "_", 1:2))) {
    CLUSTER_CELLTYPE_MAP <- c(
      "1"  = "Unknown",                    # split below; 35,911 cells, 1,584 genes, all timepoints
      "2"  = "Hepatic/Endoderm",           # HNF4G FOXA1 ONECUT1 RFX6 LINC00261 FOXA2 80 HNF4A 67; D20 56
      "3"  = "Hepatic/Endoderm",           # split below; NR5A2 CUX2 ROR1 ADAMTS19
      "4"  = "Hepatic/Endoderm",           # FOXA2 74 CLDN3 RAB25 F2; D11 65
      "5"  = "Cardiomyocyte",              # SLC8A1 73 TTN 95 MYH6 62 TNNT2 51 MYOCD 44 TECRL ACTC1; cm UMI 1.3 %; D30 76
      "6"  = "Myofibroblast",              # CDH11 81 PDGFRB 66 PDGFRA 65 ZEB2 PRDM6 COL6A3; D30 63
      "7"  = "Proliferating",              # split below; HAND1+ cycling mesenchyme, D11-D30
      "8"  = "Hepatic/Endoderm",           # ONECUT1 SERPINA1 APOC3 CDH6; D20 68
      "9"  = "Cardiac progenitor",         # GATA4 86 ZFPM2 82 TBX5 30 CHRM2 WT1 25 NKX2-5 6; D30 60 [uncertain: late for a CCP]
      "10" = "Pluripotent",                # split below; POU5F1 MIR302CHG ESRG, low depth
      "11" = "Proepicardial",              # WT1 70 TBX18 72 ALDH1A2 57 NPR3 49 BNC1 76; D30 73 (true epicardium)
      "12" = "Epicardial-like",            # HAND1 87 BNC1 70 UPK3B 22 DCN 88 LUM 94 NPNT FRZB; D20/D30
      "13" = "Proliferating",              # split below; cycling D0 hESC + D11/D20 endoderm
      "14" = "Pluripotent",                # POU5F1 78 ESRG 63, low depth; D0 86
      "15" = "Endothelial",                # CDH5 67 PECAM1 58 CD34 53 EMCN CALCRL
      "16" = "Hepatic/Endoderm",           # split below; cycling, FOXA2 65 HHEX 59
      "17" = "Hepatic/Endoderm",           # FOXA2 HNF1A-AS1 VTN APOA2 HHEX
      "18" = "Pluripotent",                # POU5F1 91 SOX2 89 NANOG 84 UTF1 FOXD3-AS1; D0 91
      "19" = "Unknown",                    # n=1,116, 1,710 genes, all timepoints; RBFOX1 CSMD1 NRG3 long-gene signature
      "20" = "Epicardial-like"             # n=347, D30 99; ITLN1 MGP BNC1 46 WT1 31 TBX18 28 (mesothelium)
    )
    CLUSTER_SUBCLUSTER_MAP <- list(
      # 1_0 (28,643): 1,584 genes, D11-D30, every lineage marker at ambient level (cm UMI
      # share 0.14 % vs 1.27 % in true CMs) — no identity; 1_1 (7,268): D0 64 POU5F1 65.
      list(cluster = "1", resolution = 0.15, labels = c("0" = "Unknown", "1" = "Pluripotent")),
      # 3_2/3_4 D0 97-99 POU5F1 86-90; 3_0/3_1/3_3 GATA4 + FOXA2 39-73 HHEX endoderm;
      # 3_5 (444) GATA4 96 TBX5 21 SLC8A1 66 MYH6 47.
      list(cluster = "3", resolution = 0.15,
           labels = c("0" = "Hepatic/Endoderm", "1" = "Hepatic/Endoderm", "3" = "Hepatic/Endoderm",
                      "2" = "Pluripotent", "4" = "Pluripotent", "5" = "Cardiac progenitor")),
      # 7_1/7_2 TOP2A 72/99 HAND1 64/79 (cycling mesenchyme); 7_0/7_3 D20/D30 CDH11 61-71
      # PDGFRB DCN (activated fibroblast); 7_4 (253) GABRP 77 ISL1 32 (amnion).
      list(cluster = "7", resolution = 0.15,
           labels = c("1" = "Proliferating", "2" = "Proliferating", "0" = "Myofibroblast",
                      "3" = "Myofibroblast", "4" = "Amnion-like epithelium")),
      # 10_0/10_1 POU5F1 37/87 (D0 plus residual undifferentiated cells); 10_2 18 cells.
      list(cluster = "10", resolution = 0.15,
           labels = c("0" = "Pluripotent", "1" = "Pluripotent", "2" = "Unknown")),
      # 13_0 D0 56 POU5F1 57 (cycling hESC); 13_1 D11/D20 FOXA2 64 HHEX 48.
      list(cluster = "13", resolution = 0.15, labels = c("0" = "Pluripotent", "1" = "Hepatic/Endoderm")),
      # 16_0 D0 78 POU5F1 78; 16_1-16_3 FOXA2 82-94 HHEX 72-77 (cycling endoderm).
      list(cluster = "16", resolution = 0.15,
           labels = c("0" = "Pluripotent", "1" = "Hepatic/Endoderm", "2" = "Hepatic/Endoderm", "3" = "Hepatic/Endoderm"))
    )
  }

  # ---- CLMC (results_CLMCD01_8samples_filtered, 14 clusters) -------------------
  # CLMC goes mesothelial/fibroblast, with one small D11 cardiomyocyte population that only
  # a finer split of cluster 5 resolves (5_4 below; sarcomeric UMI share 1.55 % vs a
  # 0.2 % D11 ambient median). The D20 "TNNT2+" cells (14 %) are TNNT2 alone, no other
  # sarcomere gene, in HAND1+ mesenchyme: not cardiomyocytes. Evidence:
  # /data/alvin/tmp/clmc_cm_look.txt, clmc_cm_sub_louvain.txt (2026-11-03).
  if (length(SAMPLE_NAMES) == 8 && isTRUE(QC$min_features == 1000) && isTRUE(CLUSTER$algorithm == 4) &&
      setequal(SAMPLE_NAMES, paste0("CLMC_D", rep(c(0, 11, 20, 30), each = 2), "_", 1:2))) {
    CLUSTER_CELLTYPE_MAP <- c(
      "1"  = "Epicardial-like",            # BNC1 84 UPK3B 37 WT1 25 TBX18 2 AQP1 MGP NPNT DCN 100; D20/D30
      "2"  = "Unknown",                    # split below; 16,776 cells, 1,394 genes, D11 84
      "3"  = "Epicardial-like",            # TENM2 KCNQ5 SOX6 CNTN5 BNC1 59 DCN 93; D30 77 (as H1 Leiden cl 8)
      "4"  = "Pluripotent",                # POU5F1 86 SOX2 81 NANOG 72 UTF1; D0 86
      "5"  = "Myofibroblast",              # split below; TGFBI HAND2 CDH11 SNAI2 MFAP4
      "6"  = "Fibroblast (low-depth)",     # 1,141 genes; DCN 93 LUM 92 MGP FRZB + AFP (ambient); D30 69
      "7"  = "Epicardial-like",            # DCN MGP FRZB SOX6 TENM2 BNC1 66 UPK3B 18; D30 79
      "8"  = "Hepatic/Endoderm",           # FOXA2 63 CLDN3 F10 MTTP APOM VTN; D11 67
      "9"  = "Hepatic/Endoderm",           # FOXA2 92 HHEX 70 TBX3 82 HNF4A 39; D11 95
      "10" = "Proliferating",              # split below; cycling D0 hESC + D20 mesenchyme
      "11" = "Proliferating",              # split below; S-phase HAND1+ cycling mesenchyme
      "12" = "Hepatic/Endoderm",           # split below; ADAMTS19 ROR1 FGF13 (as ASF5 cl 3)
      "13" = "Endothelial",                # CDH5 52 PECAM1 43 CD34 56 HOPX; D11 89
      "14" = "Hepatic/Endoderm"            # FOXA2 96 HHEX 77 FOXA3 TBX3 87; D11 99
    )
    CLUSTER_SUBCLUSTER_MAP <- list(
      # 2_0/2_1 (16,093): D11, 1,345-1,540 genes; MYL7/TTN/MYH6 no higher than the D11
      # endoderm clusters (ambient), TBX5 5 NKX2-5 4 — no identity [possibly immature CMs,
      # needs SoupX]; 2_2 (683) D0 96 POU5F1 79.
      list(cluster = "2", resolution = 0.15,
           labels = c("0" = "Unknown", "1" = "Unknown", "2" = "Pluripotent")),
      # 5_0/5_1/5_3/5_5 D20/D30 HAND1 81-96 (myofibroblast); 5_2 (1,733) D11 44 GATA4 85
      # TBX5 19 ISL1 17 NKX2-5 11; 5_4 (581) D11 96 ACTC1 90 TTN 92 MYH6 78 TNNT2 68
      # NKX2-5 63 ISL1 67, 65 % of cells with >= 6 sarcomere genes (early cardiomyocytes).
      list(cluster = "5", resolution = 0.30,
           labels = c("0" = "Myofibroblast", "1" = "Myofibroblast", "3" = "Myofibroblast",
                      "5" = "Myofibroblast", "2" = "Cardiac progenitor", "4" = "Cardiomyocyte")),
      # 10_0 D0 86 POU5F1 86 (cycling hESC); 10_1/10_2 HAND1 89-92 TOP2A 92-98.
      list(cluster = "10", resolution = 0.15,
           labels = c("0" = "Pluripotent", "1" = "Proliferating", "2" = "Proliferating")),
      # 11_0/11_1/11_3 HAND1 95-98 TOP2A 69-98; 11_2 (825) D11 86 GATA4 93 TBX5 31
      # NKX2-5 21 ISL1 31 SLC8A1 47 (Yap CCP).
      list(cluster = "11", resolution = 0.15,
           labels = c("0" = "Proliferating", "1" = "Proliferating", "3" = "Proliferating",
                      "2" = "Cardiac progenitor")),
      # 12_0 D11 FOXA2 72; 12_1 D0 93 POU5F1 79; 12_2 D11 FOXA2 56 GATA4 91 TBX5 9
      # (below the CCP TBX5 bar, kept with endoderm).
      list(cluster = "12", resolution = 0.15,
           labels = c("0" = "Hepatic/Endoderm", "2" = "Hepatic/Endoderm", "1" = "Pluripotent"))
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
