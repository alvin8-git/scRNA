# =============================================================================
# config_species_mouse.R — mouse overrides (first cohort turned out to be immune cells; see reference note)
#
# Sourced by config.R AFTER its own definitions, so everything here overwrites
# the human defaults in place. Self-guarded on .species, so it is inert for
# every other species.
#
# Mouse is NOT a validated sample type (see docs/report_product_design.md §0).
# The panels below are verified gene-by-gene against the real features file of
# E250164106_POS_A/0hCell_1 — all 74 symbols present out of 32,533 — but nobody
# has yet checked whether the resulting cell-type calls are biologically right.
# Treat output as unreviewed until the replicate-consistency check passes.
#
# THE critical override is the gene-name patterns. Mouse symbols are
# lowercase/title-case (mt-Nd1, Rps2), so the human "^MT-" matches NOTHING: 13
# mitochondrial genes match "^mt-" and 0 match "^MT-" in this data. Without the
# override percent.mt is 0 for every cell, QC$max_percent_mt silently does
# nothing, and the QC violin renders a flat line that looks like a clean
# dataset. 01_load_qc.R now warns when a pattern matches no genes.
# =============================================================================

if (.species %in% c("mouse", "mouse_cm")) {

  message("[Species] ", .species, " — applying mouse overrides (ImmGen reference)")

  # ---- gene-name patterns: the blocker ---------------------------------------
  QC$mt_pattern   <- "^mt-"      # 13 genes in this annotation; "^MT-" matches 0
  QC$ribo_pattern <- "^Rp[sl]"   # 101 genes

  # ---- QC bounds: PROVISIONAL ------------------------------------------------
  # Deliberately permissive. Cardiomyocytes are mitochondria-dense by function,
  # so a cutoff tuned for PBMCs discards healthy cells — and until percent.mt
  # actually computes there is no distribution to choose from. Over-filtering is
  # unrecoverable without reprocessing; under-filtering can be tightened later.
  # ACTION: run step 01 on two samples, read the percent.mt and nFeature
  # distributions, then set these from the data. Changing them busts the 01-03
  # sample cache, so expect one reprocess.
  QC$min_features   <- 500
  QC$max_features   <- 9000
  QC$min_counts     <- 1000
  QC$max_counts     <- 70000
  QC$max_percent_mt <- 30        # provisional — see note above

  # ---- reference -------------------------------------------------------------
  # celldex::MouseRNAseqData was fetched and inspected: its label.main includes
  # "Cardiomyocytes", which is exactly the class HumanPrimaryCellAtlas lacks and
  # the reason human cardiac cells get called "Neurons" by de-novo SingleR.
  # Gene symbols in it are title-case (Xkr4, Sox17), matching this data.
  #
  # REVISED 2026-10-25 after the first PosA run (results_0hCell1_8samples):
  # the cells are NOT cardiomyocytes. Ptprc is detected in 80-98 % of cells in
  # every cluster, Lyz2 in 83-100 %, Tnnt2 in <= 7 %, and the unbiased cluster
  # markers are Adgre1/Mertk/Msr1 (macrophage), Cxcr2/Hdc/Mmp8 (neutrophil),
  # Ciita/Flt3/Clec9a (DC), Cd79a/Pax5 (B), Itk/Cd247 (T), Nkg7/Gzma (NK),
  # Flt1/Ptprb (endothelial, 171 cells). This is a mouse immune-cell dataset,
  # so the immune-only ImmGen reference is the right one and MouseRNAseqData's
  # single "Macrophages"/"Monocytes" classes are too coarse.
  SINGLER_REF <- "ImmGen"

  # ---- marker panels ---------------------------------------------------------
  # Panel NAMES match CELLTYPE_COLORS entries and the cell-type labels, which
  # matters twice over: scType scores panels by name, and build_reference.R
  # resolves a label's defining markers by direct MARKERS[[label]] lookup.
  # "Macrophage" and "RBC" are used rather than "Immune"/"Erythrocyte" because
  # those two have no palette entry.
  # Every symbol below was confirmed present in the data.
  # Immune panels (2026-11-06): the only mouse cohort so far is immune cells (see the
  # reference note), and the cardiac panels below scored it against Pluripotent /
  # Cardiomyocyte. Named after the labels the ImmGen run produced. The cardiac panels
  # stay for the first mouse cardiac cohort (`mouse_cm`).
  MARKERS <- list(
    `Macrophage`   = c("Adgre1", "C1qa", "C1qb", "Mertk", "Msr1", "Cd68", "Lyz2"),
    `CD14+ Mono`   = c("Ly6c2", "Ccr2", "Plac8", "Cd14"),
    `Neutrophil`   = c("S100a8", "S100a9", "Retnlg", "Cxcr2", "Csf3r", "Mmp8", "Ly6g"),
    `DC`           = c("Flt3", "Clec9a", "Xcr1", "Cd209a", "H2-Ab1", "Ciita", "Itgax"),
    `B cell`       = c("Cd79a", "Cd79b", "Ms4a1", "Cd19", "Pax5"),
    `CD4 T`        = c("Cd3e", "Cd3g", "Cd4", "Il7r", "Itk"),
    `CD8 T`        = c("Cd3e", "Cd8a", "Cd8b1"),
    `NK`           = c("Nkg7", "Gzma", "Klrb1c", "Ncr1"),
    `Endothelial`  = c("Pecam1", "Cdh5", "Flt1", "Ptprb"),
    `Proliferating` = c("Mki67", "Top2a", "Ccnb1")
  )
  if (.species == "mouse_cm") MARKERS <- list(
    `Pluripotent`        = c("Pou5f1", "Sox2", "Nanog", "Lin28a", "L1td1", "Zfp42"),
    `Cardiac progenitor` = c("Nkx2-5", "Isl1", "Tbx5", "Gata4", "Mesp1", "Kdr", "Mef2c"),
    `Cardiomyocyte`      = c("Tnnt2", "Tnni3", "Tnni1", "Myh6", "Myh7", "Myl7", "Myl4",
                             "Myl2", "Actc1", "Ttn", "Actn2", "Ryr2", "Pln", "Nppa",
                             "Mybpc3", "Des"),
    `Epicardial`         = c("Upk3b", "Wt1", "Tbx18", "Aldh1a2", "Krt19", "Msln"),
    `Fibroblast`         = c("Col1a1", "Col3a1", "Dcn", "Lum", "Postn", "Pdgfra",
                             "Tcf21", "Fn1"),
    `Myofibroblast`      = c("Acta2", "Tagln", "Pdgfrb", "Cnn1"),
    `Smooth Muscle`      = c("Acta2", "Tagln", "Myh11", "Cnn1"),
    `Pericyte`           = c("Rgs5", "Abcc9", "Kcnj8", "Notch3"),
    `Endothelial`        = c("Pecam1", "Cdh5", "Tie1", "Emcn"),
    `Macrophage`         = c("Ptprc", "Cd68", "Adgre1", "Lyz2", "Itgam"),
    `Proliferating`      = c("Mki67", "Top2a", "Ccnb1"),
    # No RBC panel. With it, scType called five immune clusters (12,266 cells,
    # Ptprc 54-98 %) "RBC" on a haemoglobin UMI share of 0.00-0.22 %, and the
    # Monaco-blind override then overwrote the correct ImmGen labels. There are
    # no erythrocytes in this data to protect.
    # Off-target endoderm. Worth carrying because the human cardiac cohort
    # produced 9.8% hepatic/endoderm — a known failure mode of Wnt-modulation
    # protocols — and the same arm should be visible here if it occurs.
    `Hepatic/Endoderm`   = c("Afp", "Alb", "Ttr", "Apoa1", "Foxa2", "Hnf4a")
  )
  ALL_MARKERS <- unique(unlist(MARKERS, use.names = FALSE))

  # ---- contamination / rare types -------------------------------------------
  # Types the annotation step may override on marker evidence rather than trust
  # the reference for. None for mouse: see the RBC note in the panel block.
  CONTAMINATION_TYPES <- character(0)
  EXPECTED_TYPES <- if (.species == "mouse_cm") NULL else
    setdiff(names(MARKERS), "Endothelial")   # immune cohort: endothelial cells are carry-over

  # ---- clustering ------------------------------------------------------------
  # Same resolutions as the cm overlay; a differentiation timecourse needs the
  # finer options available for curating intermediate states.
  if (!is.null(CLUSTER$resolutions)) CLUSTER$resolutions <- c(0.1, 0.3, 0.5, 0.8)

  # ---- things deliberately NOT set ------------------------------------------
  # CLUSTER_CELLTYPE_MAP stays NULL. Cluster numbers are not stable across runs
  # and nothing has been curated for mouse yet; a map here would be a guess.
  # SUBTYPE_MARKERS / WOUND_MODULES are bat-wing concepts and stay unset.
  if (exists("WOUND_MODULES")) WOUND_MODULES <- NULL

  # ---- first analysis: replicate consistency --------------------------------
  # The design is a timecourse (0h, 6h, 24h, 72h) x arm (Cell, Med) x duplicate,
  # plus MS_Sham1/2 on PosA only. Flowcell is CONFOUNDED with timepoint — 0h and
  # 6h on PosA, 24h and 72h on PosB — so any cross-timepoint difference is
  # inseparable from a batch effect.
  # Replicate consistency is the one question that confound does NOT affect,
  # because replicates share a flowcell. That is why it is the first analysis.
  # Caveats already visible in the raw data: 72hMed_a is 3 MB against 30-360 MB
  # for every other sample (likely a failed library), and its _a/_b naming
  # breaks the _1/_2 replicate-stripping used to derive condition groups.
}
