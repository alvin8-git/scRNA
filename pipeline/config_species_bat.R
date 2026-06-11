# =============================================================================
# config_species_bat.R - bat / bat_wing marker, QC, and reference overrides.
# Sourced by config.R AFTER the base (human) definitions; mutates MARKERS, QC,
# CLUSTER, SINGLER_REF, SUBTYPE_MARKERS, CONTAMINATION_TYPES, CELLTYPE_COLORS,
# WOUND_MODULES, etc. Expects `.species` to be set by config.R. Human = no-op.
# =============================================================================
if (.species == "bat") {
  message("[Species] bat (Eonycteris spelaea) — applying whole-blood overrides")
  # Ground truth: Gamage et al. 2022 (Immunity 55, 2187-2205) — E. spelaea scRNA-seq.
  # All bat-validated markers below unless explicitly noted as ortholog inference.
  # Gene names follow human HGNC nomenclature via ESpe GTF orthology (Dec 2024).

  # ---- T cells ---------------------------------------------------------------
  # CD8_T: remove GZMK (21% bat cells; shared NK/CD4/CD8 — too broad)
  # PRF1 validated for E. spelaea effector NK & T (Gamage 2022 Fig 2C)
  MARKERS$CD8_T <- c("CD8A", "CD8B", "PRF1")

  # ---- NK cells --------------------------------------------------------------
  # KLRB1 (CD161) replaces KLRD1: KLRB1 defines CD8+CD161+ NK&T in E. spelaea
  # (Gamage 2022 Fig 2C/3C); NKG2 family (KLRD1) diverged in Chiroptera
  MARKERS$NK <- c("NKG7", "GNLY", "KLRB1")

  # ---- B cells ---------------------------------------------------------------
  # CD79B and FCMR (IgM Fc receptor) directly validated in E. spelaea
  # (Gamage 2022 Fig 2B: "high expression of canonical markers CD79B, MS4A1, FCMR")
  MARKERS$B_cell <- c("CD79A", "CD79B", "MS4A1", "FCMR")

  # ---- DC: bat-validated markers spanning cDC1 and pDC ----------------------
  # CLEC9A = cDC1 lineage marker; TCF4/IRF7 = pDC master TFs — all PDF-validated
  MARKERS$DC <- c("CLEC9A", "TCF4", "IRF7")

  # ---- Classical monocytes ---------------------------------------------------
  # S100A12 + SIRPA (CD172α) validated against non-classical mono (Gamage 2022 Fig 4F)
  MARKERS$CD14_mono <- c("CD14", "LYZ", "S100A12", "SIRPA", "CSF1R")

  # ---- Non-classical monocytes -----------------------------------------------
  # CX3CR1 replaces FCGR2A: CX3CR1 is the definitive non-classical mono marker in
  # E. spelaea (Gamage 2022 Fig 2C); FCGR3A included (expressed in bat NC mono, Fig 4F)
  MARKERS$FCGR3A_mono <- c("CX3CR1", "FCGR3A", "CDKN1C", "MS4A7")

  # ---- Neutrophils -----------------------------------------------------------
  # IDO1, ALAS1, SLC16A10: >30-fold enriched in E. spelaea neutrophils vs all
  # other immune cell types; validated by bulk RNA-seq on FACS-sorted cells
  # (Gamage 2022 Fig 6B/C). Tryptophan-degradation axis is bat-specific.
  # FCGR3B/CEACAM6 removed: FCGR and CEACAM gene families divergent in Chiroptera
  # MPO/ELANE: canonical azurophilic granule proteins; conserved across mammals (Reviewer 3)
  MARKERS$Neutrophil <- c("CSF3R", "CXCR2", "IDO1", "ALAS1", "SLC16A10", "MPO", "ELANE")

  # ---- HSPC ------------------------------------------------------------------
  # Remove AVP (arginine vasopressin — neurohypophyseal hormone, wrong context in PBMC)
  MARKERS$HSPC <- c("CD34", "GATA2")

  # ---- γδ T ------------------------------------------------------------------
  MARKERS$gamma_delta_T <- c("TRDC", "TRGC1", "TRGC2")

  # ---- Eosinophil: drop SIGLEC8 (bat orthology unconfirmed; Reviewer 3) ------
  MARKERS$Eosinophil <- c("CCR3", "EPX")

  ALL_MARKERS <- unique(unlist(MARKERS))

  # ---- Contamination: RBC and Neutrophil expected in bat whole blood ---------
  CONTAMINATION_TYPES <- c("Basophil", "Eosinophil", "Mast cell")

  # ---- SingleR: Monaco resolves CD4/CD8/γδ T cells in blood ------------------
  SINGLER_REF <- "MonacoImmune"

  # ---- Clustering: higher resolution for whole-blood diversity ----------------
  CLUSTER$resolutions <- c(0.3, 0.5, 0.8, 1.0)
  CLUSTER$default_res <- 1.0
  CLUSTER$compare_res <- c(0.5, 0.8, 1.0)

  # ============================================================================
  # Sub-type markers — bat ground truth (Gamage 2022 + ortholog confidence)
  # ============================================================================

  # CD4 T subtypes
  # naive:    SELL removed (50.4% bat cells; non-discriminating in E. spelaea)
  # effector: GZMK removed (too broad); GZMB/TNFRSF4/PRF1 retained
  # memory:   AQP3 removed (no bat validation); S100A4 removed (monocyte cross-expression)
  SUBTYPE_MARKERS[["CD4 T"]] <- list(
    "CD4 T (naive)"    = c("CCR7", "TCF7", "LEF1"),
    "CD4 T (effector)" = c("GZMB", "TNFRSF4", "PRF1"),
    "CD4 T (memory)"   = c("IL7R", "GPR183", "CD27")
  )

  # CD8 T subtypes
  # XCL1: selectively marks effector CD8+ DPP4+ NK&T in E. spelaea,
  # strongly induced upon viral infection (Gamage 2022 Fig 7I)
  # CTSW: cathepsin W; CD8-lineage-specific cysteine protease; conserved (Reviewer 3)
  SUBTYPE_MARKERS[["CD8 T"]] <- list(
    "CD8 T (effector)" = c("PRF1", "XCL1", "GZMB", "CTSW"),
    "CD8 T (memory)"   = c("IL7R", "GPR183", "CD27"),
    "CD8 T (naive)"    = c("TCF7", "LEF1", "CCR7")
  )

  # B cell subtypes — isotype genes (IGHD, IGHM, IGHG1) absent from bat GTF
  # IL4R dropped: shared with T cells, not B-specific (reviewer consensus)
  # AIM2 dropped: inflammasome sensor; not a core B cell marker
  SUBTYPE_MARKERS[["B cell"]] <- list(
    "B cell (naive)"  = c("TCL1A", "CD24", "FCER2"),
    "B cell (memory)" = c("CD27", "TNFRSF13B"),
    "Plasma"          = c("MZB1", "JCHAIN", "SDC1", "CD38", "XBP1", "PRDM1")
  )

  # Monocyte subtypes — S100A12/SIRPA validated for classical (Gamage 2022 Fig 4F)
  # FCN1: ficolin 1; validated classical monocyte marker (Reviewer 3)
  # LST1: leukocyte-specific transcript 1; marks non-classical monocytes (Reviewer 3)
  SUBTYPE_MARKERS[["Monocyte"]] <- list(
    "CD14+ Mono"   = c("CD14", "S100A12", "SIRPA", "LYZ", "FCN1"),
    "FCGR3A+ Mono" = c("CX3CR1", "FCGR3A", "CDKN1C", "MS4A7", "LST1")
  )

  # DC subtypes — Gamage 2022 Fig 4E/F resolves 3 DC populations in E. spelaea
  # cDC1: CADM1/CLEC9A/BATF3/IRF8 — validated in bat mDC sub-clustering
  #        XCR1: canonical cross-presenting cDC1 marker; conserved (Reviewer 3)
  # cDC2: FCER1A/CD74 — CD14low, high antigen-presenting genes (Gamage 2022 Fig 4F/text)
  #        NOTE: FCGR2B and S100A12 are explicitly ABSENT in bat cDC2 (p.2193) — do not use
  # pDC:  TCF4 (E2-2)/IRF7/IRF8 — constitutive IFN producers; bat immune hallmark
  SUBTYPE_MARKERS[["DC"]] <- list(
    "cDC1" = c("CADM1", "CLEC9A", "BATF3", "IRF8", "XCR1"),
    "cDC2" = c("FCER1A", "CD74"),
    "pDC"  = c("TCF4", "IRF7", "IRF8")
  )
} else if (.species == "bat_wing") {
  message("[Species] bat_wing (Eonycteris spelaea wing tissue) — applying wound-tissue overrides")

  # ---- Reference: broad tissue atlas, not blood-optimised -------------------
  SINGLER_REF <- "HumanPrimaryCellAtlas"

  # ---- Wing tissue markers --------------------------------------------------
  MARKERS$Fibroblast         <- c("COL1A1", "COL1A2", "COL3A1", "VIM", "PDGFRA", "FAP")
  MARKERS$Myofibroblast      <- c("ACTA2", "TAGLN", "MYL9", "CNN1")
  MARKERS$Keratinocyte       <- c("KRT5", "KRT14", "KRT1", "KRT10", "EPCAM")
  MARKERS$Wound_keratinocyte <- c("KRT6A", "KRT16", "KRT17", "MMP1")
  MARKERS$Endothelial        <- c("PECAM1", "CDH5", "VWF", "KDR", "FLT1")
  MARKERS$Pericyte           <- c("PDGFRB", "RGS5", "CSPG4", "NOTCH3")
  MARKERS$Macrophage         <- c("CD68", "CSF1R", "MRC1", "CD163")
  MARKERS$Melanocyte         <- c("MLANA", "DCT", "TYRP1", "MITF")
  MARKERS$FCGR3A_mono        <- NULL
  MARKERS$Neutrophil         <- c("S100A8", "S100A9", "FCGR3B", "CSF3R")
  ALL_MARKERS <- unique(unlist(MARKERS[!sapply(MARKERS, is.null)]))

  # ---- No blood contamination types in wing tissue --------------------------
  CONTAMINATION_TYPES <- c("RBC", "HSPC", "Platelet")

  # ---- Clustering: moderate resolution for tissue heterogeneity -------------
  CLUSTER$resolutions <- c(0.3, 0.5, 0.8)
  CLUSTER$default_res <- 0.5
  CLUSTER$compare_res <- c(0.3, 0.5, 0.8)

  # ---- γδ T markers ---------------------------------------------------------
  MARKERS$gamma_delta_T <- c("TRDC", "TRGC1", "TRGC2")
  ALL_MARKERS <- unique(unlist(MARKERS[!sapply(MARKERS, is.null)]))

  # ---- Sub-type markers for wing tissue -------------------------------------
  SUBTYPE_MARKERS[["Fibroblast"]] <- list(
    "Fibroblast (resting)" = c("PDGFRA", "DCN", "LUM", "CFD"),
    "Myofibroblast"        = c("ACTA2", "TAGLN", "MYL9", "POSTN"),
    "Fibroblast (wound)"   = c("COL3A1", "FN1", "SPARC", "CTGF")
  )
  SUBTYPE_MARKERS[["Macrophage"]] <- list(
    "Macrophage (M1/inflam)"  = c("IL1B", "TNF", "CXCL8", "CCL3"),
    "Macrophage (M2/repair)"  = c("MRC1", "CD163", "ARG1", "TGFB1"),
    "Macrophage (proliferat)" = c("MKI67", "TOP2A", "CDK1")
  )
  SUBTYPE_MARKERS[["Keratinocyte"]] <- list(
    "Keratinocyte (basal)"      = c("KRT5", "KRT14", "TP63", "COL17A1"),
    "Keratinocyte (suprabasal)" = c("KRT1", "KRT10", "LOR", "FLG"),
    "Keratinocyte (wound)"      = c("KRT6A", "KRT16", "MMP1", "LAMC2")
  )

  # ---- Wound-healing gene module sets (used by 11_wing_degs.R) --------------
  WOUND_MODULES <- list(
    Inflammatory     = c("IL1B", "TNF", "CXCL8", "CCL2", "IL6", "S100A8", "S100A9", "PTGS2"),
    ECM_remodeling   = c("COL1A1", "COL3A1", "MMP1", "MMP2", "MMP9", "TIMP1", "TIMP2", "LOXL2", "POSTN"),
    Angiogenesis     = c("VEGFA", "VEGFB", "FGF2", "ANGPT1", "ANGPT2", "KDR", "NRP1"),
    Proliferation    = c("MKI67", "TOP2A", "PCNA", "CDK1", "CCNB1", "CCNA2"),
    Re_epithelialize = c("KRT6A", "KRT16", "MMP1", "MMP3", "LAMC2", "ITGA3", "ITGB4"),
    Myofibroblast    = c("ACTA2", "MYL9", "TAGLN", "CNN1", "POSTN", "CTGF")
  )
}
