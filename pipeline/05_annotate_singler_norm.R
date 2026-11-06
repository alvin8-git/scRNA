# =============================================================================
# 05_annotate_singler_norm.R - SINGLER_NORM lookup mapping raw SingleR labels
# (HumanPrimaryCellAtlas + MonacoImmune) to canonical CELLTYPE_COLORS names.
# Sourced by 05_annotate.R. Pure data; no side effects.
# =============================================================================
SINGLER_NORM <- c(
  # --- HPCA labels relevant to iPSC-CM differentiation (SCRNA_SPECIES=cm) ---
  # HPCA has NO cardiomyocyte label, so nothing here can produce one; cardiomyocyte
  # identity comes from the marker route. These entries only keep the labels HPCA
  # CAN emit from splitting into near-duplicate names in the palette.
  "Embryonic_stem_cells"   = "Pluripotent",
  "iPS_cells"              = "Pluripotent",
  "Hepatocytes"            = "Hepatic/Endoderm",
  "Neuroepithelial_cell"   = "Epithelial",

  # --- HumanPrimaryCellAtlas labels ---
  # Platelet / megakaryocyte
  "Platelets"              = "Platelet",
  "Megakaryocyte"          = "Platelet",
  # Neutrophil / granulocyte precursors
  "Neutrophils"            = "Neutrophil",
  "Neutrophil_-G-CSF"      = "Neutrophil",
  "GMP"                    = "Neutrophil",
  "Pro-Myelocyte"          = "Neutrophil",
  "Myelocyte"              = "Neutrophil",
  # Basophil / eosinophil / mast
  "Basophils"              = "Basophil",
  "Eosinophils"            = "Eosinophil",
  "Mast_cells"             = "Mast cell",
  # RBC / erythroid
  "Erythrocyte"            = "RBC",
  "Erythroblast"           = "RBC",
  "BFU-E"                  = "RBC",
  "CFU-E"                  = "RBC",
  "MEP"                    = "RBC",
  # HSPC
  "CMP"                    = "HSPC",
  "HSC_-G-CSF"             = "HSPC",
  "HSC_CD34+"              = "HSPC",
  # Lymphoid
  "NK_cell"                = "NK",
  "T_cells"                = "CD4 T",
  "B_cell"                 = "B cell",
  "Pro-B_cell_CD34+"       = "B cell",
  "Pre-B_cell_CD34-"       = "B cell",
  # Stromal / tissue contamination
  "Endothelial_cells"      = "Endothelial",
  "Epithelial_cells"       = "Epithelial",
  "Fibroblasts"            = "Fibroblast",
  "Smooth_muscle_cells"    = "Smooth Muscle",
  "Keratinocytes"          = "Keratinocyte",
  "Chondrocytes"           = "Chondrocyte",
  "MSC"                    = "MSC",
  "Tissue_stem_cells"      = "MSC",
  "Macrophage"             = "Macrophage",
  "Macrophages"            = "Macrophage",
  "Osteoblasts"            = "MSC",
  # --- MouseRNAseqData labels (SCRNA_SPECIES=mouse) ---
  # Space-separated, unlike HPCA's underscores. "Cardiomyocytes" is the one that
  # matters: without it the raw plural passes through and never meets the
  # "Cardiomyocyte" marker panel or palette entry. Neurons / Astrocytes /
  # Oligodendrocytes / Adipocytes are left raw so an off-target call stays visible.
  "Cardiomyocytes"         = "Cardiomyocyte",
  "Endothelial cells"      = "Endothelial",
  "Epithelial cells"       = "Epithelial",
  "Erythrocytes"           = "RBC",
  "Granulocytes"           = "Neutrophil",
  "Monocytes"              = "CD14+ Mono",
  "Dendritic cells"        = "DC",
  "T cells"                = "CD4 T",
  "Microglia"              = "Macrophage",
  # --- ImmGenData label.main (SCRNA_SPECIES=mouse, immune datasets) ---
  # "Monocytes", "Neutrophils", "Macrophages", "Microglia", "Endothelial cells",
  # "Epithelial cells", "Fibroblasts", "B cells", "NK cells", "T cells", "DC"
  # are already covered above. NKT and ILC are left raw: neither has a palette
  # entry and folding them into T / NK would hide a real population.
  "B cells, pro"           = "B cell",
  "Basophils"              = "Basophil",
  "Eosinophils"            = "Eosinophil",
  "Mast cells"             = "Mast cell",
  "Tgd"                    = "γδ T",
  "Stem cells"             = "HSPC",
  "Stromal cells"          = "MSC",
  # --- MonacoImmuneData labels ---
  "CD4+ T cells"           = "CD4 T",
  "CD8+ T cells"           = "CD8 T",
  "T regulatory cells"     = "Treg",
  "Vd2 gd T cells"         = "γδ T",
  "Non Vd2 gd T cells"     = "γδ T",
  "MAIT cells"             = "CD8 T",
  "NK cells"               = "NK",
  "B cells"                = "B cell",
  "Plasmablasts"           = "Plasma",
  "Classical monocytes"    = "CD14+ Mono",
  "Intermediate monocytes" = "CD14+ Mono",
  "Non-classical monocytes" = "FCGR3A+ Mono",
  "Plasmacytoid DC"        = "DC",
  "Myeloid DC"             = "DC",
  "Progenitor cells"       = "HSPC",
  "Low-density neutrophils" = "Neutrophil",
  "Low-density basophils"  = "Basophil",
  # Monaco label.main catch-alls (coarser labels returned for ambiguous cells)
  "Monocytes"              = "CD14+ Mono",
  "T cells"                = "CD4 T",
  "Dendritic cells"        = "DC",
  "B cells"                = "B cell",
  "NK"                     = "NK"
)

# -----------------------------------------------------------------------------
# split_generic_t(): optional replacement for the generic "T cells" -> "CD4 T"
# collapse above (ANNOTATION$split_generic_t, default FALSE). Cells whose RAW
# SingleR label was MonacoImmune "T cells" / HPCA "T_cells" are assigned per cell:
# CD8 T when mean(CD8A,CD8B) > thr and exceeds mean(CD4,IL7R), otherwise CD4 T
# (same rule as build_reference.R step A). Human gene symbols; with no CD8A/CD8B
# in `expr` it returns `labels` unchanged (clean no-op on other species).
#   labels  character vector, already SINGLER_NORM-mapped
#   generic logical vector, TRUE where the raw label was generic T
#   expr    genes x cells log-normalised matrix, colnames aligned to labels
# -----------------------------------------------------------------------------
split_generic_t <- function(labels, generic, expr, thr = 0.5) {
  cd8g <- intersect(c("CD8A", "CD8B"), rownames(expr))
  if (!length(cd8g) || !any(generic)) return(labels)
  cd4g <- intersect(c("CD4", "IL7R"), rownames(expr))
  cells <- which(generic)
  m <- function(g, i) if (!length(g)) rep(0, length(i)) else
    as.numeric(Matrix::colMeans(expr[g, i, drop = FALSE]))
  cd8e <- m(cd8g, cells); cd4e <- m(cd4g, cells)
  labels[cells] <- ifelse(cd8e > thr & cd8e > cd4e, "CD8 T", "CD4 T")
  labels
}
