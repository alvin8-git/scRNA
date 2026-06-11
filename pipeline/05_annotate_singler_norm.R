# =============================================================================
# 05_annotate_singler_norm.R - SINGLER_NORM lookup mapping raw SingleR labels
# (HumanPrimaryCellAtlas + MonacoImmune) to canonical CELLTYPE_COLORS names.
# Sourced by 05_annotate.R. Pure data; no side effects.
# =============================================================================
SINGLER_NORM <- c(
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
