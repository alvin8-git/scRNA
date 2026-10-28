#!/usr/bin/env Rscript
# mouse_reint.R <run_dir> <out_txt> — evidence for re-interpreting the mouse design:
# tissue of origin (heart vs peritoneum vs blood), post-MI macrophage/neutrophil states,
# T-cell content (nude mice are athymic), ambient cardiac signal, and any human-CVP trace
# (mouse-only reference: conserved cardiac genes are the only possible footprint).
suppressMessages({ library(Seurat); library(Matrix) })
a <- commandArgs(trailingOnly = TRUE); RUN <- a[1]; OUT <- a[2]
sink(OUT); options(width = 240)
o <- readRDS(file.path(RUN, "integrated/integrated_annotated.rds")); md <- o@meta.data
cm <- SeuratObject::LayerData(o, assay = "RNA", layer = "counts"); rm(o); gc()
cl <- as.character(md$seurat_clusters); cls <- as.character(sort(as.integer(unique(cl))))
smp <- as.character(md$sample); sms <- unique(smp)
lab <- sapply(cls, function(k) names(which.max(table(md$cell_type[cl == k]))))
P <- list(
  tissue_heart_resident_mac = c("Timd4","Lyve1","Folr2","Cd163","Mrc1","Igf1"),
  tissue_peritoneal_LPM     = c("Gata6","Prg4","Icam2","Saa3","Ltbp1"),
  tissue_alveolar_or_spleen = c("Siglecf","Marco","Car4","Vcam1","Spic"),
  recruited_mono_mac        = c("Ccr2","Ly6c2","Plac8","Chil3","Fn1","Arg1"),
  post_MI_repair_mac        = c("Spp1","Trem2","Gpnmb","Fabp5","Lgals3","Cd9"),
  neutrophil                = c("S100a8","S100a9","Retnlg","Cxcr2","Ly6g","Mmp8"),
  T_cell_core               = c("Cd3e","Cd3g","Cd4","Cd8a","Cd8b1","Trac","Trbc2","Lck"),
  NK_ILC_gdT                = c("Ncr1","Klrb1c","Gzma","Trdc","Tcrg-C1","Il7r","Rorc","Gata3"),
  stroma_heart              = c("Col1a1","Dcn","Pdgfra","Tcf21","Pecam1","Cdh5","Kdr"),
  cardiomyocyte_ambient     = c("Mb","Myl2","Myl3","Tnnt2","Actc1","Myh6","Nppa","Ttn"),
  cvp_conserved             = c("Gata4","Myl7","Hand1","Tbx5","Isl1","Mdk"),
  blood_rbc                 = c("Hba-a1","Hbb-bs","Alas2"))
P <- lapply(P, function(g) g[g %in% rownames(cm)])
det <- function(i, g) round(100 * Matrix::rowMeans(cm[g, i, drop = FALSE] > 0))
cat("cells", ncol(cm), "\n\n=== detection % per CLUSTER (cluster label on top) ===\n")
for (nm in names(P)) { if (!length(P[[nm]])) next
  m <- sapply(cls, function(k) det(cl == k, P[[nm]])); if (is.null(dim(m))) m <- t(m)
  colnames(m) <- paste0(cls, ":", substr(lab, 1, 6)); cat("\n--", nm, "\n"); print(m) }
cat("\n\n=== detection % per SAMPLE ===\n")
for (nm in names(P)) { if (!length(P[[nm]])) next
  m <- sapply(sms, function(s) det(smp == s, P[[nm]])); if (is.null(dim(m))) m <- t(m)
  cat("\n--", nm, "\n"); print(m) }
cat("\n\n=== UMI share (%) per sample: cardiac ambient, RBC, ECM ===\n")
sets <- list(cardiac = P$cardiomyocyte_ambient, rbc = P$blood_rbc, ecm = c("Col1a1","Col3a1","Dcn"))
print(round(t(sapply(sms, function(s) { i <- which(smp == s); tot <- sum(cm[, i])
  sapply(sets, function(g) 100 * sum(cm[g[g %in% rownames(cm)], i]) / tot) })), 3))
cat("\n=== cells co-expressing >= 3 of Tnnt2/Actc1/Myl7/Gata4/Myl4/Ttn (possible human-CVP or CM footprint) ===\n")
cg <- intersect(c("Tnnt2","Actc1","Myl7","Gata4","Myl4","Ttn"), rownames(cm))
n3 <- Matrix::colSums(cm[cg, ] > 0) >= 3
print(data.frame(sample = sms, cells = as.integer(table(smp)[sms]),
                 n_cardiac3 = sapply(sms, function(s) sum(n3 & smp == s)),
                 pct = round(100 * sapply(sms, function(s) mean(n3[smp == s])), 2)), row.names = FALSE)
cat("\n=== T-cell cluster: is it real T (Cd3e+Trac+) or NK/ILC labelled CD4 T? ===\n")
tk <- cls[grepl("T", lab) & !grepl("Mono", lab)]
for (k in tk) { i <- cl == k; cat(sprintf("cl %s (%s, n=%d): Cd3e %d Trac %d Trbc2 %d Cd4 %d Cd8a %d Trdc %d Ncr1 %d Il7r %d Rorc %d Gata3 %d\n", k, lab[k], sum(i),
  det(i,"Cd3e"), det(i,"Trac"), det(i,"Trbc2"), det(i,"Cd4"), det(i,"Cd8a"), det(i,"Trdc"), det(i,"Ncr1"), det(i,"Il7r"), det(i,"Rorc"), det(i,"Gata3"))) }
cat("\nREINT_DONE\n"); sink()
