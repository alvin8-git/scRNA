# Cardiac cohorts: fresh look at the labels against Yap et al. 2023 (npj Regen Med 8:26)

Date: 2026-10-24. Evidence scripts: `/data/alvin/tmp/fresh_look.R`, `fresh_look2.R`, `raw_qc_cm.R`;
outputs in `/data/alvin/tmp/fresh_look/`. Paper text: `/data/alvin/tmp/paper/paper.txt`.

## What the paper (same lab, same H1 line, same LN521+221 protocol) actually says

| Fact | Consequence for our labels |
|---|---|
| Day 0 = seeding of pluripotent H1 on LN521+221; CHIR 10 µM given at **day 4** (confluence) for 24 h; IWP2 day 7–9. | D0 is pre-induction pluripotent (TBXT/MIXL1/EOMES 0% in both D0 samples — confirmed). **D11 = 7 days post-CHIR**, D20 = 16 dpc, D30 = 26 dpc. We had been treating the D-number as days of differentiation. |
| Day 7 = ISL1/MESP1-peak "CVPs"; day 11 = NKX2-5/GATA4-rising **committed cardiac progenitors (CCPs)**, non-contracting; first action potentials day 15–17. | The lab's "CCP" is what our pipeline calls Cardiomyocyte (immature) plus the TBX5+/GATA4+ progenitor cluster. Their D11 product is a progenitor by function, a CM by transcript. |
| CCP signature: MYH6, SLC8A1, ACTN2, ACTC1, ANKRD1, CRHBP, TNNT2, IGFBP7, CCDC80, TNNI1, MYL4. | Cluster 12 (7-sample) / cluster 9 (8-sample) carry exactly this set (MYH6 81%, SLC8A1 83%, TTN 96%, ANKRD1 56%). |
| Day-11 scRNA-seq (10x v2, ~3,000 cells × 2 replicates, Seurat v3 tSNE): **CCP 33.2%, mesenchymal + fibroblast 56.7%, epithelial 8%, endothelial 1.8%**. No endoderm, epicardial, amnion or pluripotent class was named. | Our D11 mesenchymal fraction matches theirs; our "epithelial" is three times larger and resolves as endoderm; our CCP is lower (see below). |
| D20/D30 were never characterised in vitro — the paper transplants D11 cells. | There is no lab vocabulary for the late time points; those labels are ours alone. |

## What was assumed or missed

1. **The 8-sample cohort has never been curated.** `results_H1D01_8samples_filtered` carries raw SingleR
   labels: Neurons 72–93% of every differentiated sample, "Pluripotent" 20% at D11 (= endoderm clusters 7/13),
   "Smooth Muscle" = fibroblasts. The 11 MB `Customer_report.html` built on 2026-10-24 contains 28,455
   "Neurons" strings and zero "Cardiomyocyte". Both archived maps in `config_species_cm.R` are guarded on the
   earlier 7-name set that *excluded* H1D11_1, so nothing matches the 8-name run.
2. **"Cardiac progenitor" (clusters 0/4, 42% of D11) is not a cardiac progenitor.** NKX2-5 3%, TBX5 5%,
   ISL1 6%. It is HAND1 89–92% / MSX1 65% / TMEM88 89% / BMP4 70% / LUM 85% lateral-plate mesenchyme — the
   paper's "mesenchymal" population and the precursor of the D20/D30 fibroblast/epicardial compartment.
3. **"Proepicardial" (cluster 10, 10.5% of D11) is the actual cardiac progenitor.** WT1 16%, TBX18 11%,
   ALDH1A2 2% (no proepicardial signal); GATA4 83%, TBX5 37%, NKX2-5 13%, MYL7 82%, SLC8A1 41%, CCDC80 79%;
   nearest clusters by pseudobulk are 0 (0.963) and 12 (0.941). This is the lab's CCP, first-heart-field
   flavour. True CCP at D11 = clusters 10 + 12 ≈ 18–20% vs the paper's 33.2%.
4. **"Unknown" at D20/D30 (cluster 2, 18–28% of D30) is low-depth fibroblast-lineage cells, not a cell type
   and not hepatocytes.** Hepatic genes are detected in 70–97% of cells but carry only 1.1% of UMIs per cell
   (IQR 0.8–1.3, uniform = ambient); the genuine endoderm cluster 13 carries 3.4% (IQR 1.9–5.9) and 7.8% at
   pseudobulk. Fibroblast UMI share 1.37% matches clusters 1/7; nearest cluster is 5 (0.94). Label it
   Fibroblast (low-depth) rather than hiding a quarter of D30 under Unknown. Clusters 20/21 = doublets/dying.
5. **"Epithelial" (cluster 16) is amnion-like epithelium**: GABRP 54%, VTCN1 29%, TFAP2A 40%, GATA3 35%,
   ISL1 37%, WFDC2 69%, GRHL2 53%, EPCAM 74%. A classic BMP4-driven off-target, 2% of D11 and D30.
6. **Genuine epicardium is tiny.** Only cluster 18 (215 cells, D30) is WT1 60% / TBX18 71% / ALDH1A2 82%.
   Clusters 1/7 (12k cells, 50% of D20, 31–37% of D30) are WT1 14–24%, TBX18 1–2%, UPK3B 41–69%, BNC1 63–66%,
   KRT19 95%, CDH1 56–64%: epicardial-like mesothelium / epicardium-derived cells, HAND1+ and TMEM88+.
   "Epicardial" is defensible; "Epicardial-like (WT1-low)" is honest.
7. **Cardiomyocyte yield at D20/D30 is 1–3% and that is real, not QC loss.** Pre-QC matrices: 11–22 strong
   CM barcodes (TNNT2≥5 & MYH6≥3) per D30 sample out of 13–14k; the mito gate removes none (median mito
   <1.5% in every sample — this chemistry barely captures mt reads, so `max_percent_mt` is inert here);
   the only CM-like barcodes removed fall under `min_features` 1000. SoupX reduces D30 CM-like barcodes
   from 633 to 194, i.e. two thirds of the raw TNNT2 signal at D30 was ambient.
8. **Nothing neural, haematopoietic, skeletal-muscle or trophoblast exists** (PAX6/SOX10/HBG1/MYOD1 ≤1%).
   SingleR's "Neurons" was always an artefact of a reference with no cardiomyocyte class.
9. **Endoderm is two populations**, both already labelled Hepatic/Endoderm: cluster 6 is hepatic-specified
   foregut (FOXA2 75%, HHEX 55%, TBX3 54%, HNF1B 20%), cluster 13 is yolk-sac/visceral-like (APOB 86%,
   MTTP 82%, RBP4 87%, FGB 86%, HNF4A 54%). Cluster 19 is its proliferating fraction. 20% of D11, 3% of D30.
10. The fibroblast resting/activated split inside cluster 5 is thin (TAGLN 72%, ACTA2 70% cluster-wide;
    CDH11 86%, SNAI2 69%, FBN1 72%). The whole cluster is culture-activated mesenchymal fibroblast.

## Actual cell types — 7-sample corrected cohort (`results_H1D01sx_7samples_filtered`, 22 clusters)

| Cl | n | Current | Evidence-based identity | Decisive genes (% det) |
|---|---|---|---|---|
| 3, 9, 11 | 13,845 | Pluripotent | Pluripotent (11 = low-UMI fraction) | POU5F1 92/90/36, SOX2 83/63/32, ESRG, XACT |
| 0, 4 | 12,961 | Cardiac progenitor | **Lateral-plate mesenchyme (HAND1+)** | HAND1 89/92, MSX1 65/67, TMEM88 89/93, NKX2-5 3/3 |
| 14 | 2,208 | Cardiac progenitor | Mesenchyme, low depth | HAND1 46, GATA4 56, CDH11 58, nearest cl 7 |
| 10 | 3,737 | Proepicardial | **Cardiac progenitor / CCP (TBX5+ GATA4+)** | GATA4 83, TBX5 37, NKX2-5 13, MYL7 82, SLC8A1 41, WT1 16 |
| 12 | 3,703 | Cardiomyocyte (immature) | Cardiomyocyte, immature = CCP signature | MYH6 81, SLC8A1 83, TNNT2 69, NKX2-5 56, NPPA 22, MYL2 1 |
| 8 | 3,878 | Proliferating | Proliferating mesenchyme | TOP2A 93, HAND1 79, BMP4 80 |
| 6, 19 | 4,915 | Hepatic/Endoderm | Hepatic-specified foregut endoderm (19 proliferating) | FOXA2 75, HHEX 55, TBX3 54, TTR 74 |
| 13 | 2,540 | Hepatic/Endoderm | Yolk-sac / visceral endoderm-like | APOB 86, MTTP 82, RBP4 87, FGB 86, HNF4A 54 |
| 15 | 1,832 | Endothelial | Endothelial (endocardial-leaning) | CDH5 65, PECAM1 45, CD34 52, NPR3 8, HAPLN1 82 |
| 16 | 1,364 | Epithelial | **Amnion-like epithelium** | GABRP 54, VTCN1 29, TFAP2A 40, GATA3 35, ISL1 37 |
| 1, 7 | 11,846 | Epicardial | Epicardial-like mesothelium (WT1-low) | UPK3B 69/41, BNC1 63/66, KRT19 95/80, WT1 14/24, TBX18 1/2 |
| 18 | 215 | Proepicardial | Epicardium (true) | WT1 60, TBX18 71, ALDH1A2 82, UPK3B 73 |
| 5 | 5,229 | Myofibroblast + 2 subtypes | Fibroblast, culture-activated | CDH11 86, SNAI2 69, FBN1 72, PDGFRB 65, ACTA2 70 |
| 17 | 1,106 | Myofibroblast/SMC/Prolif | Smooth-muscle-like / stressed | ACTA2 89, TAGLN 91, CNN1 64, NPPB 46, ANKRD1 62 |
| 2 | 6,668 | Unknown | **Fibroblast-lineage, low depth** (ambient hepatic reads) | DCN 93, LUM 96, hepatic UMI 1.1%/cell, nearest cl 5 |
| 20, 21 | 267 | Unknown | Doublets (20) / dying (21, mito 10.7%) | — |

Re-mapped composition: D11 ≈ 18–20% CCP, 53% mesenchyme (incl. proliferating), 20% endoderm, 5% endothelial,
2% amnion. D30 ≈ 35–45% fibroblast (cl 5 + 2 + 17), 31–37% epicardial-like, 3–7% CM/progenitor,
3% endoderm, 2% amnion, <1% endothelium, 2% true epicardium (D30 only).

## Actual cell types — 8-sample raw cohort (`results_H1D01_8samples_filtered`, 20 clusters, uncurated)

| Cl | n | SingleR label now | Identity (same evidence) |
|---|---|---|---|
| 4, 8, 11 | 15,761 | Pluripotent | Pluripotent |
| 1, 3 | 17,673 | Neurons | Lateral-plate mesenchyme (HAND1 87/74, MSX1 62/39, NKX2-5 5/6) |
| 10 | 4,449 | Neurons | Cardiac progenitor / CCP (GATA4 79, TBX5 31, NKX2-5 11, WT1 17) |
| 9 | 4,717 | Neurons | Cardiomyocyte, immature (MYH6 73, SLC8A1 76, NKX2-5 48) |
| 12 | 4,157 | Neurons | Proliferating mesenchyme |
| 7 | 5,336 | Pluripotent | Hepatic-specified endoderm (FOXA2 75, HHEX 55) |
| 13 | 2,813 | Pluripotent | Yolk-sac endoderm-like (APOB 85, MTTP 79) |
| 14 | 2,034 | Endothelial | Endothelial |
| 15 | 1,600 | Neurons | Amnion-like epithelium (GABRP 55, TFAP2A 40, ISL1 38) |
| 2, 5 | 16,869 | Neurons | Epicardial-like mesothelium (UPK3B 71/41, BNC1 66/61, WT1 15/20) |
| 18 | 274 | Smooth Muscle | Epicardium, true (WT1 59, TBX18 66, ALDH1A2 75) |
| 6 | 6,660 | Smooth Muscle | Fibroblast (CDH11 88, FBN1 71, PDGFRB 67) |
| 16 | 1,428 | MSC | Smooth-muscle-like / stressed (CNN1 67, NPPB 43) |
| 0 | 19,340 | Neurons | Fibroblast-lineage, low depth + ambient (1,220 genes; 56% D20 incl. H1D20_1; hepatic UMI 0.9%/cell) |
| 17 | 372 | Epithelial | Dying (mito 6.6%) |
| 19 | 207 | Fibroblast (resting) | Doublets / ambient-heavy (hepatic UMI 1.1%, uniform) |

Cluster 0 is 32% of the cohort and is where the uncorrected ambient RNA (and H1D20_1, rho 0.795) collects;
in the SoupX cohort the equivalent cluster 2 is 9%. That difference, not biology, is the main reason the
two cohorts look different.

## Recommendation

- Pull or regenerate the 8-sample customer report; write an 8-name-guarded `CLUSTER_CELLTYPE_MAP` from the
  table above before any customer sees it.
- In the 7-sample map: rename 0/4/14 → "Mesenchyme (lateral plate)", 10 → "Cardiac progenitor", 16 →
  "Amnion-like epithelium", 2 → "Fibroblast (low-depth)", 1/7 → "Epicardial-like"; consider collapsing the
  resting/activated fibroblast split. "Cardiomyocyte (immature)" can be reported to the lab as "CCP" at D11.
- Report D-numbers as days from seeding with CHIR at day 4, matching the paper.
- Frozen-reference verdict is unchanged: a reference built from mislabelled anchors cannot transfer; fix the
  labels first, then revisit.
