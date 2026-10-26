# cm cohorts re-run with the 2026-10-26 pipeline — what changed

Pipeline: commits `96c6047` (replicate-aware statistics, Leiden, per-layer HVGs, no %mt
regression for cm, adaptive QC, SoupX step 00), `bd3d673` (map guard, runner result dir),
`dd92a8d` (per-contrast descriptive rule). Cluster maps for the Leiden runs are in
`config_species_cm.R` (uncommitted at the time of writing), guarded on
`CLUSTER$algorithm == 4`; the Louvain maps are now guarded on `== 1`.

Both cohorts were re-run from step 01 into new directories. First with the curated maps
**off** (SingleR only, section 2), then re-curated against the Yap et al. 2023 rules in
`docs/cm_fresh_look_2026-10-24.md` (section 3 onward). Old runs are untouched; their
`reports/` and `export/` were copied to `*_backup_20261026/`.

| | Old run (Louvain, curated) | New run (Leiden, re-curated) |
|---|---|---|
| 7-sample SoupX | `Results/results_H1D01sx_7samples_filtered` | `Results/results_H1D01sx_7samples_leiden` |
| 8-sample raw | `Results/results_H1D01_8samples_filtered` | `Results/results_H1D01_8samples_leiden` |

Numbers come from `/data/alvin/tmp/compare_runs.R` (cell-level join of the two
`integrated_annotated.rds`; outputs `cmp3_7sx_*.csv`, `cmp4_8_*.csv`; the maps-off pass is
`cmp_7sx_*`, `cmp_8_*`), the per-cluster
evidence in `look_{7sx,8}.txt` and `sub*_*.txt`, and each run's `composition_test.csv`
and `differential/DE_summary.csv`.

## 1. Verdict

The clustering changed, not only its numbering. Cluster-assignment agreement (adjusted
Rand index) between old and new is 0.40 (7-sample) and 0.33 (8-sample) on identical cell
sets. Lineages with sharp programs (pluripotent, endoderm, endothelium, cardiomyocyte,
amnion, true epicardium) come back as whole clusters. The HAND1+ mesenchyme / epicardial
mesothelium / culture-activated fibroblast / cardiac-progenitor continuum is cut
differently, and Leiden also pulls cycling cells of several lineages into shared clusters.
After re-curation with Yap-rule subcluster splits, a third of cells carry a different label
than before (33 % and 36 %), the per-type means across samples move by at most 10 pp, and
every timepoint trend the old reports showed is reproduced with the new replicate-aware
statistics. The biggest honest change is at D11: cardiac progenitor is up (10 → 18 %),
lateral-plate mesenchyme down (41 → 15 %), with the balance moving into proliferating
mesenchyme and culture-activated fibroblast. The new customer reports are deliverable,
with the caveats in section 7.

## 2. Maps off: what SingleR alone does (both cohorts)

Without the curated map 84 % (7-sample) and 94 % (8-sample) of labels changed, because
HumanPrimaryCellAtlas has no cardiomyocyte, epicardial or mesenchyme class: it called
45–52 % of cells "Neurons" and the endoderm "Pluripotent". This is the documented HPCA gap
(`05_annotate_singler_norm.R` header), not a new defect, and it is why the maps exist.
Endothelial, the one type both systems agree on, was unchanged to 0.1 pp. Details in
`cmp_7sx.md` / `cmp_8.md`.

## 3. Clusters, old vs new

| Cohort | Cells | Clusters old → new (res 0.5) | ARI clusters | Cells whose curated label changed | ARI labels |
|---|---|---|---|---|---|
| 7-sample SoupX | 76,314 | 22 → 22 | 0.404 | 33.2 % | 0.459 |
| 8-sample raw | 103,690 | 20 → 19 | 0.334 | 32.8 % | 0.448 |

QC is unchanged (fixed thresholds), so the comparison is cell-for-cell.

**New clusters against old curated labels** (purity = share of the cluster carrying its
majority old label). 7-sample: median 78 %, 9 of 22 clusters ≥ 80 %, 6 below 60 %.
8-sample: median 66 %, 5 of 19 ≥ 80 %, 7 below 60 %. The impure clusters are the
mesenchyme/epicardial/fibroblast continuum and the cycling clusters; full tables in
`cmp3_7sx_newcluster_purity.csv`, `cmp4_8_newcluster_purity.csv`.

## 4. How the new maps were built (Yap et al. 2023 rules)

Vocabulary and thresholds are those of `docs/cm_fresh_look_2026-10-24.md`: committed
cardiac progenitor (CCP) = GATA4-high, TBX5 ≈ 35 %, NKX2-5 rising, WT1 low; lateral-plate
mesenchyme = HAND1/MSX1/TMEM88 with NKX2-5 < 10 %; cardiomyocyte = MYH6/SLC8A1/TNNT2/MYOCD;
culture-activated fibroblast (CDH11/FBN1/PDGFRB/ACTA2) = "Myofibroblast"; epicardial-like
mesothelium = UPK3B/BNC1 with WT1 low and TBX18 ≈ 0; true epicardium = WT1+ TBX18+ ALDH1A2+;
"Proliferating" is reserved for HAND1+ cycling mesenchyme; D0 is pluripotent, cell cycle is
not a lineage.

Where a Leiden cluster mixed two of these, it was split with `CLUSTER_SUBCLUSTER_MAP`
(FindSubCluster on the saved SNN graph, deterministic):

| Cohort | Cluster | Split into | Why |
|---|---|---|---|
| 7-sample | 2 (7,640) | CCP ×2, mesenchyme ×2 | CCP subclusters read GATA4 84, TBX5 36, NKX2-5 14, WT1 19 vs the paper's 83/37/13/16 |
| 7-sample | 1 (11,088) | epicardial-like ×2, mesenchyme | a D11 HAND1 88 UPK3B 11 subcluster inside a D30 mesothelium cluster |
| 7-sample | 3, 5 (12,861) | pluripotent, proliferating, CCP, epicardial-like | cycling cells from D0 hESC, D11 mesenchyme/CCP and D30 mesothelium |
| 8-sample | 5 (7,275) | pluripotent ×2, mesenchyme ×2 | D0 hESC and D11 mesenchyme co-clustered |
| 8-sample | 14 (3,402) | epicardial-like, CCP, myofibroblast, cardiomyocyte | four lineages |
| 8-sample | 2 (13,611) | epicardial-like, mesenchyme, pluripotent, **proepicardial** | recovers the WT1 41 / TBX18 38 / ALDH1A2 44 true-epicardium population (293 cells) the run otherwise lacked |
| 8-sample | 3, 10 (13,218) | pluripotent, proliferating, CCP, epicardial-like | cycling clusters, as 7-sample 3/5 |
| 8-sample | 1 (15,119) | fibroblast (low-depth), **pluripotent** | 3,323 shallow D0 hESC (POU5F1 70) co-clustered with the ambient-heavy fibroblast bulk on depth |
| 8-sample | 6 (6,925) | pluripotent ×4, proliferating, CCP | cycling hESC plus D11 mesenchyme and a GATA4 82 TBX5 31 progenitor subcluster |

Not split: the two endoderm clusters in each cohort hold 11–20 % D0 cells that are spread
across every subcluster and express FOXA2 themselves (spontaneous endoderm in the hESC
culture, or ambient). No graph cut separates them; they stay Hepatic/Endoderm and that
is why D0 now shows 7.5 % endoderm (7-sample) where the old run folded those cells into a
pluripotent cluster.

## 5. Proportions, old vs new — 7-sample SoupX cohort (final)

Per-type mean across samples, percentage points:

| Type | Old | New | Δ |
|---|---|---|---|
| Pluripotent | 28.8 | 24.5 | −4.4 |
| Epicardial-like | 17.0 | 20.0 | +3.0 |
| Mesenchyme (lateral plate) | 14.9 | 6.2 | −8.6 |
| Hepatic/Endoderm | 7.2 | 9.6 | +2.4 |
| Myofibroblast | 5.7 | 11.6 | +5.9 |
| Proliferating | 4.8 | 7.5 | +2.7 |
| Cardiac progenitor | 3.6 | 5.9 | +2.3 |
| Cardiomyocyte (immature) | 3.7 | 3.1 | −0.6 |
| Endothelial | 1.8 | 1.7 | 0.0 |
| Amnion-like epithelium | 1.6 | 1.7 | +0.2 |
| Proepicardial | 0.3 | 0.3 | +0.1 |

Sample-level moves over 5 pp (17 of 112): D0 Pluripotent 98.6 → 79.5 (the endoderm and
cycling-cell effects above); D11 Mesenchyme 42.5/40.1 → 15.5/14.4; D11 Cardiac progenitor
10.5/10.7 → 17.1/18.9; D11 Myofibroblast 1.8/1.6 → 10.9/10.1; D11 Proliferating 7.8/7.3 →
15.9/14.5; D30 Myofibroblast 16.1/14.1 → 22.4/20.2.

**Against Yap's day-11 numbers** (CCP 33 %, mesenchymal + fibroblast 57 %, epithelial 8 %,
endothelial 1.8 %): new D11 reads CCP + immature CM 21–25 %, mesenchyme + fibroblast +
proliferating ≈ 40 %, endoderm 20 %, endothelial 5 %. Closer on CCP than the old run
(14–16 %), still below the paper; the endoderm share remains our main divergence from
the paper's "epithelial 8 %".

Composition test (sample-level, D20 dropped as single-sample; FDR < 0.05, 15 rows):
Pluripotent 80 % → 2 % by D11 (FDR 2e-8); Cardiac progenitor 0.2 % → 18 % at D11 then
0.9 % at D30 (FDR 1e-6, 2e-6); Mesenchyme peaks D11 (15 %) and is gone by D30; Endothelial
peaks D11 (5 %); Epicardial-like 2 % → 37 % and Myofibroblast 2 % → 21 % by D30; Hepatic/
Endoderm 21 % at D11 down to 4 % at D30. Pseudobulk DE: 9 of 18 contrasts tested
(D20 descriptive), thousands of FDR < 0.05 genes per tested contrast as expected across
a differentiation time course.

## 6. Proportions, old vs new — 8-sample raw cohort (final)

Per-type mean across samples, percentage points:

| Type | Old | New | Δ |
|---|---|---|---|
| Pluripotent | 25.1 | 24.4 | −0.6 |
| Epicardial-like | 16.8 | 19.1 | +2.3 |
| Fibroblast (low-depth) | 16.3 | 17.3 | +1.0 |
| Mesenchyme (lateral plate) | 13.6 | 3.9 | −9.7 |
| Myofibroblast | 7.4 | 9.9 | +2.5 |
| Hepatic/Endoderm | 6.2 | 7.1 | +0.8 |
| Cardiac progenitor | 3.5 | 5.9 | +2.4 |
| Cardiomyocyte | 3.9 | 3.1 | −0.8 |
| Proliferating | 3.4 | 5.4 | +2.0 |
| Endothelial | 1.5 | 1.5 | 0.0 |
| Amnion-like epithelium | 1.4 | 1.8 | +0.4 |
| Proepicardial | 0.2 | 0.3 | 0.0 |

Sample-level moves over 5 pp (16 of 104): D11 Mesenchyme 38.3/35.3 → 13.9/12.2; D0
Pluripotent 98.7 → 82.6/81.9 (endoderm 3.5 %, cycling 3 %, and small ambient-driven
calls of 1–3 % each); D11 Cardiac progenitor 9.5/11.9 → 17.6/19.8; D11 Proliferating
6.9/6.3 → 13.8/13.2; D11 Myofibroblast 3.8/3.5 → 10.5/9.7; H1D20_1 Pluripotent 0.9 → 9.5
(the ambient-heavy sample; these cells sit in the cycling cluster 6_0 with the D0 cells).
The same D11 shift as the 7-sample cohort, and the true-epicardium class is back (0.3 %).

Composition test (all four timepoints, 2 replicates each; FDR < 0.05, 13 rows):
Pluripotent 82 % → 5 % at D11 (FDR 9e-7); Cardiac progenitor 0.2 % → 19 % at D11 then
2 % at D20 (FDR 5e-6, 1e-5); Mesenchyme 1 % → 13 % at D11 then < 1 % (FDR 5e-6); Hepatic/
Endoderm peaks D11 at 20 % (FDR 5e-6); Endothelial peaks D11 (FDR 1e-5). Pseudobulk DE:
24 of 24 contrasts tested.

D11 against Yap (CCP 33 %, mesenchymal + fibroblast 57 %): CCP + cardiomyocyte 21 %,
mesenchyme + fibroblast + proliferating 42 %, endoderm 19 %, endothelial 5 %, residual
pluripotent 6 %.

## 7. Caveats carried into the new reports

- D0 shows 7.5 % Hepatic/Endoderm and 4 % Proliferating (7-sample) because of the
  unsplittable endoderm clusters and the cycling subclusters; the old reports showed 98.6 %
  Pluripotent at D0. State this in the delivery note.
- The 8-sample cohort still has no ambient correction; its two low-depth fibroblast
  clusters (27 % of cells) are where ambient RNA and H1D20_1 collect, as before.
- "Myofibroblast" means culture-activated fibroblast (Yap doc), not a smooth-muscle call.
- D20 in the 7-sample cohort has one SoupX-corrected sample, so every D20 contrast in DE
  and the composition test is descriptive.
- The per-sample SoupX ambient table (`qc/soupx_rho.csv`: H1D20_1 rho 0.795, the
  exclusion criterion; the other seven 0.05–0.27) is written by the SoupX step into the
  run directory it ran in. The fresh `_leiden` directories did not have it, so the first
  rebuild of their reports lacked the Ambient RNA section. It was copied over from the
  old run directories (same samples, rho is a per-sample quantity) and the reports and
  exports rebuilt. A re-run into a new directory must carry this file across until the
  SoupX step stores it beside the corrected matrices.
- The "low QC" tag on raw H1D20_2 was wrong, twice. Its retention is 45.7 % because all
  17,834 called barcodes pass the 200-gene load gate and 8,677 of them then fail the cm
  1,000-gene floor; the after-load denominator (commit 790c5f5) therefore changed
  nothing. The surviving 9,157 cells have 2,863 median genes and rho 0.09: a fine
  library behind an over-permissive caller. H1D20_1 (rho 0.795, 1,300 median genes)
  retained 92 % and was never tagged. The tag now uses median genes < 500 or ambient
  rho > 0.3, and retention is shown but not judged (regression guard T13).

## 8. Fixes made during the run

Committed in `dd92a8d`: per-contrast descriptive rule in 06b and 09; 08b no longer crashes
on a DE table without tested contrasts. Uncommitted: the Leiden maps and subcluster specs
in `config_species_cm.R`, and the `algorithm == 1` guards on the Louvain maps.
