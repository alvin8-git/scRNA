# cm cohorts re-run with the 2026-10-26 pipeline — what changed

Commits `96c6047` (replicate-aware statistics, Leiden, per-layer HVGs, no mito regression
for cm, adaptive QC, SoupX step 00) and `bd3d673` (map guard, runner result dir), plus the
per-contrast fixes made during this run (see "Fixes made during the run").

Both cohorts were re-run from step 01 into new directories with the curated cluster maps
**disabled** (`SCRNA_NO_CLUSTER_MAP=1`), so labels come from SingleR (HumanPrimaryCellAtlas)
plus the Yap 2023 normalisation only. Conditions = timepoint. Old runs are untouched;
their `reports/` and `export/` were copied to `*_backup_20261026/` first.

| | Old run | New run |
|---|---|---|
| 7-sample SoupX | `Results/results_H1D01sx_7samples_filtered` | `Results/results_H1D01sx_7samples_leiden` |
| 8-sample raw | `Results/results_H1D01_8samples_filtered` | `Results/results_H1D01_8samples_leiden` |

Numbers below come from `/data/alvin/tmp/compare_runs.R` (cell-level join of the two
`integrated_annotated.rds` objects; outputs `cmp_7sx_*.csv`, `cmp_8_*.csv`) and from each
run's `composition_test.csv` and `differential/DE_summary.csv`.

## Verdict

The clustering changed, not only its numbering: the adjusted Rand index between old and
new cluster assignments is 0.40 (7-sample) and 0.33 (8-sample), with 9 of 22 and 5 of 19
new clusters mapping cleanly (≥ 80 %) onto one old curated label and 6 and 7 clusters
mixing two old labels. The cell-type labels changed far more (84 % and 94 % of cells) but
that is the disabled curated map, not the algorithm: without it SingleR's human atlas has
no cardiomyocyte, epicardial or mesenchyme class and files roughly half the cells as
"Neurons" and the endoderm as "Pluripotent", exactly the failure the curated map existed
to fix. The new runs are therefore a valid test of the clustering and statistics changes
and **not** a deliverable; a re-curated map on the Leiden clusters is needed before either
customer report can replace the old one.

## 1. Cells and clusters

| Cohort | Cells (identical old/new) | Clusters old → new (res 0.5) | ARI clusters | ARI labels | Cells with changed label |
|---|---|---|---|---|---|
| 7-sample SoupX | 76,314 | 22 → 22 | 0.404 | 0.248 | 83.5 % |
| 8-sample raw | 103,690 | 20 → 19 | 0.334 | 0.196 | 93.9 % |

QC was unchanged (fixed thresholds; adaptive mode is off by default), so cell sets are
identical and the comparison is cell-for-cell.

## 2. Do the new clusters still carve the old biology?

For each new cluster, the majority old (curated) label and its share of the cluster.

**7-sample:** median purity 78 %. Clean (≥ 80 %): Pluripotent ×4 (clusters 10, 11, 18,
19), Hepatic/Endoderm ×2 (9, 13), Endothelial (15, 99 %), Unknown (21), Proepicardial (22).
Mixed (< 60 %): cluster 3 Mesenchyme 36 % / Pluripotent 24 %; cluster 5 Proliferating 45 % /
Pluripotent 43 %; cluster 6 Myofibroblast 44 % / Mesenchyme 15 %; cluster 7 Mesenchyme 47 % /
Myofibroblast 39 %; cluster 2 Mesenchyme 60 % / Cardiac progenitor 26 %; cluster 17
Mesenchyme 58 % / Epicardial-like 18 %. Full table: `cmp_7sx_newcluster_purity.csv`.

**8-sample:** median purity 66 %. Clean: Hepatic/Endoderm ×3 (9, 12, 18), Endothelial (16,
99 %), Cardiomyocyte (15, 81 %). Mixed: cluster 3 Mesenchyme 37 % / Pluripotent 22 %;
cluster 5 Mesenchyme 45 % / Pluripotent 44 %; cluster 6 Pluripotent 54 % / Proliferating
34 %; cluster 14 Epicardial-like 47 % / Cardiac progenitor 30 %; cluster 10 Mesenchyme 51 % /
Pluripotent 31 %; cluster 11 Myofibroblast 55 % / Mesenchyme 35 %; cluster 4 Fibroblast
(low-depth) 58 % / Mesenchyme 19 %. Full table: `cmp_8_newcluster_purity.csv`.

Reading: the lineages with sharp expression programs (endoderm, endothelium, pluripotent,
cardiomyocyte) are recovered as whole clusters. The mesenchyme / epicardial / myofibroblast
/ cardiac-progenitor continuum is cut differently. Three changes act on exactly that
region: no %mt regression (the old regression removed part of the CM/progenitor axis),
per-sample HVG selection (fewer D0-dominated genes), and Leiden refinement. Which of the
three dominates was not isolated in this run [not computed].

## 3. Cell-type labels

Largest flows of shared cells, old curated label → new SingleR-only label (7-sample):

| old | new | cells | % of old |
|---|---|---|---|
| Epicardial-like | Neurons | 11,452 | 96.7 |
| Mesenchyme (lateral plate) | Neurons | 11,356 | 74.9 |
| Hepatic/Endoderm | Pluripotent | 7,314 | 98.1 |
| Fibroblast (low-depth) | Epithelial | 4,866 | 73.0 |
| Myofibroblast | Smooth Muscle | 4,000 | 87.5 |
| Cardiac progenitor | Neurons | 2,648 | 70.9 |
| Cardiomyocyte (immature) | Neurons | 2,504 | 69.8 |

Labels gone from the new runs: Amnion-like epithelium, Cardiac progenitor, Cardiomyocyte,
Epicardial-like, Fibroblast (low-depth), Hepatic/Endoderm, Mesenchyme (lateral plate),
Myofibroblast, Proepicardial, Proliferating. New: Neurons, Epithelial, Smooth Muscle. This
is the documented HPCA gap (`05_annotate_singler_norm.R` header), not a new defect.

## 4. Proportions

Per-type mean across samples, old → new, percentage points (only the big movers):

| Type | 7-sample old → new | 8-sample old → new |
|---|---|---|
| Neurons | 0 → 44.6 | 0 → 51.6 |
| Epicardial-like | 17.0 → 0 | 16.8 → 0 |
| Mesenchyme (lateral plate) | 14.9 → 0 | 13.6 → 0 |
| Fibroblast (low-depth) | 8.9 → 0 | 16.3 → 0 |
| Epithelial | 0 → 7.7 | 0 → 21.9 |
| Smooth Muscle | 0.3 → 11.9 | 0 → 9.2 |
| Pluripotent | 28.8 → 33.9 | 25.1 → 15.7 |
| Endothelial | 1.8 → 1.7 | 1.5 → 1.5 |

Endothelial is the one type both label systems agree on, and its proportion is unchanged
to 0.1 pp in both cohorts. Its time course is also the same in both label systems: absent
at D0, peak at D11, declining by D30. Sample-by-type changes over 5 pp: 44 of 126 (7-sample), 53 of
136 (8-sample); all are label-system effects (`cmp_*_composition.csv`).

## 5. New statistics (first use on cm)

**Composition test** (sample-level, limma on arcsin-sqrt proportions, FDR < 0.05).
8-sample, all four timepoints with 2 replicates: Pluripotent falls from 28 % (D0) to 4 %
(D20) and 5 % (D30), FDR 6e-6 to 2e-5; Endothelial peaks at D11 (4.8 % vs 0.04 % at D0,
FDR 9e-5) and drops again by D20/D30; Smooth Muscle rises D0 → D30 (2.7 % → 17.6 %, FDR
0.002); Epithelial rises D11 → D20 (4.7 % → 39 %, FDR 0.008). 7-sample (D20 dropped,
single sample): the same picture with sharper statistics because SoupX-corrected
replicates agree better. Pluripotent 78 % (D0) → 29 % (D11) → 7 % (D30), FDR 3e-7 and
3e-8; Endothelial 0.03 % → 5.0 % at D11 (FDR 1e-5) back to 0.6 % at D30; Smooth Muscle
2.3 % → 22 % by D30 (FDR 3e-6); "Neurons" (the SingleR bucket for cardiac lineages)
19 % → 53 % at D11 (FDR 2e-5).

**Pseudobulk DE** (edgeR QL, all pairwise timepoint contrasts). 8-sample: 24 of 24
contrasts tested; 7-sample: 9 of 18 (the 9 involving D20 are descriptive). Counts of
FDR < 0.05 and |log2FC| > 1 genes are in the thousands for every tested contrast (e.g.
"Neurons" D30 vs D0: 11,783 genes in the 8-sample run). That is expected when the label is
a grab-bag across a differentiation time course, and it is one more reason the SingleR-only
labels are not usable as delivered.

## 6. Fixes made during the run (uncommitted)

- `06b_differential.R`: a single-sample condition made the whole run descriptive; now only
  contrasts that touch it are, dispersion is estimated from the remaining replicates.
- `09_bootstrap_proportions.R`: same rule for the composition test (single-sample conditions
  are dropped from the design and named in the log).
- `08b_html_report.R`: a DE table with no tested contrasts crashed the report; the section
  is now omitted.

## 7. What to do next

1. Re-curate `CLUSTER_CELLTYPE_MAP` for both Leiden runs from
   `annotation/canonical_markers_dotplot.pdf` and `logs/05_annotate.log` (the log prints a
   paste-ready map), guarded on `CLUSTER$algorithm == 4`, then re-run `05 06 06b 07 09` and
   the reports. Only then compare proportions against the old customer report.
2. If the mixed mesenchyme/epicardial clusters are to be split, isolate which of the three
   normalisation changes is responsible by re-running step 04 with `NORM$vars_to_regress`
   set back to `"percent.mt"` for the cm overlay (cheap: 03 is cached per sample).
3. Decide whether the 1-sample D20 arm of the 7-sample cohort should be reported at all in
   the customer DE and composition sections, since every D20 contrast is descriptive.
