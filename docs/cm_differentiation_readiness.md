# hiPSC/hESC → cardiomyocyte differentiation: data assessment and `cm` mode

Assessment of `Samples/Cardio` and the design of `SCRNA_SPECIES=cm`. Written 2026-09-10
against the 7-sample H1 time course, before the first pipeline run completed. Companion to
`docs/bat_wing_readiness.md`, which does the same job for wing tissue.

---

## 1. The dataset

| sample | timepoint | cells | median genes | median UMI |
|---|---|---|---|---|
| H1D0_1  | D0  | 10,558 | 2,124 | 5,014 |
| H1D0_2  | D0  | 11,004 | 2,700 | 7,528 |
| H1D11_2 | D11 | 17,109 | 3,208 | 10,143 |
| H1D20_1 | D20 | 19,279 | **1,300** | **2,977** |
| H1D20_2 | D20 | 17,834 | **1,019** | **2,140** |
| H1D30_1 | D30 | 16,715 | 2,822 | 7,531 |
| H1D30_2 | D30 | 15,159 | 3,508 | 11,351 |

**107,658 cells**, 37,488 genes, identical gene universe across samples, GRCh38 with Ensembl
IDs + symbols. MGI DNBelab C4 (consistent with the platform noted on `hiPSC.jpg`). D11 has no
replicate. Cell line H1 is the hESC line WA01 — note the name collides visually with the
repo's bundled `Samples/H1` human PBMC fixture, but the derived sample names (`H1D0_1`, …)
do not clash.

---

## 2. Two data-quality facts that shape everything

### 2.1 D20 is 2.5–3× shallower than the rest

Median genes per cell: **1,019–1,300 at D20** vs 3,208 (D11) and 3,508 (D30_2). Every
cardiomyocyte marker appears to *dip* at D20 and recover at D30 — TNNT2 30.1 → 15.2 → 37.5,
TTN 44.6 → 21.4 → 65.2, ACTC1 83.1 → 19.2 → 52.5. **That dip is at least partly library
depth, not loss of cardiomyocytes.** Do not describe D20 as a regression in differentiation
without first running step 09 (depth-normalised bootstrap proportions) and step 10
(rarefaction). This is the single most likely source of a wrong conclusion from this cohort.

### 2.2 Heavy ambient RNA at D20/D30

Detection rate in **essentially every cell** at D30: COL1A1 100.0%, COL3A1 100.0%, LUM 99.8%,
DCN 99.1%, **AFP 99.4%**, FN1 96.7%. Near-universal detection of secreted or highly abundant
transcripts is the ambient-RNA ("soup") signature, not biology.

This is the same failure mode that made bat neutrophils unrecoverable from raw detection
rates: shared markers plus a raised background collapses one cell type into another. Score
these populations on **soup-corrected** signal — per-cluster mean minus the mean in a
compartment that cannot express the gene (e.g. collagen inside pluripotent cells at D0) —
never on detection rate alone.

---

## 3. What the markers say about the biology

Per-cell detection rates, D0 → D11 → D20 → D30 (replicates averaged):

| state | evidence |
|---|---|
| **D0 is cleanly pluripotent** | POU5F1 89.2 → 1.4 → 0.2 → 0.5; SOX2 69.4 → 0.1; LIN28A 78.3, DNMT3B 88.0 |
| **D11 is the progenitor window** | NKX2-5 0.0 → **18.5** → 2.1 → 6.7; TBX5 0.0 → **19.7**; ISL1 0.5 → **15.5**; GATA4 0.3 → **71.2**; TBX3 → 32.3. Largely gone by D20 |
| **CM programme engages and holds** | MYL7 0.9 → 89.8; TNNI1 0.5 → 72.5; TTN 8.4 → 65.2; MYL4 → 54.5; MYH6 0.1 → 27.3; MYH7 0.0 → 21.7 |
| **The CMs are immature / non-ventricular** | **MYL2 1.4% at D30**, IRX4 5.0%, HEY2 8.3% — ventricular identity essentially absent. Atrial/immature markers dominate: MYL7 89.8, NR2F2 23.6, NPPA 11.4 |
| **Maturation is beginning** | PLN 0.2 → 33.2, MYBPC3 0.1 → 16.2, CACNA1C 13.9 → 34.8 |
| **Epicardium emerges** | WT1 0.1 → **18.4**, ALDH1A2 0.9 → 9.7, TBX18 → 3.8 — a population absent from the source figures |
| **Endothelium is minimal** | PECAM1 1.1%, CDH5 1.4% at D30 (APLNR peaks 15.5% at D11) — unlike the Day-21 figure, which shows a clear endothelial cluster |
| **True smooth muscle is scarce** | MYH11 4.1%, ANPEP 2.9%. PDGFRB 42.0% is more plausibly fibroblast than pericyte here |
| **Substantial hepatic/endoderm off-target** | AFP 6.4 → 12.8 → **95.6** → **99.4**; FOXA2 peaks 23.4% at D11; SOX17 2.0% at D11. Subject to §2.2 — magnitude unresolved until soup-corrected |

### Anomalies worth a second look

- **TNNI3 runs backwards**: 38.5% at D0 → 12.1% at D30. TNNI3 is the adult cardiac troponin I
  and should *rise* with maturation; 38.5% in undifferentiated hESCs is not credible. Suspect
  ambient or mis-assignment; do not use TNNI3 as a maturation marker in this cohort.
- **SLN peaks at D20** (0.1 → 27.3 → 34.9 → 4.7) then collapses. Non-monotonic.
- **EPCAM falls** 51.2 → 17.5, i.e. most of the "epithelial" signal is the pluripotent
  compartment (hESCs are EPCAM+), not a distinct late epithelium.

---

## 4. Marker panel: what was supplied vs what the data supports

Supplied: ACTC1, TNNT2, TTN, MYH6, MYL4, CCDC80, CACNA1C, plus the figure panel (ANPEP,
PDGFRA, PDGFRB, APLNR, ENG, PECAM1, TTN, TNNT2, TNNI3, TNNI1, TNNC1, ANKRD1, ACTN2, ACTC1,
MYH7, MYH6, MYL4, NKX2-5, MESP1, GATA6, GATA4, ISL1, EPCAM, ACTA2, FN1, LUM, SOX2, NANOG).
All present in this annotation.

**Added, with the evidence:**

| added | why |
|---|---|
| **MYL7** | 89.8% at D30 — the strongest CM discriminator here, better than the supplied MYL4 (54.5%) |
| **TNNI1** | 72.5% — the correct troponin I for immature CM, and TNNI3 is unreliable (§3) |
| **MYL2, IRX4, HEY2** | ventricular identity. Included *so their absence is visible*, which is the finding |
| **NPPA, NR2F2, SLN, GJA5, KCNJ3** | atrial/immature identity, which is what these cells actually are |
| **POU5F1, LIN28A, DNMT3B** | POU5F1 89.2% at D0 vs NANOG 27.4% — the supplied SOX2/NANOG pair understates the pluripotent fraction |
| **PLN, MYBPC3, RYR2, ATP2A2, SCN5A** | maturation / Ca-handling axis; PLN rises 0.2 → 33.2 |
| **TBX5, MESP1, KDR** | resolves the D11 progenitor window |
| **WT1, TBX18, ALDH1A2, UPK3B** | epicardium — a real population missing from the figures |
| **AFP, FOXA2, SOX17, TTR, APOA1** | the endoderm off-target, the largest unexplained signal in the data |
| **MKI67, TOP2A, CCNB1** | cycling cells (TOP2A 60.8% at D0 → 14.2% at D20 → 33.4% at D30) |
| **RGS5, KCNJ8, NOTCH3** | pericyte markers that are specific, unlike PDGFRB alone |
| **COL1A1, COL3A1, DCN, POSTN, TCF21** | fibroblast compartment, the dominant non-CM population |

**Recategorised:** **CCDC80** was supplied as a cardiac marker. It tracks the fibroblast
compartment here (4.9 → 61.5 → 28.9 → 68.3) and is grouped under Fibroblast, not
Cardiomyocyte.

**Checked and rejected as absent from the biology:** MYOD1/MYOG/PAX7 (skeletal muscle),
PAX6/SOX10 (neural crest/neuroectoderm), RAX/SIX6/VSX2 (retina) — all present in the
annotation but not expressed, so the usual iPSC off-target lineages are *not* a problem here.
The off-target is endodermal, not neural.

---

## 5. Cell types to expect beyond the source figures

The figures label Cardiomyocytes, Epithelial, Fibroblast, Fibroblast/Epi, Pericytes,
Endothelial. On this data also expect:

1. **Pluripotent / residual undifferentiated** — dominates D0; any persistence past D11 is a
   quality signal
2. **Cardiac progenitor** (NKX2-5⁺/ISL1⁺/TBX5⁺) — expected as a transient D11 state, but **no
   such cluster was found** in the min_features 1000 run. The cluster first called this
   (cluster 10) is proepicardial — see §6.4
2b. **Proepicardial** (TBX18⁺/WT1⁺/TCF21⁺/TBX5⁺/LHX2⁺/SFRP5⁺, GATA4-high, NKX2-5 low) — 50% D11,
   persists to D30; sits at the tip of the stromal branch
3. **Cardiomyocyte (atrial) vs (ventricular) vs (immature)** — the figures lump these; the
   data says atrial/immature with essentially no ventricular
4. **Epicardial** (WT1⁺/ALDH1A2⁺) — rising to D30
5. **Hepatic/Endoderm** (AFP⁺/FOXA2⁺) — potentially large, magnitude pending soup correction
6. **Proliferating** — cycling fraction varies 5–60% across timepoints
7. **Myofibroblast vs resting fibroblast** — ACTA2/POSTN vs DCN/LUM/TCF21

`hiPSC.jpg` is an hiPSC → **NK cell** differentiation (HSC_pro, hiPSC_NK_main, PBMC-NK,
stromal). Its cell-type vocabulary does not transfer to this cardiac series; treated here as
an example of presentation style and of the `?uncertain` labelling convention only.

---

## 6. `cm` mode as implemented

`pipeline/config_species_cm.R`, sourced by `config.R` when `SCRNA_SPECIES=cm`; keyword `cm`
accepted by `run_pipeline.sh`. Kept in its own file rather than added to
`config_species_bat.R`.

| setting | value | rationale |
|---|---|---|
| `QC$max_features` | 9000 | human default 5000 discards 12.7% of cells (21.2% of D0_2/D30_2); 9000 sits just above p99 for every sample and drops loss to 0.2% |
| `QC$max_counts` | 70000 | p99 reaches 54,422; one D30_2 cell hits 435,190 |
| `QC$max_percent_mt` | **unchanged (20)** | expected to need raising for mito-rich cardiac cells; measured median 0.1–1.0%, p99 4.6%, only 13/107,658 cells above 20% — the gate is already inert. Same C4 property as the bat data |
| `SINGLER_REF` | HumanPrimaryCellAtlas | supplies Embryonic_stem_cells, iPS_cells, Fibroblasts, Endothelial_cells, Epithelial_cells, Smooth_muscle_cells, Hepatocytes, MSC |
| `CLUSTER$default_res` | 0.5 (sweep 0.1/0.3/0.5/0.8) | the source figures used res 0.1 (5–6 clusters), too coarse to separate atrial/ventricular CM, epicardium and endoderm across 107k cells |
| `CONTAMINATION_TYPES` | Endothelial, Pericyte, Pluripotent, Proliferating | small-but-real populations that cluster majority vote would swallow |
| `SUBTYPE_MARKERS` | CM → ventricular/atrial/immature; Fibroblast → activated/resting | |
| `CELLTYPE_COLORS` | +10 entries in the base palette | CM lineage reds/oranges, progenitor+pluripotent purples, endoderm brown |
| `SINGLER_NORM` | +4 HPCA entries | Embryonic_stem_cells/iPS_cells → Pluripotent, Hepatocytes → Hepatic/Endoderm, Neuroepithelial_cell → Epithelial |

### 6.1 The blocking limitation: HPCA has no cardiomyocyte label

Verified 2026-09-10 — of HPCA's 36 `label.main` values, none is a cardiomyocyte; the nearest
cardiac-adjacent label is `Smooth_muscle_cells`. **SingleR therefore cannot call the main
population of this experiment** and will assign cardiomyocytes to Smooth_muscle_cells,
Fibroblasts or MSC.

This is structurally identical to bat neutrophils under MonacoImmune, and the same three
options apply:

1. **Marker route** (in effect now): scType scores the `MARKERS` panel, then a guarded
   `CLUSTER_CELLTYPE_MAP` built by reading `annotation/canonical_markers_dotplot.pdf`. Follow
   the `setequal(SAMPLE_NAMES, …)` guard pattern — cluster numbers are not stable across runs,
   and a partial map shatters every cluster left out of it.
2. **Frozen reference** — `build_reference.R` on a curated CM run would give run-independent
   `cell_type_ref` labels including Cardiomyocyte. This is the durable fix and is what
   ultimately resolved the bat case. It requires one hand-curated run first.
3. A cardiac-specific external reference (e.g. a published fetal-heart atlas) via celldex.

Expect a high `Unassigned` rate on the first pass and do **not** read raw SingleR labels as
the answer.

### 6.2 Harmony over a differentiation time course — open question

The pipeline integrates on `sample`. D0 pluripotent cells and D30 cardiomyocytes are
genuinely different cell states, not batch variants of one another, so per-sample correction
risks regressing out the differentiation itself — the same concern raised for dissection depth
in `docs/bat_wing_readiness.md` §.

The first run keeps Harmony on (replicates do need batch correction) with timepoints carried
as `SCRNA_CONDITION`. **Check `integrated/harmony_before_after.pdf` and the UMAP before
trusting the composition:** if D0 and D30 collapse onto each other, re-run 04 with a reduced
`HARMONY$theta`, or integrate replicates only within timepoint. Step 04's mixing-quality
warning fires on under-correction (>80% one sample per cluster), not over-correction, so this
one has to be judged by eye.

---

## 6.3 Steps 09/10 result, and a correction to the endoderm finding (2026-09-10)

Steps 09 (bootstrap proportions) and 10 (rarefaction) were run against
`Results/results_H1D01_7samples_filtered`. Outputs: `bootstrap_proportions_report.pdf`,
`bootstrap_summary.csv`, `rarefaction_report.pdf`, `rarefaction_summary.csv`.

**Step 09 says the proportions are statistically precise.** Bootstrapping all samples
to the smallest (9,168 cells, 1,000 draws) gives CI widths under 1.5 pp throughout, and
replicates agree closely — D20 cardiomyocyte 1.98% [1.79–2.20] vs 2.12% [1.93–2.31],
D30 4.49% vs 5.59%. The differences between timepoints are far larger than the CIs, so
none of the composition shifts are sampling noise.

**But step 09 does NOT settle the D20 depth question, and an earlier note in this file
implied it would.** Step 09 resamples *cells* to a common count; it corrects for
differing cell numbers and sampling noise, not for differing library depth. Depth bias
acts earlier — on which genes are detected and therefore which label a cell receives —
and resampling cells cannot undo it. The test that does address it is depth
stratification within a timepoint:

| timepoint | median genes per quartile | CM lineage % | Hepatic/Endoderm % |
|---|---|---|---|
| D11 | 1,422 / 2,431 / 3,286 / 4,311 | 25.2 → 19.7 → 16.1 → 11.5 | 25.6 → 35.0 → 30.1 → 23.2 |
| D20 | 839 / 1,123 / 1,337 / 3,771 | 1.4 → 1.8 → 4.5 → 4.0 | **75.4 → 71.1 → 54.7 → 2.9** |
| D30 | 1,188 / 2,105 / 3,540 / 5,053 | 9.4 → 14.5 → 5.5 → 2.3 | **65.3 → 18.9 → 5.4 → 2.4** |

The endoderm fraction collapses from ~70% in the lowest-depth quartile to ~3% in the
highest, at both D20 and D30. A genuine cell type does not vanish that completely in the
deepest-sequenced cells.

**Per-cluster depth confirms it.** Cluster 0 — 23,498 cells, 24.8% of the run, the
single largest cluster and the basis of the "Hepatic/Endoderm 29.2%" figure — has a
median of **1,141 genes / 2,430 UMI**, the lowest of any substantial cluster. The other
endoderm clusters sit at normal depth: cluster 11 at 2,784 genes, cluster 13 at 2,880
(11,340 UMI), cluster 20 at 3,104.

Cluster 0's top `FindAllMarkers` genes are RBP4, TTR, FGB, AHSG, APOC3, APOA1, APOA2,
AFP, CST3 — **all secreted, highly abundant plasma proteins** — plus **MALAT1**, a
ubiquitous nuclear lncRNA that is a classic signature of ambient-dominated or
low-quality droplets. In a low-content cell, ambient transcripts make up a larger
*share* of the transcriptome, so soup genes appear "enriched" relative to high-content
clusters. Its median %MT is only 0.55, so these are not dying by the mitochondrial
criterion — they are simply low-content.

**Revised reading of the endoderm result:**

| | cells | % of run | confidence |
|---|---|---|---|
| Definitive hepatic endoderm (clusters 11, 13, 20) | 4,193 | **4.4%** | high — normal depth, definitive programmes: ALB, APOB, MTTP, F2, CEBPA, HNF4A, HNF1A, ONECUT1, HHEX, FOXA2 |
| Cluster 0 | 23,498 | 24.8% | **low — ambiguous.** Low-content cells whose profile matches the ambient plasma-protein soup. Either degraded/ruptured cells that passed the 200-gene floor, or genuinely low-RNA endoderm |

So **"hepatic endoderm is 29.2% and the dominant product" is not supported.** The
defensible statement is: definitive hepatic endoderm is ~4.4%, and a further ~25% of
barcodes are low-content and unresolvable without ambient correction. The cardiac
numbers are less affected — cluster 10 sits at 1,561 genes and cluster 12 at 3,038, and
its markers (MYH6, TTN, MYOCD, SLC8A1, LDB3, TECRL) are not abundant secreted
transcripts.

**What would resolve cluster 0**, in increasing order of effort:

1. Raise `QC$min_features` from 200 to ~800–1,000 for this cohort and re-run. Cheap, and
   would show immediately how much of the composition depends on those barcodes.
2. Run ambient correction — SoupX or CellBender — on the raw (unfiltered) matrices, which
   is the principled fix. Neither is currently in the pipeline or `setup_env.sh`.
3. Ask whether the D20 preparation had a viability or over-loading problem (question 2 in
   `docs/cm_questions_for_data_owner.md`); 70% of cluster 0 is D20.

Until then, quote the cardiac and pluripotent numbers, and treat the endoderm fraction
as a range (4.4% definitive, up to ~29% if cluster 0 is real).

---

## 6.4 min_features 1000 re-run, relabels, and trajectory (2026-09-14)

Option 1 above was taken: `QC$min_features` 1000 (in `config_species_cm.R`). The run now has
82,399 cells in 20 clusters; the 200-floor run is preserved as
`Results/results_H1D01_7samples_filtered_minfeat200/`. Cluster numbers differ between the two,
so the §6.3 cluster numbers refer to the old run only.

**Curation changes in the 1000 run**

- **Cluster 0 is still there** (17,581 cells, culture-average soup: cargo genes without identity
  TFs) and is labelled `Unknown`. Raising the floor did not dissolve it; SoupX on raw matrices
  remains the fix.
- **Cluster 1 was mixed** and is split with the new `CLUSTER_SUBCLUSTER_MAP` hook in
  `05_annotate.R` (FindSubCluster at res 0.2): 3,530 POU5F1-high cells (87% D0) → Pluripotent,
  the rest → Fibroblast. D0 went from 76.2% to **96.1% pluripotent** (replicates 95.9 / 96.3).
- **Cluster 10 is Proepicardial, not Cardiac progenitor.** % of cells expressing:

| | TBX18 | WT1 | TCF21 | TBX5 | GATA4 | PDGFRA | NKX2-5 |
|---|---|---|---|---|---|---|---|
| cluster 10 | 17 | 25 | 27 | 38 | 88 | 49 | **13** |
| cardiomyocytes (cl 7) | 3 | 4 | 4 | 35 | 63 | 19 | 38 |

  Top markers are C7, SCN7A, LHX2, SFRP5, TBX18, HGF, COLEC11 — a proepicardial /
  epicardium-precursor programme. 50% of it is D11 but it persists to D30.

**Trajectory** — `pipeline/projects/cm/trajectory.R` (monocle3), outputs in `trajectory/`.

- `trajectory.R cardiac` (24,435 cells; Pluripotent, Cardiomyocyte, Proepicardial, Epicardial,
  Fibroblast, Myofibroblast; each type capped at 6,000): rooted on D0 pluripotent cells, it gives
  **two fates** — a cardiomyocyte branch, and a stromal branch that splits into epicardial and
  proepicardial. Median scaled pseudotime per sample rises with day and replicates agree
  (D0 0.15/0.15, D11 0.64, D20 0.64/0.71, D30 0.72/0.72). Pseudotime is distance from the root
  along the tree, so values on different branches are not comparable.
- `trajectory.R cm` (Pluripotent + Cardiomyocyte, 9,310 cells): the two populations are separate
  islands joined by one forced graph edge — **no cells bridge D0 and D11**, because nothing
  between those days was sampled (MESP1 is flat throughout). Early pseudotime is not meaningful.

**Pseudotime is not a maturation axis here; the MYH7 fraction is.** Within cardiomyocytes, D11
cells get the *highest* pseudotime (median 36.9 vs 32.6 at D20 and 34.5 at D30). D11
cardiomyocytes are sequenced ~2× deeper (median 5,752 UMI vs 2,512–3,502), so they show more of
every gene and the graph orders them by expression amplitude (pseudotime vs genes detected,
Spearman −0.28). Within-cell isoform ratios cancel depth:

| per sample | D11 | D20_1 | D20_2 | D30_1 | D30_2 |
|---|---|---|---|---|---|
| MYH7 / (MYH6 + MYH7) | 0.01 | 0.25 | 0.25 | 0.39 | 0.41 |
| TNNI3 / (TNNI1 + TNNI3) | 0.06 | 0.04 | 0.04 | 0.04 | 0.03 |
| MYL2 / (MYL2 + MYL7) | 0.00 | 0.00 | 0.00 | 0.00 | 0.01 |

- **MYH6 → MYH7 switch** is clean, stepwise, and replicate-tight: genuine maturation D11 → D30.
- **No TNNI1 → TNNI3 switch** — the cardiomyocytes are still fetal-like at D30.
- **MYL2 absent** — not committed ventricular cardiomyocytes, consistent with §3.

Use the cardiac tree for lineage structure and the MYH7 fraction for maturation. The Moran's I
gene lists (`pseudotime_genes*.csv`, ~15–16k genes at q < 0.05) are dominated by ribosomal
genes, MALAT1 and lncRNAs; filter those and rank by `morans_I` before interpreting.

---

## 7. Run

```bash
bash pipeline/run_pipeline.sh cm \
  Samples/Cardio/H1D0_1 Samples/Cardio/H1D0_2 Samples/Cardio/H1D11_2 \
  Samples/Cardio/H1D20_1 Samples/Cardio/H1D20_2 Samples/Cardio/H1D30_1 \
  Samples/Cardio/H1D30_2 \
  condition=H1D0_1=D0,H1D0_2=D0,H1D11_2=D11,H1D20_1=D20,H1D20_2=D20,H1D30_1=D30,H1D30_2=D30
```

Output: `Results/results_H1D01_7samples_filtered/` — 6 PDFs plus the interactive HTML report,
the same deliverable shape as `results_ES01_8samples_filtered`.

Steps 09 (bootstrap proportions) and 10 (rarefaction) matter more than usual here because of
§2.1 and should be run against the finished run dir. Trajectory (§6.4) runs against the
finished run dir too:

```bash
SCRNA_SPECIES=cm SCRNA_RESULTS_DIR=Results/results_H1D01_7samples_filtered \
  Rscript pipeline/projects/cm/trajectory.R cardiac   # or: cm
```
