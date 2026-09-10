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
2. **Cardiac progenitor** (NKX2-5⁺/ISL1⁺/TBX5⁺) — a distinct transient D11 state
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
§2.1 and should be run against the finished run dir. A trajectory analysis is genuinely
appropriate for this design, unlike the blood cohorts, but steps 11–14 currently live under
`pipeline/projects/bat_wing/` with wing-specific labels and would need a `projects/cm/`
sibling.
