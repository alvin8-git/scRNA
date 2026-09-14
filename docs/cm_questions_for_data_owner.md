# Questions for the data owner — H1 cardiomyocyte differentiation time course

Ready to send. Each question states the observation behind it, so the answer can be
specific rather than general. Observations come from the 7-sample H1 series in
`Samples/Cardio` (107,658 cells, D0/D11/D20/D30); detail in
`docs/cm_differentiation_readiness.md`.

---

## 1. Which cardiomyocyte subtype was the target — ventricular, atrial, or nodal?

We ask because the data looks **immature and atrial-like, with essentially no ventricular
specification**:

| marker | identity | detection at D30 |
|---|---|---|
| MYL2 | ventricular | **1.4%** |
| IRX4 | ventricular | 5.0% |
| HEY2 | ventricular | 8.3% |
| MYL7 | immature/atrial | **89.8%** |
| TNNI1 | immature | 72.5% |
| NR2F2 | atrial | 23.6% |
| NPPA | atrial | 11.4% |

If ventricular cardiomyocytes were the goal, this is the headline result rather than a
detail. If the target was atrial or unspecified/generic CM, the outcome looks as expected.
Either way it changes how we frame the report.

Related: is D30 the intended endpoint, or an interim one? The reference figures you shared
are Day 21 and Day 42.

## 2. Was anything different about D20 — loading, dissociation, or library prep?

Both D20 samples stand out, but **not simply as "shallower sequencing"**. They captured the
*most* cells yet have by far the lowest per-cell content, and their UMI distribution is
heavily skewed:

| sample | cells | total UMI | median UMI | median genes | mean/median UMI |
|---|---|---|---|---|---|
| H1D11_2 | 17,109 | 202M | 10,143 | 3,208 | 1.17 |
| H1D30_2 | 15,159 | 199M | 11,351 | 3,508 | 1.16 |
| H1D30_1 | 16,715 | 152M | 7,531 | 2,822 | 1.21 |
| **H1D20_1** | **19,279** | 136M | **2,977** | **1,300** | **2.37** |
| **H1D20_2** | **17,834** | 129M | **2,139** | **1,019** | **3.39** |

D20_1 has only 13% more cells than D11 but 3.4× lower median UMI, so cell number alone does
not account for it. A mean/median ratio of 2.4–3.4 (vs ~1.2 elsewhere) means a large
population of **low-content barcodes** — debris, dying cells, or an over-loaded chip — rather
than an even reduction in depth.

Specifically:
- Were the D20 libraries prepared or sequenced in a **separate run / different kit lot**?
- Was **more cell suspension loaded** for D20?
- Was **viability or dissociation** noticeably worse at that timepoint? (D20 is often when
  the monolayer is densest and hardest to dissociate.)

This matters because **every cardiomyocyte marker appears to dip at D20 and recover at D30**
(TNNT2 30.1 → 15.2 → 37.5; TTN 44.6 → 21.4 → 65.2). We suspect that dip is largely technical,
and we do not want to report a differentiation "regression" that is really a capture artefact.

## 3. What are your definitions for the non-cardiomyocyte types?

Not asking for marker lists — those are in the literature, and we have a validated panel
(appendix below). We are asking about **your conventions**, because these types share markers
and your prior figures make specific choices we would like to match:

- Your Day-21/42 panels group **ANPEP + PDGFRA + PDGFRB** as *Pericytes*. PDGFRA is more
  commonly a fibroblast marker — is that grouping deliberate for this system?
- Your Day-42 figure has a cluster labelled **"Fibroblast/Epi"**. Is that a category you want
  preserved as an honest hedge, or resolved into one or the other?
- Do you count **ISL1⁺/NKX2-5⁺** cells as *cardiac progenitor*, *second heart field*, or fold
  them into cardiomyocytes?
- Do you distinguish **myofibroblast** from *fibroblast*, or report them together?

The aim is that our labels line up with what you have already published, so the two are
comparable.

## 4. Do you have the processed object with cluster assignments from the Day-21 / Day-42 work?

A table of `barcode → cell_type` (or the Seurat/Scanpy object) would be the single most
useful thing you could send.

Reason: the standard reference we classify against (HumanPrimaryCellAtlas) has **no
cardiomyocyte label at all** — of its 36 cell types the nearest cardiac one is
"Smooth_muscle_cells". So automated annotation cannot name the main population of your
experiment, and currently assigns cardiomyocytes to smooth muscle, fibroblast or MSC. With
labelled cells from a run you have already curated, we can train a cardiac-specific reference
and get consistent labels across every future run.

If none exists, we will curate one pass by hand from the marker dot plot and build the
reference from that — it just costs a round of manual annotation.

## 5. Is a hepatic / endoderm population expected from this protocol?

**AFP** is detected in 6.4% of cells at D0, 12.8% at D11, then **95.6% at D20 and 99.4% at
D30**. FOXA2 peaks at 23.4% at D11.

Two possible readings, and you may know immediately which:
- a genuine yolk-sac / hepatic endoderm off-target, a known failure mode of some Wnt-modulation
  protocols, or
- ambient RNA from a smaller number of such cells (COL1A1, COL3A1, LUM and DCN are likewise at
  99–100%, which is a soup signature rather than biology).

Have you seen AFP⁺ cells by flow or immunostaining in these differentiations?

## 6. Was any cardiomyocyte purification attempted?

The data shows no sign of it — fibroblast, epicardial and endoderm populations are all
present. Was **lactate selection** or **SIRPA⁺/CD90⁻ sorting** used and ineffective, or
deliberately omitted so the whole culture could be profiled? This determines whether the
substantial non-cardiomyocyte fraction is a problem to report or the point of the experiment.

## 7. Experimental design

- Are the replicates **independent differentiations** (biological) or **split wells from one
  differentiation** (technical)? This decides which statistical comparisons are legitimate.
- **Why does D11 have no replicate** — lost, failed QC, or never planned? D11 is the only
  window on the cardiac progenitor state (NKX2-5 18.5%, TBX5 19.7%, ISL1 15.5%, all gone by
  D20), so we cannot separate biology from sample idiosyncrasy there.
- Why these timepoints (D0/D11/D20/D30)? D11 lands neatly on the progenitor window — was that
  deliberate?
- Are these the **same series** as the Day-21/Day-42 figures, or a separate experiment? Those
  figures show clear endothelial and pericyte clusters; here endothelium is ~1% (PECAM1 1.1%,
  CDH5 1.4%).

## 8. Protocol and literature

Could you point us to:
- the **differentiation protocol** you followed (or your own paper/preprint describing it) —
  small-molecule Wnt modulation (CHIR/IWP), monolayer vs EB, media schedule;
- any **recent papers you consider the benchmark** for this system, including your lab's own
  work;
- any **orthogonal data** on these same samples — flow, immunofluorescence, bulk RNA-seq,
  electrophysiology or contractility. Anchoring even one timepoint would resolve several
  ambiguities at once.

## 9. What decision does this dataset support?

Protocol development, a baseline for a later perturbation (disease line, drug, gene edit), or
a standalone characterisation? And do you want **pseudotime/trajectory** analysis, or
per-timepoint composition? The design supports trajectory properly — unlike most of our other
datasets — but it is worth building only if it answers a question you actually have.

---

## Appendix — the marker definitions we are currently using

So you can diff against yours. Detection rates are from your data (D0 → D11 → D20 → D30).

| type | markers | note from this data |
|---|---|---|
| Cardiomyocyte | TNNT2, TTN, ACTC1, ACTN2, **MYL7**, MYL4, **TNNI1**, MYH6, DES, MYBPC3 | MYL7 (89.8%) and TNNI1 (72.5%) discriminate better than MYL4 (54.5%) |
| CM (ventricular) | MYL2, IRX4, MYH7, HEY2 | included so the near-absence is visible |
| CM (atrial) | NPPA, NR2F2, SLN, GJA5, KCNJ3 | SLN peaks at D20 (34.9%) then falls to 4.7% |
| Cardiac progenitor | NKX2-5, ISL1, TBX5, GATA4, MESP1, KDR | D11-specific detection, but no distinct progenitor cluster in the min_features 1000 run |
| Proepicardial | TBX18, WT1, TCF21, TBX5, LHX2, SFRP5, GATA4 | NKX2-5 only 13%; 50% D11, persists to D30 (added 2026-09-14) |
| Pluripotent | POU5F1, SOX2, NANOG, LIN28A, DNMT3B | POU5F1 (89.2% at D0) far more sensitive than NANOG (27.4%) |
| Fibroblast | COL1A1, COL3A1, DCN, LUM, POSTN, TCF21, CCDC80 | ⚠️ COL1A1/COL3A1/LUM/DCN all at 99–100% detection — soup-compromised; POSTN and TCF21 are the reliable ones |
| Myofibroblast | ACTA2, TAGLN, FN1, POSTN | |
| Smooth muscle | MYH11, CNN1, TAGLN, ACTA2 | MYH11 only 4.1% — little true SMC |
| Pericyte | **RGS5, KCNJ8, NOTCH3**, PDGFRB, ANPEP | PDGFRB (42%) alone cannot separate pericyte from fibroblast here |
| Endothelial | PECAM1, CDH5, APLNR, ENG | small: PECAM1 1.1%, CDH5 1.4% at D30 |
| Epithelial | EPCAM, KRT8, KRT18, KRT19 | EPCAM falls 51.2 → 17.5% — most signal is pluripotent cells, not late epithelium |
| Epicardial | WT1, TBX18, ALDH1A2, UPK3B | absent from your figures; WT1 rises 0.1 → 18.4% |
| Hepatic/Endoderm | AFP, FOXA2, SOX17, TTR, APOA1 | see Q5 |
| Proliferating | MKI67, TOP2A, CCNB1 | TOP2A 60.8 → 14.2 → 33.4% |

**Two markers we would flag:**
- **CCDC80** was suggested as a cardiac marker. In this data it tracks the fibroblast
  compartment (4.9 → 61.5 → 28.9 → 68.3), so we have grouped it there.
- **TNNI3** behaves anomalously — 38.5% in undifferentiated D0 cells, *falling* to 12.1% by
  D30. It should rise with maturation and should not be present in hESCs at all. We are not
  using it as a maturation marker; is it reliable in your hands?

**Checked and absent** (so the usual iPSC off-targets are not a concern here): MYOD1/MYOG/PAX7
(skeletal muscle), PAX6/SOX10 (neural crest/neuroectoderm), RAX/SIX6/VSX2 (retina). The
off-target in this series is endodermal only.
