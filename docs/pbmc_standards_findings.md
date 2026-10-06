# PBMC loading-series standards — what the data supports

Run `Results/results_PBMC5KAQL_6samples_filtered` (49,348 cells). Six libraries, one flowcell
(E250126075_L01, dnbc4tools 3.1), differing only in barcode index — **sequencing batch is not a
confounder here**. Numbers below come from `qc/cell_fate.csv`, `Mirxes/Mirxes_PBMC_metrics.xlsx`,
`integrated/integrated_annotated.rds` and `/data/alvin/tmp/soupx_PBMC_rho.csv`.

| Sample | Material | Operators | Loaded | Recovered | % of loaded | Final cells |
|---|---|---|---|---|---|---|
| PBMC_5K_A_QL | StemCellTech | Alex+QiuLing | 5,000 | 2,985 | 59.7 | 2,760 |
| PBMC_10K_K | StemCellTech | Kumar | 10,000 | 4,568 | 45.7 | 4,222 |
| PBMC_15K_A_QL | StemCellTech | Alex+QiuLing | 15,000 | 7,909 | 52.7 | 7,080 |
| PBMC_30K_K | StemCellTech | Kumar | 30,000 | 13,436 | 44.8 | 11,666 |
| PBMC_30K_QL | CryoStor 1 cycle | QiuLing | 30,000 | 14,658 | 48.9 | 12,541 |
| PBMC_30K_S_K | CryoStor 1 cycle | Sabrina+Kumar | 30,000 | 12,893 | 43.0 | 11,079 |

Viability at load: **82.5%** (StemCellTech) vs **65.3%** (CryoStor), from the lab record.

---

## 1a. Fresh (StemCellTech) vs CryoStor freeze-thaw

At matched 30K loading, CryoStor material differs consistently from fresh:

| Measure | 30K_K (fresh) | 30K_QL (cryo) | 30K_S_K (cryo) |
|---|---|---|---|
| Median genes/cell | 2,588 | 2,172 | 2,043 |
| Median UMI/cell | 6,987 | 5,330 | 4,994 |
| Median %MT | 0.50 | 0.68 | 0.75 |
| Sequencing saturation | 70.9% | 81.8% | 84.7% |
| Doublet rate | 8.52% | 12.54% | 12.57% |
| Ambient rho | 0.043 | 0.036 | 0.041 |
| CD4 T (naive) | 14.2% | 8.2% | 7.8% |
| NK | 10.3% | 12.1% | 12.8% |
| FCGR3A+ Mono | 4.8% | 7.2% | 7.0% |

**Freeze-thaw costs library complexity, not capture.** Recovery is unchanged (43.0–48.9% vs 44.8%),
but each cell yields ~16–21% fewer genes and ~24–29% fewer UMIs, with saturation 11–14 points
higher at comparable mean reads per cell (69–77k across all six). Fewer distinct molecules per cell
is the signature of degraded input, not of under-sequencing.

**CryoStor did not recover more cells, and the higher saturation is an effect rather than a cause.**
Only one of the two CryoStor libraries exceeds fresh: 30K_QL recovers 12,541 final cells (+7.5%)
while 30K_S_K recovers 11,079 (−5.0%). The gap between the two CryoStor replicates is larger than
either's gap to fresh, so there is no material-level effect on cell number. Nor were they sequenced
more: total cDNA reads are within 3% across all three (994M–1,026M) and mean reads per cell is
*lowest* in 30K_QL (69,553 vs 76,365 fresh). Saturation is higher (81.8/84.7% vs 70.9%) because each
cell holds fewer distinct molecules (median UMI 5,509–5,858 vs 7,758) — the same reads resample the
same library. One asymmetry is worth noting: the fresh sample loses 4.7% of cells at QC against
1.5–2.0% for CryoStor, plausibly because the upper thresholds (`QC$max_counts` 25,000) clip more of
its high-UMI cells. [uncertain — not broken down by which QC criterion fired.]

**Naive CD4 T cells are selectively lost, and nothing else actually gains.** All three 30K libraries
loaded the same number of cells, so absolute counts are comparable:

| Population | Fresh 30K_K | Cryo 30K_QL | Cryo 30K_S_K | ratio vs fresh |
|---|---|---|---|---|
| **CD4 T (naive)** | 1,661 | 1,034 | 863 | **0.62 / 0.52** |
| CD14+ Mono | 3,105 | 3,657 | 3,302 | 1.18 / 1.06 |
| CD4 T (memory) | 2,657 | 2,893 | 2,501 | 1.09 / 0.94 |
| B cell (naive) | 2,016 | 2,001 | 1,795 | 0.99 / 0.89 |
| NK | 1,199 | 1,518 | 1,420 | 1.27 / 1.18 |
| FCGR3A+ Mono | 562 | 898 | 773 | 1.60 / 1.38 |
| Total cells | 11,666 | 12,541 | 11,079 | 1.08 / 0.95 |

Naive CD4 T is the only population that falls — ~45% fewer cells. Everything else is flat within the
spread of the two replicates.

**The monocyte "increase" is displacement, not expansion.** CD14+ monocytes track the total cell
count (×1.18 and ×1.06 against totals of ×1.08 and ×0.95); normalised, that is ×1.10–1.12, which is
exactly the 26.6% → 29.2/29.8% share shift. Proportions are compositional and must sum to 100%, so
losing ~700 naive T cells from ~11,700 inflates every surviving population's share by ~6%. The same
arithmetic explains the NK, B memory and CD4 memory wobbles (all ≤0.5 points). Mechanistically the
loss is one-directional: naive T cells are small, with a high nuclear-to-cytoplasmic ratio and little
anti-apoptotic reserve (BCL2 family, heat-shock proteins), so they take the brunt of ice formation
and osmotic shock and continue apoptosing after thaw; monocytes simply survive it.

**One change is NOT displacement: FCGR3A+ monocytes, ×1.45 after normalising.** Total monocytes
(CD14+ plus FCGR3A+) stay within ±10% of fresh, so the compartment is not growing — the subsets
redistribute. Two readings fit and this data cannot separate them: a phenotype shift (freeze-thaw is
an activating stress; monocytes upregulate CD16/FCGR3A and downregulate CD14 on activation, which
would relabel the same cells), or preferential survival of the non-classical subset. The first is
more parsimonious. [uncertain — needs surface protein or an activation score.]

**Ambient RNA does not increase** (0.036–0.041 vs 0.043). This is worth noting because it is the
opposite of the intuitive expectation: freeze-thaw damage shows up as degraded cells and doublet
calls, not as free-floating mRNA.

Operator cannot be blamed for any of the above: both CryoStor libraries were run by different
operators (QiuLing; Sabrina+Kumar) and agree closely with each other, while differing from Kumar's
fresh 30K sample.

## 1b. Comparison with prior human PBMC runs

| Comparator | Cells | Valid benchmark? |
|---|---|---|
| H1 / H2 (`results_H1-H2_filtered`) | 359 / 848 | **No.** Bundled quick-start fixture. At n=359, RBC reads 14.2% and NK 29.8% — sampling noise. Also the only samples with hardcoded doublet rates (`DOUBLET$doublet_rate` H1=0.031, H2=0.077), so their doublet numbers were set by config, not estimated. |
| H1_pre / H2_post / H3_post | 2,211 / 22,489 / 9,267 | **No.** A pre- vs post-SORT experiment; sorting reshapes composition by design (CD4 memory 36.6–42.1%, CD14+ Mono 10.2–13.2%). |
| DemoScRNA | 8,255 | **Partially.** The only plain PBMC comparator. |

Against DemoScRNA, the fresh samples agree on monocytes (21.3% vs 23.1–26.6%), CD4 memory
(27.1% vs 22.8–24.0%) and NK (9.7% vs 9.2–10.9%), but **naive CD4 T differs by 14 points**
(29.2% vs 14.2–15.7%) and naive B by 4–6. A naive/memory split of that size is ordinary donor
variation (it tracks donor age), so this is most readily read as a different donor, not a platform
discrepancy. [uncertain — no donor metadata was available to confirm.]

**The new libraries are materially deeper than anything in the archive**: median 2,043–2,823 genes
and 4,994–8,013 UMIs, against 1,385–1,945 genes and 3,238–4,478 UMIs in every prior human run.
Median %MT (0.45–0.75) sits in the same low band as the prior runs (0.47–1.13), consistent with the
DNBelab C4 platform behaving normally.

## 1c. Doublets

Rate = doublets removed / cells after QC.

| Loaded | Sample | After QC | Doublets | Rate | ~0.8% per 1,000 recovered |
|---|---|---|---|---|---|
| 5,000 | 5K_A_QL | 2,884 | 124 | 4.30% | 2.3% |
| 10,000 | 10K_K | 4,415 | 193 | 4.37% | 3.5% |
| 15,000 | 15K_A_QL | 7,454 | 374 | 5.02% | 6.0% |
| 30,000 | 30K_K | 12,753 | 1,087 | 8.52% | 10.2% |
| 30,000 | 30K_QL | 14,339 | 1,798 | 12.54% | 11.5% |
| 30,000 | 30K_S_K | 12,672 | 1,593 | 12.57% | 10.1% |

Rates scale monotonically with loading and sit close to the rule-of-thumb expectation at high load,
running above it at low load. Prior runs for context: DemoScRNA 11.79%, H1_pre 5.35%, H2_post
19.03%, H3_post 11.66% — the new series sits inside that envelope and, unlike those, scales cleanly.
**Any doublet comparison against the H1/H2 fixture is invalid** (hardcoded rates). These six used
scDblFinder's own estimate, since `DOUBLET$doublet_rate[[sample]]` returns NULL for unlisted samples.

The ~47% excess in the CryoStor pair has two readings: genuine multiplets from clumping of
freeze-thaw-damaged cells, or scDblFinder calling degraded low-complexity cells as doublets. Their
lower genes/UMI favours a contribution from the second. The data cannot separate them; genetic or
hashtag demultiplexing would.

## 1d. scDblFinder (this pipeline) vs scrublet (dnbc4tools)

The dnbc4tools logs (`<sample>/logs/<date>.txt` on the mount) report a doublet rate per sample.
It is **not** comparable to ours as printed: scrublet's `Estimated` is an *inferred total* doublet
burden, while our figure is the proportion actually removed.

| Sample | scrublet flagged | scrublet "Estimated" | ours flagged | ours ÷ scrublet flagged |
|---|---|---|---|---|
| PBMC_5K_A_QL | 25/2,760 = 0.91% | 2.2% | 124/2,884 = 4.30% | 4.7× |
| PBMC_10K_K | 70/4,193 = 1.67% | 3.6% | 193/4,415 = 4.37% | 2.6× |
| PBMC_15K_A_QL | 178/7,302 = 2.44% | 5.2% | 374/7,454 = 5.02% | 2.1× |
| PBMC_30K_K | 490/12,392 = 3.95% | 7.8% | 1,087/12,753 = 8.52% | 2.2× |
| PBMC_30K_QL | 887/13,541 = 6.55% | 10.8% | 1,798/14,339 = 12.54% | 1.9× |
| PBMC_30K_S_K | 755/11,873 = 6.36% | 10.8% | 1,593/12,672 = 12.57% | 2.0× |

**Denominators are not the explanation.** scrublet ran on dnbc4tools' own cell calls
(2,760–13,541), ours on post-QC cells (2,884–14,339) — within ~5%.

**scrublet's "Estimated" is flagged ÷ detectable.** The log records
`Estimated detectable doublet fraction` = 41.1% (5K) rising to 60.7% (30K_QL); 0.91/0.411 = 2.2%.
scrublet detects only *heterotypic* doublets and scales its flagged count up to guess the total
including homotypic ones it cannot see. scDblFinder performs no such inflation — its calls are its
estimate. So the like-for-like comparison at 5K is 4.30% vs **0.91%**, not vs 2.2%.

**The priors differ, and that drives the pattern.** scrublet logs `Expected = 5.0%` for all six
samples — flat, regardless of loading — with an auto threshold falling from 0.32 to 0.16 as loading
rises. Our call is `scDblFinder(dbr = NULL)` (`02_doublets.R:42`), so v1.23.4 derives the
expectation from `dbr.per1k = 0.008`, i.e. 0.8% per 1,000 cells, which scales with cell count
(2.3% at 5K to 11.5% at 30K). The excess is roughly proportional (~2×) except at 5K (4.7×), where
scrublet's flat 5% prior is least suited to a small library.

**Both tools agree on the biology.** Each rate rises monotonically with loading, and the CryoStor
excess at matched 30K loading is 1.38× by scrublet (10.8 vs 7.8) and 1.47× by ours (12.55 vs 8.52).
Two independent methods, same conclusion.

**Which to quote.** Use the pipeline figure when describing what the analysed data contains — those
cells are genuinely absent from the 49,348. Use scrublet's `Estimated` only as an independent guess
at true doublet burden, stating that it is inferred rather than observed. Nothing elsewhere in this
document changes under either tool.

## 1e. Does ambient correction change the cell types detected? (tested, not assumed)

All six samples were SoupX-corrected (`adjustCounts`, counts removed 1.60–4.30%, matching rho
exactly, no cell lost to the correction) and re-run as a separate cohort,
`Results/results_PBMC5KAQLsx_6samples_filtered`. Corrected vs original, overall:

| Cell type | Original % | Corrected % | Δ |
|---|---|---|---|
| **CD14+ Mono** | 27.57 | 22.66 | **−4.91** |
| **RBC** | 0.25 | 3.79 | **+3.54** |
| CD4 T (naive) | 11.58 | 11.99 | +0.41 |
| B cell (naive) | 16.72 | 17.09 | +0.37 |
| NK | 11.27 | 11.63 | +0.36 |
| CD4 T (memory) | 22.98 | 23.16 | +0.18 |
| all others | | | ≤0.10 |

**The biology does not change.** No cell type appears or disappears, and every genuine population
moves ≤0.41 points — which is what rho 0.016–0.043 should produce. Cells drop 49,348 → 47,883 (−3%).

**But the CD14+ Mono → RBC swap is an artefact, and it is the useful finding.** The corrected run's
1,815 "RBC" cells are not red cells: HBB is detected in **4.2%** of them (HBA1 1.4%, HBA2 1.7%,
ALAS2 1.3%), against **50.4%** in the original run's 121 genuine RBCs. They have the lowest
complexity in the run — **1,442 median genes vs 2,294** for everything else — and retain myeloid
character (LYZ 34.7%, S100A8/9 ~26–28%) with low CD3E (5.9%) and MS4A1 (3.4%). They occur in all six
samples including the cleanest (rho 0.016).

Subtracting 1.6–4.3% of counts pushed a tail of already-low-complexity monocytes below the point
where SingleR can identify them, and HumanPrimaryCellAtlas's RBC profile — low gene count, few
distinguishing transcripts — became their nearest match.

**Conclusion: do not ambient-correct at this contamination level.** It costs ~3% of cells and
mislabels ~4% of them, in exchange for sub-0.5-point changes to real populations. Contrast the bat
ES49 cohort at rho ≈ 0.2, where correction recovered 2.9 points of neutrophils and cut RBC from
10.55% to 2.04% — a real gain. Correction earns its keep when ambient is the dominant distortion;
below roughly rho 0.05 it degrades borderline cells faster than it cleans real ones. Measure rho
first, then decide.

## 2. Further conclusions

**Loading does not distort composition.** Across a 6-fold loading range on identical material, the
four fresh samples agree within ~3 points on every cell type (CD4 memory 22.8–24.0, naive B
16.7–18.6, CD14+ Mono 23.1–26.6). For a loading standard this is the key result: pick the loading
from the cell yield you need, not from fear of compositional bias.

**Ambient RNA scales with loading, nearly linearly** — rho 0.016 → 0.021 → 0.032 → 0.043 for 5K →
10K → 15K → 30K. All values are low (1.6–4.3%; the bat whole-blood cohorts measured 12.5–27.6%),
so no correction is warranted here, but the dose-response is a usable rule: higher loading buys
cells at a modest cost in soup.

**Operator effects cannot be separated from loading in this design.** Alex+QiuLing ran 5K and 15K;
Kumar ran 10K and 30K — the two operators are perfectly interleaved across the loading series, and
at 30K operator is additionally confounded with material. One suggestive observation: both
Alex+QiuLing libraries have higher median genes than the adjacent Kumar libraries (2,614 vs 2,440;
2,823 vs 2,588). That is a zigzag against a smooth loading trend, which is what an operator effect
would look like — but with n=2 per operator and loading differing in every pair, it is **not
evidence of operator bias** [uncertain].

**Neutrophils are absent, as a PBMC prep requires.** 21 of 49,348 cells (0.043%) carry the label,
and those fail the definitive markers (FCGR3B 0%, CXCR2 0%, ELANE/MPO/LTF/CAMP/CEACAM8 0%) — they
are LYZ⁺/S100A8⁺ monocyte-adjacent cells given a nearest-available label. FCGR3B is 0.39–0.73% of
all cells, i.e. background, and **0.000% of the ambient pool** in all six samples; even lysed
neutrophils would leave debris in the soup, and there is none. The four archive runs behave
identically (45 labelled cells across 42,222, FCGR3B in 2.2% of them).

Three independent reasons, all pointing the same way:

1. **The prep removes them by design.** PBMC means *mono*nuclear; neutrophils are polymorphonuclear.
   Ficoll-Paque (ρ ≈ 1.077 g/mL) works because granulocytes are denser — their granule content
   sediments them to the pellet with the erythrocytes, while lymphocytes, monocytes, DC and NK band
   at the interface. The neutrophils are in the discarded pellet.
2. **They barely survive handling.** The shortest-lived leukocyte (circulating half-life of hours),
   apoptosing and degranulating soon after draw. The two CryoStor samples carry a second filter:
   neutrophils essentially do not survive freeze-thaw.
3. **Droplet scRNA is biased against them anyway** — ~10–20× less mRNA than a lymphocyte plus high
   RNase content, so survivors tend to fall below the UMI knee at cell calling or below
   `QC$min_features` afterwards.

**The control that proves it is the prep, not the pipeline:** the same pipeline on bat whole blood
(no density separation) reports 9.3–11.9% neutrophils, up to 29.9% in one animal. When neutrophils
are present, this pipeline finds them in quantity. For a PBMC standard the near-zero count is a
positive QC signal: the gradient worked cleanly in all six preparations, with no operator- or
material-linked carryover. Granulocytes from these donors would need whole blood with RBC lysis, or
the Ficoll pellet, and ideally fresh rather than cryopreserved material.

## Limitations

- Operator is confounded with loading throughout, and with material at 30K. No operator claim is
  supportable from this design.
- One donor per material; donor and material are confounded, so "CryoStor loses naive CD4 T" is
  strictly "this CryoStor sample differs from this fresh sample".
- No fresh control at 5K/10K/15K for the CryoStor material, so the material effect is established
  only at 30K.
- Doublet calls are algorithmic, not genetically confirmed.

## Recommended follow-ups

1. **Same operator, both materials, matched loading** — settles whether the naive-CD4-T loss is the
   freeze-thaw or the donor/operator.
2. **Two operators on a split of the same tube at one loading** — the only way to measure operator
   bias here; answers the median-genes zigzag.
3. **Multiplexed (hashtag or genetic) 30K run** — separates true multiplets from degraded cells
   called as doublets, and would confirm whether the CryoStor doublet excess is real.
4. **A fresh sample at 5K and a CryoStor sample at 5K** — tests whether the loading/ambient
   dose-response is material-independent.
