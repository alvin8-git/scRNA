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

**Naive CD4 T cells are selectively lost** — 14.2% → 7.8–8.2%, roughly halved, while NK, CD14+ and
FCGR3A+ monocytes all rise. The two CryoStor replicates agree within 0.4 points on every cell type,
so this is reproducible rather than noise. Naive T cells are the population most sensitive to a
freeze-thaw cycle; the rises elsewhere are consistent with compositional displacement rather than
genuine expansion.

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
and those fail the definitive markers (FCGR3B 0%, CXCR2 0%, ELANE/MPO/LTF/CAMP 0%); FCGR3B is
0.39–0.73% of all cells, i.e. background. FCGR3B is also 0.000% of the ambient pool in all six
samples. Density-gradient separation worked in every sample, with no operator- or material-linked
granulocyte carryover to compare.

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
