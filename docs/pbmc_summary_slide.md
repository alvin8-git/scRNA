# PBMC loading standards — summary of findings

- **Loading does not distort composition.** The four fresh samples span a 6-fold loading range and agree within ~3 points on every cell type (CD4 memory 22.8–24.0%, naive B 16.7–18.6%, CD14+ Mono 23.1–26.6%).
- **Freeze-thaw costs library complexity, not capture.** CryoStor (65.3% viability) vs fresh (82.5%) at 30K: 16–21% fewer genes, 24–29% fewer UMIs, saturation +11–14 pts — but recovery unchanged (43.0–48.9% vs 44.8%).
- **Naive CD4 T cells roughly halve** with freeze-thaw (14.2% → 7.8–8.2%); NK, CD14+ and FCGR3A+ monocytes rise correspondingly. Both CryoStor replicates agree within 0.4 pts.
- **Doublets scale with loading** (4.30% → 12.57%). The CryoStor pair runs ~47% above fresh at matched 30K loading — genuine clumping and/or degraded cells called as doublets; not separable here.
- **Ambient RNA is low and loading-dependent** (SoupX rho 0.016 → 0.043, vs 0.125–0.276 in bat whole blood) and is NOT raised by freeze-thaw.
- **Contamination <1.3% in every sample.** Neutrophils effectively absent: 21 of 49,348 cells, and those lack FCGR3B/CXCR2/ELANE entirely; FCGR3B is 0.000% of the ambient pool in all six.
- **Operator effects cannot be measured in this design** — operators are interleaved with loading, and confounded with material at 30K.

| Sample | Loaded | Doublet % | SoupX rho | Contamination % |
|---|---|---|---|---|
| PBMC_5K_A_QL | 5,000 | 4.30 | 0.016 | 0.94 |
| PBMC_10K_K | 10,000 | 4.37 | 0.021 | 1.04 |
| PBMC_15K_A_QL | 15,000 | 5.02 | 0.032 | 1.24 |
| PBMC_30K_K | 30,000 | 8.52 | 0.043 | 1.19 |
| PBMC_30K_QL (cryo) | 30,000 | 12.54 | 0.036 | 1.19 |
| PBMC_30K_S_K (cryo) | 30,000 | 12.57 | 0.041 | 1.08 |
