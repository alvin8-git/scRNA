# Bat Wing Readiness Assessment — Samples T1–T6

Assessment date: 2026-07-27. Scope: can the existing pipeline handle `Samples/T1`–`T6` as bat wing
tissue? Evidence is from the raw matrices and the code, not from a pipeline run.

**Summary:** a `bat_wing` species mode already exists (`config_species_bat.R`, steps 11–14), so this
is tuning, not new architecture. Steps 01–07 need four small config fixes. Steps 11–14 need three
missing R packages plus condition labels. The larger issue is that **T1–T6 are not one tissue** — the
data splits them into three groups, and they should not all go into one Harmony run.

---

## 1. The six samples are three different sample types

All six are DNBelab `filter_matrix/`, 30,158 genes, identical annotation (same md5), HGNC-style
symbols in column 2, `MT-`-prefixed mito genes (13, correctly matched by `^MT-`). So the *format* is
fine and identical to the blood samples. The *content* is not.

Percent of cells expressing (≥1 UMI), from the raw matrices:

| marker | T1 | T2 | T3 | T4 | T5 | T6 |
|---|---|---|---|---|---|---|
| COL1A1 | 30.6% | 97.2% | 0.4% | 0.2% | 3.5% | 69.3% |
| DCN | 80.3% | 99.9% | 2.0% | 1.1% | 2.0% | 96.5% |
| KRT14 | 99.5% | 50.1% | 0.2% | 0.1% | 0.1% | 1.1% |
| ACTA2 | 9.5% | 43.5% | 0.3% | 1.7% | 0.5% | 43.4% |
| PTPRC (CD45) | 7.6% | 27.6% | 91.3% | 98.7% | 66.9% | 25.6% |
| CD3E | 4.0% | 25.7% | 93.5% | 89.1% | 41.7% | 23.7% |
| MS4A1 | 0.0% | 0.4% | 6.8% | 67.4% | 7.6% | 1.6% |
| S100A8 | 5.4% | 5.3% | 4.2% | 11.3% | 99.9% | 17.6% |
| MLANA | 4.9% | 3.1% | 0.1% | 0.1% | 0.0% | 0.1% |

Corroborated by share of total UMIs: T1 KRT14 1.01%, T2 DCN 1.27%, T6 DCN 0.48%, T5 HBB 6.40% and
S100A8 0.90%, while T3/T4 have COL1A1 at 0.0001% and 0.0000%.

- **T1 — wing tissue, epidermis-dominant.** KRT14/KRT10 near-universal, melanocytes present (MLANA
  ~5%).
- **T2 — wing tissue, dermis/stromal-dominant.** Fibroblast (COL1A1, DCN), smooth muscle (ACTA2),
  pericyte (RGS5 20.4%), endothelium (PECAM1 43.3%).
- **T6 — wing tissue, stromal, heavily blood-contaminated.** DCN 96.5% but HBB in 100% of cells and
  the worst mito fraction (see §2).
- **T3, T4 — not wing tissue.** Essentially zero collagen or keratin; CD45 in 91–99% of cells. Most
  plausibly **CD45-sorted immune infiltrate** (from wing or elsewhere), T4 being B-cell-rich
  (MS4A1 67%).
- **T5 — whole blood.** S100A8 in 99.9% of cells plus HBB 6.4% of UMIs: neutrophils + RBC.

**Confirm this against the experimental design before running.** The grouping matters because
`HARMONY$group_by_vars <- "sample"` corrects *per sample*. Putting wing tissue, sorted immune cells,
and whole blood into one integration would regress out the tissue difference — i.e. the biology —
as if it were batch. Run tissue samples as one cohort; treat blood/sorted fractions as a separate
run, or as `SCRNA_CONDITION` groups within a deliberately designed comparison.

## 2. QC thresholds are the blood profile and clip real tissue cells

`config_species_bat.R` does **not** touch `QC` in either bat mode (contrary to what the docs claimed
before 2026-07-27). Per-cell distributions from the raw matrices:

| S | cells | med UMI | med genes | p90 genes | p99 genes | med %MT | p90 %MT | >5000 genes | >25000 UMI | >20% MT |
|---|---|---|---|---|---|---|---|---|---|---|
| T1 | 4,865 | 7,047 | 2,154 | 4,297 | 7,215 | 0.5 | 2.0 | 6.1% | 10.4% | 0.2% |
| T2 | 7,224 | 6,513 | 2,195 | 4,098 | 6,903 | 2.2 | 6.3 | 4.5% | 5.0% | 0.8% |
| T3 | 27,138 | 4,269 | 2,130 | 4,663 | 6,411 | 0.5 | 1.8 | 7.2% | 4.2% | 0.0% |
| T4 | 21,055 | 5,587 | 2,277 | 3,795 | 5,583 | 0.6 | 1.6 | 2.2% | 1.5% | 0.0% |
| T5 | 9,844 | 5,566 | 1,693 | 4,968 | 7,731 | 0.4 | 1.1 | 9.9% | 11.4% | 0.1% |
| T6 | 20,253 | 4,291 | 1,642 | 3,695 | 6,086 | 7.7 | 20.1 | 3.0% | 1.8% | 10.0% |

Caps at the time of this audit were `max_features = 5000`, `max_counts = 25000`,
`max_percent_mt = 20`, `min_counts = 500` — the base blood values, since `bat_wing` did not touch
`QC`. **[fixed 2026-07-28]** the overlay now sets `max_features = 8000` and `max_counts = 60000`;
the other two are unchanged. On the T1/T2/T6 run this cut QC losses to 1–2% for T1 and T2.

- `max_counts = 25000` discards **10.4% of T1 and 11.4% of T5**, and `max_features = 5000` another
  6.1%/9.9%. These are the high-complexity cells — keratinocytes and fibroblasts are large and
  transcriptionally rich, so the cap removes the cell types of interest, not doublets. Suggest
  `max_features ≈ 8000`, `max_counts ≈ 60000` for tissue, and let scDblFinder handle doublets.
- `min_counts = 500` is inert: 0.0% of cells fall below it in any sample (`filter_matrix` is already
  cell-called). Note `01_load_qc.R` also applies `min.features = 200` at `CreateSeuratObject`.
- `max_percent_mt = 20` costs **10.0% of T6** and almost nothing elsewhere.

### Low mito is normal for this platform — T6 is the only real outlier

These are whole cells, not nuclei. The low mitochondrial fractions are a property of this
platform/workflow, not of these samples: across **40 samples from 22 prior runs**, median %MT is
0.04–2.68 (q1 0.47, median 0.80, q3 1.26). The human PBMC quick-start samples sit at the very bottom
(H1 0.04%, H2 0.06%), so this is not a bat-annotation artifact — it spans both species on the same
DNBelab C4 workflow.

| sample set | median %MT |
|---|---|
| 40 prior samples (bat + human blood) | 0.04 – 2.68 (median 0.80) |
| T1, T3, T4, T5 | 0.52, 0.45, 0.60, 0.42 — squarely in range |
| T2 | 2.24 — high end, still within range |
| **T6** | **7.65 — ~3× the highest value ever recorded here** |

Mito reads are being counted correctly, so this is not a quantification failure: the per-cell
distributions have proper right tails (T1 p99 = 7.6%, max = 63.3%; T2 max = 86.1%). The bulk of cells
is simply healthy.

Two consequences:

1. **Don't raise `max_percent_mt` for wing tissue.** Relative to a ~0.5–1% baseline, 20% is already
   extremely permissive and functions as a near-inert filter (≤1% of cells in a normal sample). The
   thresholds actually worth changing are `max_features`/`max_counts` above. Note the pre-2026-07-27
   docs claimed bat_wing raised `max_percent_mt`; it never did, and it should not.
2. **Flag T6 as a sample-quality concern, not a threshold problem.** At 7.65% median with 10% of cells
   over 20%, it is a genuine outlier against 40 historical samples and is also the most
   blood-contaminated tissue sample (HBB in 100% of cells). Treat its stromal proportions with
   caution rather than loosening QC to accommodate it.
- **Pervasive ambient RNA.** HBB is detected in 98–100% of cells in *every* sample and CD68 in
  85–91%, including T3/T4 which have no myeloid structure. Treat low-level marker positivity as
  soup, not biology; consider DecontX/SoupX before annotation.

## 3. Annotation gaps for tissue

`bat_wing` correctly sets `SINGLER_REF <- "HumanPrimaryCellAtlas"` (the blood-only MonacoImmune would
be useless here). HPCA `label.main` has 36 labels; `SINGLER_NORM` (56 entries) now maps 25 and passes
11 through unchanged. The gaps below were found by this audit; the tissue-label ones were fixed on
2026-07-28 and are marked **[fixed]**.

- **[fixed] `Keratinocytes` (plural) was unmapped.** Config keys are singular `Keratinocyte`, so
  `SUBTYPE_MARKERS[["Keratinocyte"]]` never fired and the basal/suprabasal/wound split silently
  never happened — on T1, the most keratinocyte-rich sample. `Keratinocytes → Keratinocyte` is now
  in `SINGLER_NORM`, alongside `Chondrocytes`, `MSC`, `Tissue_stem_cells → MSC`, `Osteoblasts → MSC`
  and `Macrophage(s)`. `Fibroblasts → Fibroblast` was already mapped, so fibroblast refinement
  always worked.
- **HPCA cannot emit `Pericyte`, `Melanocyte`, `Myofibroblast`, or `Wound_keratinocyte` at all.**
  `MARKERS` defines them, but they can only be assigned via `CLUSTER_CELLTYPE_MAP` or
  `SUBTYPE_MARKERS` refinement. Real populations are present (RGS5 in 20.4% of T2, MLANA in 4.9% of
  T1), so expect them to hide inside `Smooth_muscle_cells`/`Fibroblasts` until mapped by hand.
- **[fixed]** `Chondrocytes`, `MSC`, `Tissue_stem_cells`, `Osteoblasts` were unmapped and are
  plausible for wing membrane; all four now normalise and have colours.
- **HPCA swaps fibroblast and smooth muscle in wing tissue.** Confirmed on the T1/T2/T6 run: the
  cluster HPCA calls `Smooth_muscle_cells` is DCN 176 / COL1A1 26 / LUM 37 with ACTA2 1.0 and
  MYH11 0.6 (a fibroblast), while the one it calls `MSC` is ACTA2 60 / MYH11 22 / TAGLN 42 (the real
  smooth muscle). This is not fixable in `SINGLER_NORM` — the reference itself is wrong for this
  tissue — so it needs a `CLUSTER_CELLTYPE_MAP`. See the guarded example in
  [reference-config.md](reference-config.md#manual-cluster-annotation-cluster_celltype_map).
- **`T_cells → "CD4 T"`** in `SINGLER_NORM` coerces every generic HPCA T cell into CD4, then triggers
  the blood-specific naive/effector/memory refinement. Wrong for tissue T cells.
- The ScType consensus pass is explicitly blood-scoped — `05_annotate.R:253-255` already says "Best
  for PBMC / whole blood … novel tissue still needs manual review." Expect
  `consensus_annotation.csv` to be low-value here, and the scType override that propagates
  RBC/Platelet/Eosinophil/Mast is a blood construct.

**[fixed] `CELLTYPE_COLORS` had no colour for any wing label.** At 35 entries it covered only
`Fibroblast` and `Endothelial`; `Myofibroblast`, `Keratinocyte`, `Wound_keratinocyte`, `Pericyte`,
`Macrophage`, `Melanocyte` and every wing sub-type (`Fibroblast (resting)`, `Keratinocyte (basal)`,
`Macrophage (M1/inflam)`, …) were missing. That made plot colours unstable run-to-run and was a hard
blocker on the annotation loop, because `validate_config.R` errors when a `CLUSTER_CELLTYPE_MAP`
type has no colour — so the documented loop (read the dotplot → fill the map → re-run `05 06 07`)
failed on the second iteration. The palette is now 52 entries covering all wing labels plus
`Chondrocyte`, `MSC` and `Melanocyte/Schwann`.

## 4. Marker genes absent from this bat annotation

Config markers checked against the T1 feature list. Everything resolves except the genes below.
The four that affect the wing overlay (`CTGF`, `KRT5`, `LOR`, `TRGC2`) were applied to
`config_species_bat.R` on 2026-07-28; the bat whole-blood substitutions were already in place.

| config gene | status | fix |
|---|---|---|
| `CTGF` (SUBTYPE_MARKERS Fibroblast, WOUND_MODULES Myofibroblast) | absent | renamed — use `CCN2` (present) |
| `KRT5` (MARKERS/SUBTYPE_MARKERS Keratinocyte) | genuinely absent from the annotation | `KRT14`, `KRT15`, `KRT17` present and sufficient for basal |
| `LOR` (Keratinocyte suprabasal) | absent | `IVL`, `FLG` present |
| `TRGC2` (γδ T) | absent | `TRGC1`, `TRDC` present |
| `FCGR3A` (Monocyte sub-type) | absent | `FCGR3B` present |
| `CST3` (CD14_mono) | absent | drop; `CD14`/`LYZ`/`S100A8` cover it |
| `IGHD`, `IGHM`, `IGHG1` (B cell sub-types) | absent | `IGHA1`, `JCHAIN`, `MZB1` present |

Wing `MARKERS` are otherwise complete: Fibroblast 6/6, Myofibroblast 4/4, Endothelial 5/5,
Pericyte 4/4, Macrophage 4/4, Melanocyte 4/4, Wound_keratinocyte 4/4. `WOUND_MODULES` are complete
except `CTGF`.

## 5. Steps 11–14 and 06b

- **`CellChat`, `monocle3`, and `SeuratWrappers` are not installed** in `scrna_seurat`, and steps 13
  and 14 `library()` them unconditionally — both will hard-fail. `setup_env.sh` intends to install
  them from GitHub; that evidently did not complete.
- **11, 13 need `SCRNA_CONDITION` with ≥2 levels** (`condition` is derived from
  `SAMPLE_CONDITIONS[merged$sample]`, never stored in the RDS); 11 and 13 skip cleanly with one
  level, **14 does not guard** and assumes exactly two condition colours.
- **`14_trajectory.R` hardcodes cell-type names** (`Fibroblast`, `Myofibroblast`,
  `Fibroblast (resting)`, `Fibroblast (wound)`, `Macrophage`, `Macrophage (M1/inflam)`, …). These
  only exist if §3's label mapping is fixed first, otherwise every lineage is skipped.
- **`12_pathways.R` hardcodes the legends** `"Up in recovering"` / `"Up in healthy"` regardless of
  the actual condition levels — a mislabelling risk on wing figures. Its `org.Hs.eg.db` + KEGG
  `"hsa"` usage is defensible (this annotation uses human symbols) and documented at `12:5`.
- **`06b_differential.R` exits without doing anything unless there are exactly 2 samples**
  (`06b:40-44`), so it contributes nothing to a six-sample run. Its `SAMPLE_COLORS` also only
  defines H1/H2, giving NA fills in the module-score plot for other sample names.

## 6. Practical notes

- ~90K cells and ~210M non-zero entries across the six samples — several times the blood cohorts.
  Sparse counts alone are ~2.5 GB; the RAM governor reported 110 GB available and 5 merge workers, so
  it fits, but per-cell SingleR against HPCA on 90K cells is the runtime bottleneck.
- `DOUBLET$doublet_rate` is `list(H1=…, H2=…)`; `DOUBLET$doublet_rate[["T1"]]` returns `NULL`, which
  `02_doublets.R:42` treats as auto-estimate. Correct behaviour, no change needed.
- Six samples would exceed the 4-sample naming threshold and collapse to
  `Results/results_T1_6samples_filtered/`. Moot in practice: the cohort split (§8) means two
  three-sample runs, which keep the readable dash-joined names.
- Minor: `config.R` cannot be sourced from a script outside `pipeline/` — `.pipeline_dir` falls back
  to `"."` and `source("./pdf_helpers.R")` fails. Ad-hoc analysis scripts must `setwd("pipeline")`
  first. This is the path-resolver issue noted in `docs/eng-audit-2026-06-11.md`.

## 7. Suggested order of work

1. **[done]** Confirm what T1–T6 actually are, and decide the cohort split. Nothing else is worth
   doing until the integration grouping is right. Result: T1/T2/T6 are wing tissue and were run
   under `bat_wing`; T3/T4 are sorted CD45+ leukocytes and T5 is whole blood, run separately under
   `bat`. See §8.
2. **[done]** Config-only, small: raise `max_features`/`max_counts` for tissue in the `bat_wing`
   block and leave `max_percent_mt` at 20; add wing entries to `CELLTYPE_COLORS` (unblocks the
   annotation loop); add `Keratinocytes`/`Chondrocytes`/`MSC`/`Tissue_stem_cells` to `SINGLER_NORM`;
   apply the §4 gene substitutions.
3. **[done]** Run `01`–`07` on the tissue cohort and do the manual `CLUSTER_CELLTYPE_MAP` pass —
   pericytes, melanocytes and myofibroblasts will need it since HPCA cannot name them. The map that
   shipped fixes the fibroblast/smooth-muscle swap and names the melanocyte/Schwann and epithelial
   clusters; pericytes and myofibroblasts did not separate at res 0.5 and remain unresolved.
4. **Not done.** Only if the wound-healing analysis is wanted: install
   `CellChat`/`monocle3`/`SeuratWrappers`, define `SCRNA_CONDITION`, and fix the hardcoded labels in
   `12` and `14`. Steps 11–14 were skipped on the T1/T2/T6 run for exactly these reasons.

## 8. What was actually run (2026-07-28)

Two cohorts, split on the §1 finding that the six samples are three different sample types.

| Run | Species | Samples | Cells retained | Results dir |
|---|---|---|---|---|
| Wing tissue | `bat_wing` | T1, T2, T6 | 25,924 | `Results/results_T1-T2-T6_filtered/` |
| Blood / leukocyte | `bat` | T3, T4, T5 | see run log | `Results/results_T3-T4-T5_filtered/` |

Wing composition after the manual map (% of sample, T1 / T2 / T6): Fibroblast 15.6 / 50.3 / 32.5 ·
Keratinocyte (basal) 55.2 / 3.9 / 7.9 · RBC 7.5 / 3.4 / 18.6 · Endothelial 2.8 / 12.5 / 15.6 ·
CD4 T (memory) 3.9 / 11.7 / 9.3 · Epithelial 0.1 / 0.2 / 9.0 · Mast cell 11.4 / 5.5 / 0.4 ·
Smooth Muscle 0.5 / 5.1 / 3.3 · CD14+ Mono 1.0 / 4.3 / 2.9 · Melanocyte/Schwann 1.5 / 3.1 / 0.3.

Three things worth carrying forward:

- **The three wing samples are not equivalent layers.** T1 is epidermis-dominant (KRT14 in 99.5% of
  barcodes, collagen 31%), T2 is full-thickness, T6 is deep dermis with the epidermis essentially
  absent (KRT14 1.1%). Cross-sample proportion differences are dissection depth first, biology
  second.
- **Labels are run-relative.** No frozen reference exists for wing tissue, so `05r`/`08c` self-skip
  and `options("scrna.label_source")` reads `denovo`. Adding a sample changes the clusters *and*
  silently disables the guarded `CLUSTER_CELLTYPE_MAP`.
- **Two clusters remain uncertain.** The `Epithelial` cluster (KRT18/PAX9/PSMB11, 99% T6) is
  glandular but not pinned to a gland type, and `Melanocyte/Schwann` is two populations merged at
  res 0.5.
