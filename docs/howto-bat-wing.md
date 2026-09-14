# How to Run on Bat Wing Tissue

This guide covers running the pipeline on *Eonycteris spelaea* wing (wound-healing) tissue. The
`bat_wing` species keyword activates a separate set of overrides from `bat` (whole blood) — different
reference, QC thresholds, marker panels, and four extra steps (`11`–`14`) for condition-contrast
wound-healing analysis. The worked example throughout is the real `T1`/`T2`/`T6` run in
`Results/results_T1-T2-T6_filtered/`.

See [How to Run on Bat Whole Blood](howto-bat-whole-blood.md) for the blood-side `bat` keyword; the
two overlays share nothing except the source file (`config_species_bat.R`) and should not be applied
to the same sample.

---

## Prerequisites

- Conda environment active: `conda activate scrna_seurat`
- 10x/DNBelab matrices with HGNC-style gene names (same orthology mapping as `bat`)
- For steps `13` (CellChat) and `14` (monocle3 trajectory): `CellChat`, `monocle3` and
  `SeuratWrappers`. As of 2026-09-14 `monocle3` is installed in `scrna_seurat`, `CellChat` is being
  installed, and **`SeuratWrappers` is not installable** in this env (its Bioconductor 3.18
  dependency chain conflicts). Step `14` calls `SeuratWrappers::as.cell_data_set()`, so it will
  hard-fail until it builds the monocle3 object directly, as `pipeline/projects/cm/trajectory.R`
  does. Check with `Rscript -e 'requireNamespace("CellChat")'` before running `13`.

---

## When to use `bat_wing` vs `bat`

`bat_wing` is for dissected wing tissue (skin/dermis/epidermis); `bat` is for whole blood or
sorted blood fractions. **Do not assume a batch of samples is one tissue type** — verify first. The
readiness audit for `Samples/T1`–`T6` found the six samples split into three groups by raw marker
expression (`COL1A1`, `DCN`, `KRT14` for tissue; `PTPRC`/`CD3E`/`S100A8` for blood):

| Group | Samples | Evidence | Species keyword |
|---|---|---|---|
| Wing tissue | T1, T2, T6 | COL1A1/DCN/KRT14 positive, PTPRC 8–28% | `bat_wing` |
| CD45-sorted immune infiltrate | T3, T4 | PTPRC 91–99%, ~0% collagen/keratin | `bat` |
| Whole blood | T5 | S100A8 99.9%, HBB 6.4% of UMIs | `bat` |

This matters because `HARMONY$group_by_vars <- "sample"` corrects *per sample* — putting tissue,
sorted cells, and whole blood into one Harmony run regresses out the tissue difference as if it were
batch effect. Run tissue samples as their own cohort.

---

## What the `bat_wing` overlay changes

`config_species_bat.R` is sourced by `config.R` after the human base definitions; the `bat_wing`
branch mutates:

| Parameter | Human/blood default | `bat_wing` override |
|-----------|---------------------|----------------------|
| `SINGLER_REF` | `"HumanPrimaryCellAtlas"` (human) / `"MonacoImmune"` (bat) | `"HumanPrimaryCellAtlas"` — MonacoImmune is blood-only |
| `QC$max_features` / `QC$max_counts` | 5000 / 25000 | **8000 / 60000** — blood caps clip 5–11% of real keratinocytes/fibroblasts, which are larger and more transcript-rich than blood cells |
| `QC$max_percent_mt` | 20 | unchanged — stays 20 for every species; do not raise it for tissue |
| `CLUSTER$resolutions` / `default_res` | human 0.5, bat 1.0 | `c(0.3, 0.5, 0.8)`, default 0.5 |
| `MARKERS` | PBMC/blood panels | Fibroblast, Myofibroblast, Keratinocyte, Wound_keratinocyte, Endothelial, Pericyte, Macrophage, Melanocyte, tissue-specific Neutrophil, γδ T (`TRDC`, `TRGC1`) |
| `SUBTYPE_MARKERS` | blood subsets | Fibroblast (resting/myofibroblast/wound), Macrophage (M1/M2/proliferating), Keratinocyte (basal/suprabasal/wound) |
| `WOUND_MODULES` | not defined | 6 gene modules (Inflammatory, ECM_remodeling, Angiogenesis, Proliferation, Re_epithelialize, Myofibroblast) — consumed by step `11` |
| `CONTAMINATION_TYPES` | includes Neutrophil (blood) | reduced to `RBC`, `HSPC`, `Platelet` — no blood-construct types in tissue |
| `CELLTYPE_COLORS` | — | not touched by the overlay itself, but wing labels must exist in the base palette or `validate_config.R` hard-errors on any `CLUSTER_CELLTYPE_MAP` entry lacking a colour |

**Cache-invalidation consequence:** step 01's cache key (`cache_hash()`) hashes `QC` along with
`DOUBLET`/`NORM`/`DIM`/`CLUSTER` + species. Because `bat_wing`'s `QC$max_features`/`max_counts` differ
from both `human` and `bat`, running the same sample under `bat_wing` after it was cached under `bat`
(or vice versa) invalidates `sample_cache/<sample>/` and forces steps 01–03 to reprocess.

---

## Running the pipeline: the T1/T2/T6 example

```bash
conda activate scrna_seurat
bash pipeline/run_pipeline.sh bat_wing Samples/T1 Samples/T2 Samples/T6 01 02 03 04 05 06 07
```

No transcript of the original invocation survives; this command is reconstructed from the run's
`logs/` directory, which holds `01_01_load_qc`, `02_02_doublets`, `03_03_individual`,
`04_04_integrate`, `05_05_annotate`, `05r_reference_transfer`, `06_06_visualize`,
`07_07_finalize_reports`, `08b_html_report` and `08c_benchmark_concordance` — no `06b` and no
`11`–`14` logs, so explicit step ids were evidently passed. Either way:

- `06b_differential.R` does nothing unless the run has exactly 2 samples (it has 3 here).
- `11`–`14` need a two-level `SCRNA_CONDITION`, and none was defined for this cohort (see
  [Steps 11–14](#steps-11-14-wound-healing-analysis)).

Left un-stepped (`bash pipeline/run_pipeline.sh bat_wing <samples>`), the default step list for
`bat_wing` is `01 02 03 04 05 06 06b 07 11 12 13 14` — wider than the blood default, so pass explicit
step ids if you want to skip the wound-healing steps or 06b.

The run retained **25,924 cells** across T1/T2/T6.

---

## The manual annotation loop

Same shape as the blood doc, with a tissue-specific twist: leave `CLUSTER_CELLTYPE_MAP <- NULL` on
the first run and let SingleR (HumanPrimaryCellAtlas) auto-annotate. After step 05:

```bash
grep -A 30 "CLUSTER_CELLTYPE_MAP" Results/results_T1-T2-T6_filtered/logs/05_annotate.log
```

Review `annotation/canonical_markers_dotplot.pdf`. HPCA has 36 `label.main` labels; it does not know
tissue-specific labels at all — `Pericyte`, `Melanocyte`, `Myofibroblast`, and `Wound_keratinocyte`
can only be assigned via `CLUSTER_CELLTYPE_MAP` or `SUBTYPE_MARKERS` refinement, never emitted
directly by HPCA. Confirmed on the T1/T2/T6 run: HPCA also **swaps fibroblast and smooth muscle** —
the cluster it calls `Smooth_muscle_cells` is the real fibroblast population (DCN 176 / COL1A1 26 /
LUM 37, ACTA2 1.0 / MYH11 0.6), and the one it calls `MSC` is the real smooth muscle (ACTA2 60 /
MYH11 22 / TAGLN 42). This is not a `SINGLER_NORM` fix — the reference itself is wrong for wing
tissue — so it requires a `CLUSTER_CELLTYPE_MAP`.

The map that shipped for this run is guarded so it cannot leak into a run with different samples or
cluster numbering:

```r
if (length(SAMPLE_NAMES) == 3 && setequal(SAMPLE_NAMES, c("T1", "T2", "T6"))) {
  CLUSTER_CELLTYPE_MAP <- c(
    "0"  = "Fibroblast",
    "1"  = "RBC",
    "2"  = "Keratinocyte (basal)",
    ...
    "19" = "RBC"
  )
}
```

Copy this `setequal(SAMPLE_NAMES, ...)` pattern for any new sample set rather than assigning a map at
top level — cluster numbers are not stable across runs. All 20 clusters must be listed: unmapped
clusters fall back to **per-cell** SingleR labels (`05_annotate.R:379`), so a partial map shatters
every cluster left out of it, not just the ones you intended to leave alone.

Final composition after the map (from `docs/bat_wing_readiness.md`; `cluster_annotation_table.csv`
holds only per-cluster counts), % of sample T1 / T2 / T6: Fibroblast 15.6 / 50.3 / 32.5 · Keratinocyte (basal) 55.2 / 3.9 / 7.9 ·
RBC 7.5 / 3.4 / 18.6 · Endothelial 2.8 / 12.5 / 15.6 · CD4 T (memory) 3.9 / 11.7 / 9.3 ·
Epithelial 0.1 / 0.2 / 9.0 · Mast cell 11.4 / 5.5 / 0.4 · Smooth Muscle 0.5 / 5.1 / 3.3 ·
CD14+ Mono 1.0 / 4.3 / 2.9 · Melanocyte/Schwann 1.5 / 3.1 / 0.3.

Two clusters remain uncertain even after the map: `Epithelial` (KRT18/PAX9/PSMB11, 99% T6) is
glandular but not pinned to a gland type, and `Melanocyte/Schwann` is two populations merged at
resolution 0.5. Pericytes and myofibroblasts did not separate from smooth muscle/fibroblast at this
resolution and remain unresolved.

Re-run after editing the map:

```bash
bash pipeline/run_pipeline.sh bat_wing Samples/T1 Samples/T2 Samples/T6 05 06 07
```

---

## Steps 11–14 (wound-healing analysis)

These are `bat_wing`-only, live in `pipeline/projects/bat_wing/`, and all source the core
`config.R`. They compare two `SCRNA_CONDITION` groups (e.g. healthy vs. wound/recovering) — set with
`condition=<sample1>=<label1>,<sample2>=<label2>,...` as an arg to `run_pipeline.sh`. None of them ran
on the T1/T2/T6 cohort (no condition contrast was defined for that run).

| Step | Requires | Outputs | Precondition |
|---|---|---|---|
| `11_wing_degs.R` | `integrated_annotated.rds` (steps 04+05) | DEG CSVs per cell type, volcano plots, `WOUND_MODULES` score UMAPs/violins, top-DEG heatmap → `DIRS$differential` | **Skips cleanly** if `length(CONDITION_LEVELS) < 2` |
| `12_pathways.R` | `all_DEG_combined.csv` from step 11 | GO/KEGG enrichment CSVs, bar plots, `pathways_report.pdf` → `DIRS$pathways` | Needs step 11's output; uses `org.Hs.eg.db` + KEGG `"hsa"` (defensible — bat annotation uses human gene symbols) |
| `13_cellchat.R` | `integrated_annotated.rds` (steps 04+05) | One CellChat object per condition, chord diagrams, bubble plots, differential signalling heatmap, `cellchat_report.pdf` | **Skips cleanly** if `length(CONDITION_LEVELS) < 2` (same guard as 11) |
| `14_trajectory.R` | `integrated_annotated.rds` (steps 04+05) | Pseudotime UMAPs, trajectory graph, top pseudotime DEGs, `trajectory_report.pdf` | **Does not guard** condition count — assumes exactly two condition colours; also hardcodes cell-type names (`Fibroblast`, `Myofibroblast`, `Fibroblast (resting)`, `Fibroblast (wound)`, `Macrophage`, `Macrophage (M1/inflam)`, …), so it only produces output once the `CLUSTER_CELLTYPE_MAP` pass above has named those clusters |

`12_pathways.R` also hardcodes its GSEA legend as `"Up in recovering"` / `"Up in healthy"` regardless
of the actual `SCRNA_CONDITION` labels used — check the legend against your real condition names
before trusting the figure.

---

## Gotchas

- **T1–T6 are not one tissue.** Check raw marker percentages (`COL1A1`, `DCN`, `KRT14` for tissue vs.
  `PTPRC`, `CD3E`, `S100A8` for blood) before deciding the cohort — see
  [When to use](#when-to-use-bat_wing-vs-bat) above.
- **Don't raise `max_percent_mt` for wing tissue.** This DNBelab C4 workflow runs a ~0.5–1% mito
  baseline (median 0.04–2.68% across 40 prior samples, both species); 20% is already a near-inert
  filter. `T6`'s 7.65% median mito is a genuine sample-quality outlier (~3× the highest value ever
  recorded on this platform) and should be flagged, not accommodated by loosening QC.
- **Ambient RNA is pervasive.** `HBB` shows in 98–100% of cells and `CD68` in 85–91%, in every
  sample including ones with no myeloid or blood structure. Treat low-level marker positivity as
  ambient signal, not biology.
- **Labels are run-relative.** No frozen SingleR reference exists for wing tissue, so `05r` and `08c`
  self-skip and `options("scrna.label_source")` reads `"denovo"`. Adding or removing a sample changes
  the clustering *and* silently disables the `setequal`-guarded `CLUSTER_CELLTYPE_MAP` — expect to
  redo the manual annotation loop for any cohort other than exactly T1/T2/T6.
- **The ScType consensus pass is blood-scoped** (see the "closed cell-type vocabulary" comment in
  `05_annotate.R`) — expect
  `consensus_annotation.csv` to be low-value on tissue, and its RBC/Platelet/Eosinophil/Mast override
  is a blood construct that doesn't apply here.
- **`06b_differential.R`** silently contributes nothing unless the run has exactly 2 samples.

---

## Related

- [How to Run on Bat Whole Blood](howto-bat-whole-blood.md) — the `bat` (blood) counterpart
- [Configuration Reference](reference-config.md) — species overrides table and the manual
  `CLUSTER_CELLTYPE_MAP` pattern
- `docs/bat_wing_readiness.md` — the full readiness audit this guide is drawn from, including the
  T1–T6 sample-identity investigation and per-sample QC distributions
- [Annotation Strategy](explanation-annotation.md) — understanding contamination-type overrides and
  the SingleR normalisation layer
