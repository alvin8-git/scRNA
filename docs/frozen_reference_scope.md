# Frozen-reference benchmark — build scope

**Goal:** run-independent cell-type labels so shared anchor samples (Aksh1, ES332) get
identical proportions regardless of co-samples → a real new-vs-previous benchmark.
Proven viable: a SingleR transfer resolved ES17 cluster 0 (86.8% Neutrophil) where the
within-run consensus deadlocked. See [memory: project_es17_neutrophil_unresolvable].

## Locked design decisions

1. **Engine: SingleR `trainSingleR()` once → `classifySingleR()` per run.**
   The frozen artifact *is* the trained model. Earlier 20-min runs were untrained SingleR
   recomputing markers every call; a trained model precomputes de.genes → classify is
   minutes, `fine.tune=TRUE` affordable. Already a pipeline dependency; gives score + delta
   + pruned labels for confidence handling.
2. **Reference source = in-house curated run, NOT an external atlas.**
   Data is MGI DNBelab C4, not 10x/microarray — a platform-matched reference transfers far
   better than HPCA/Monaco. Build from a high-quality run; **hold OUT Aksh1 + ES332** so the
   anchor benchmark isn't circular.
3. **Curation:** keep only cells in AUTO-consensus clusters (SingleR=scType) with mean
   SingleR delta above threshold; drop REVIEW clusters + low-confidence cells. Clean labels in,
   clean transfer out.
4. **Additive, never replace.** Transfer writes `cell_type_ref` / `ref_score` / `ref_pruned`
   alongside the existing de-novo `cell_type`. De-novo annotation stays intact.
5. **Confidence gate:** pruned.labels + score/delta threshold; below → `Unassigned (low-conf)`;
   report % unassigned per sample (the ES17 cDC2 over-call in the quick test shows why this matters).

## Components to build

| File | Role | When run |
|------|------|----------|
| `pipeline/build_reference.R` | filter run → curated cells → `trainSingleR` → save model+metadata | rarely, manual |
| `pipeline/05r_reference_transfer.R` | if `REFERENCE_MODEL` set: `classifySingleR` → add ref columns | every run (gated) |
| `pipeline/08c_benchmark_concordance.R` | per-anchor cell-type % vs baseline → delta table + drift flags + report | every run (gated) |

**config.R additions:** `REFERENCE_MODEL` (path, default NULL), `ANCHOR_SAMPLES`
(`c("Samples/Aksh1","Samples/ES332")`), `REF_SCORE_MIN`, `DRIFT_FLAG_PP` (default 5).

**run_pipeline.sh:** call `05r` after `05` (gated on `REFERENCE_MODEL`); `08c` after `08b`
(gated on `REFERENCE_MODEL` + `ANCHOR_SAMPLES`). Both additive; never fail the pipeline
(mirror the existing 08b warn-don't-die pattern).

**Model artifact bundles:** trained model, label vocabulary, per-label n, gene set, source run,
held-out anchors, git commit, date, **baseline anchor proportions**. Versioned filename.

## Acceptance test (go/no-go)

Transfer the model onto BOTH the ES17 and ES03 runs; Aksh1 + ES332 cell-type % must match
within `DRIFT_FLAG_PP` (≤5 pp) across the two runs. That proves run-independence — the whole point.

## Decisions (locked 2026-06-25)

1. **Reference source = ES03 run**, with Aksh1 + ES332 held out of the build.
2. **Drift tolerance = 5 pp** to flag anchor drift.
3. **Trained model** (`trainSingleR`), not the quick untrained path.

## Species / sort status (verified from metadata)

- ES03, ES332 (and the whole reference + ES17 runs) are **bat whole blood — _Eonycteris spelaea_**
  (`E250183788_metrics.xlsx` Species field), NOT human PBMC. Annotated via human references +
  bat marker overrides (`config_species_bat.R`); gene symbols are human orthologs.
- **ES03/ES14/ES251/ES258/ES332/ES35/ES407/ES459 are all PRESORT** (before sorting = whole blood,
  population-complete) — confirmed by user, 2026-06-25. Matches the H-series `presort` convention
  and the ~9,880 neutrophils retained in the ES03 run (a postsort PBMC prep depletes granulocytes).
- This is the ideal reference property: a presort/whole-blood reference carries ALL populations
  (incl. neutrophils), so it can label any query — even a more-depleted one — without missing types.
- Reinforces ES03 as reference: **species- and tissue-matched** to the ES17 query AND population-complete.
- ES17 main batch (ES17/18/122/129/171/354): sort status unknown. User leans postsort; however
  the RBC/platelet/neutrophil populations PRESENT in that run usually indicate presort/whole blood
  (standard PBMC sorting depletes them) — pending lab confirmation. Does NOT affect Phase 1: a
  presort reference is a population superset and labels any query (incl. postsort) without gaps.
- ES = _Eonycteris spelaea_ (bat); all ES + Aksh samples same species as anchors. Species-matched.
- Phase 3-4 design feed (NOT Phase 1): `08c` must report RBC/platelet/neutrophil fractions per
  batch (direct sort-status readout) and flag that the ES17 main batch transfers across a batch
  (and possibly sort/composition) gap from the presort reference.
- CAVEAT: the frozen reference inherits the ES03 run's bat-override annotation — labels are
  consistent, not ground-truth bat cell types. Run `build_reference.R` with `SCRNA_SPECIES=bat`.

## Batch structure caveat (important)

ES17 run = one batch (ES17/18/122/129/171/354) + the anchors (Aksh1/ES332) from a *different*
batch. This does NOT change the decisions, but splits the output into two confidence tiers:

- **Anchors are batch-matched to the ES03 reference** (they live in ES03's batch). And the anchor
  cells are the *same input matrix* re-analyzed in both runs — so anchor drift between runs is
  **purely pipeline artifact** (no batch or biology difference). That makes the benchmark a clean,
  best-case test: the frozen reference should drive Aksh1/ES332 drift to ~0. Strong.
- **The ES17 main batch transfers ACROSS a batch gap** from the ES03-batch reference. SingleR does
  not batch-correct, so main-batch labels (including the cluster-0 86.8% Neutrophil call) carry
  cross-batch uncertainty — could be real neutrophils or batch-shifted myeloid cells. **Treat
  main-batch labels as provisional; confirm the neutrophil identity within-batch.**

Consequences for the build:
- Keep ES03 as reference (anchors batch-matched = reliable benchmark). Do NOT exclude the anchor
  batch — we want anchors well-matched; we are measuring run-composition reproducibility, not
  cross-batch robustness.
- `08c` report must **split prediction score / delta by batch** (anchor vs main) — a direct
  readout of whether main-batch transfers are systematically lower-confidence (batch effect).
- Production note: if accurate labeling of *new batches* (not just the anchor benchmark) becomes a
  goal, a single-batch ES03 reference will transfer less reliably. Options then: a multi-batch
  reference, or Seurat `MapQuery` (integrates query into reference space, batch-aware). Keep SingleR
  for the anchor benchmark regardless.

## Effort (CC-paced)

Phase 1 `build_reference.R` + build model · Phase 2 `05r` + config · Phase 3 `08c` + report ·
Phase 4 validate on ES17/ES03 anchors. Each phase minutes of CC + one model-build run
(~loads the source rds once). No change to steps 01–05 de-novo path.
