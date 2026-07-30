# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A single-cell RNA-seq pipeline in R (Seurat v5 + Harmony + SingleR), driven entirely by shell/Rscript
— no package, no build step. Input is count matrices from MGI DNBelab C Series (C4) or 10x CellRanger;
output is per-step PDFs plus one self-contained interactive HTML report. Supported species/tissues:
human PBMC, human whole blood, bat (*Eonycteris spelaea*) whole blood, and bat wing tissue.

## Commands

Everything runs inside the `scrna_seurat` conda env (R 4.4.3, Seurat 5.4.0). `run_pipeline.sh`
auto-activates it; standalone `Rscript` calls need `conda activate scrna_seurat` first.

```bash
bash pipeline/setup_env.sh                       # create the conda env (mamba, 10-20 min, one time)

bash pipeline/run_pipeline.sh Samples/H1 Samples/H2        # full run, bundled human example data
bash pipeline/run_pipeline.sh /path/A                      # one sample -> no integration
bash pipeline/run_pipeline.sh bat /path/ES03 /path/ES12    # species keyword: bat | human | bat_wing
bash pipeline/run_pipeline.sh bat /path/A /path/B 05 06 07 # only these steps (resume from checkpoints)
bash pipeline/run_pipeline.sh 07                           # steps only; samples come from config defaults
bash pipeline/run_pipeline.sh --help
```

Arg parsing is positional-agnostic: anything absolute / a directory / `./…` is a sample path, `bat`
`human` `bat_wing` set species, `condition=…` exports `SCRNA_CONDITION`, `--no-report` skips the HTML,
everything else is a step id. Valid steps: `01 02 03 04 05 06 06b 07 08 11 12 13 14`.

Checks (there is no CI, no Makefile, no test framework):

```bash
Rscript pipeline/validate_config.R                 # pre-flight; run_pipeline.sh runs this first
Rscript pipeline/validate_config.R --strict-paths   # escalate missing sample dirs from warning to error
Rscript pipeline/tests/test_regressions.R           # the whole suite (one grep-based regression guard)
```

`test_regressions.R` has no per-test selection — it is a flat list of "this old bug must not come
back" greps over the pipeline sources. Add new guards to it rather than creating new files.

Steps not in the default sets, run directly against a finished run dir:

```bash
bash pipeline/build_report.sh Results/results_A-B_filtered [--samples=A,B] [--max-cells=N]
SCRNA_RESULTS_DIR=Results/results_A-B_filtered Rscript pipeline/09_bootstrap_proportions.R
SCRNA_RESULTS_DIR=Results/results_A-B_filtered Rscript pipeline/10_rarefaction.R
SCRNA_SPECIES=bat Rscript pipeline/build_reference.R <ref_run_dir> --holdout=Aksh1,ES332
```

## Architecture

**`pipeline/config.R` is the whole configuration system.** Every step begins with
`source(file.path(.pipeline_dir, "config.R"))`, which defines ~30 loose globals (`QC`, `DOUBLET`,
`NORM`, `DIM`, `CLUSTER`, `HARMONY`, `MARKERS`, `SINGLER_REF`, `CLUSTER_CELLTYPE_MAP`,
`CELLTYPE_COLORS`, `PARALLEL`, `PLOT`, `DIRS`, …) — not a single config object. Steps are therefore
individually runnable; each resolves its own directory via a `sys.frame(1)$ofile` / `commandArgs()`
fallback. Changing pipeline behaviour almost always means editing `config.R`, not a step script.

**Species overlay is in-place mutation.** `config.R` defines human values, then sources
`config_species_bat.R`, which rewrites `MARKERS`, `ALL_MARKERS`, `CONTAMINATION_TYPES`, `SINGLER_REF`,
`CLUSTER$resolutions`, `SUBTYPE_MARKERS`, `WOUND_MODULES` under `if (.species == "bat")` /
`bat_wing`. `bat_wing` additionally raises `QC$max_features` (8000) and `QC$max_counts` (60000) for
tissue-sized cells — the only overlay that touches `QC`, and it invalidates the 01–03 sample cache.
`CELLTYPE_COLORS` is never overlaid; wing labels live in the base palette. Human is a no-op. Human *whole blood* has no keyword — set `SINGLER_REF <- "MonacoImmune"`
by hand.

**Run directory.** `RESULTS_DIR` = `Results/results_<sample1>-<sample2>-…_<filtered|raw>`, derived
identically by `run_pipeline.sh` (bash) and `config.R` (R) — keep those two in sync. More than 4
samples collapses to `results_<firstSample>_<N>samples_…` to stay under the Windows 260-char path
limit. `DIRS` maps subdirs (`qc/ doublets/ individual/ integrated/ annotation/ differential/ logs/`
…); note `DIRS$reports` is `RESULTS_DIR` itself, so step 07's PDFs land at the run root while 08b's
HTML goes in `reports/`.

**Checkpoint / resume contract.** Each step hands off a `.rds`:

| Step | reads | writes |
|---|---|---|
| 01 | `Samples/<n>/{filter_matrix,filtered_feature_bc_matrix,raw_*}` | `individual/<n>/<n>_filtered.rds` |
| 02 | `<n>_filtered.rds` | `<n>_singlets.rds` |
| 03 | `<n>_singlets.rds` | `individual/<n>/…` + `sample_cache/<n>/` |
| 04 | `sample_cache/<n>/…` | `integrated/integrated_seurat.rds` |
| 05 | `integrated_seurat.rds` | `integrated/integrated_annotated.rds` |
| 06, 06b, 07, 08, 08b, 09, 10 | `integrated/integrated_annotated.rds` | plots / PDFs / HTML |

`sample_cache/` at the repo root is a **run-independent** per-sample cache written by 01–03 and read
by 04, so a sample already processed in one run is free in the next. Its key (`cache_hash()`) is a
matrix-file fingerprint (size + mtime) plus `QC`/`DOUBLET`/`NORM`/`DIM`/`CLUSTER` + species —
so changing `HARMONY`, `MARKERS`, `SINGLER_REF`, or `CLUSTER_CELLTYPE_MAP` does **not** bust it.
Force reprocessing with `rm -rf sample_cache/<sample>/`.

**Cell-type labels have two independent sources.** Step 05 produces de-novo labels in
`merged$cell_type` (SingleR → `SINGLER_NORM` normalisation in `05_annotate_singler_norm.R` →
contamination overrides → manual `CLUSTER_CELLTYPE_MAP` → marker refinement → T-cell subclustering),
finalised as `final_cell_type`. Because SingleR is run-relative, labels shift with the sample mix; so
`05r_reference_transfer.R` classifies cells against a **frozen** SingleR model built once by
`build_reference.R`, adding run-independent `cell_type_ref`. `config.R`'s `apply_reference_labels()`
then promotes `cell_type_ref` over the de-novo call when that output exists, preserving the original
as `cell_type_denovo` and recording which is live in `options("scrna.label_source")`. `08c` benchmarks
anchor-sample composition against the frozen baseline and flags drift over `DRIFT_FLAG_PP` (5 pp).
Design notes: `docs/frozen_reference_scope.md`, `docs/howto-frozen-reference.md`.

**Additive stages never fail a run.** `run_pipeline.sh` runs `05r` (after 05), `08b` (HTML), and
`08c` outside the step loop, logging a warning instead of exiting on error, and they self-skip when
`REFERENCE_MODEL` is unset or `integrated_annotated.rds` is absent. Keep that property when editing
them. Bat-wing-only steps 11–14 live in `pipeline/projects/bat_wing/` and source the core `config.R`.

**Parallelism.** `PARALLEL$workers` (steps 01–03, ~8 GB each) vs `PARALLEL$merge_workers` (04–06b,
~16 GB each), capped by a `MemAvailable` governor. `run_pipeline.sh` pins BLAS/OMP/MKL to 1 thread to
stop future/BiocParallel workers from oversubscribing.

## Env vars

`SCRNA_SAMPLE1..N` (sample paths, exported by the runner), `SCRNA_SPECIES`, `SCRNA_CONDITION`
(`name=label,…`), `SCRNA_BASE_DIR`, `SCRNA_RESULTS_DIR` (point a step at an existing run dir instead
of re-deriving one from samples), `SCRNA_REFERENCE_MODEL` (enables 05r/08c), `SCRNA_ANCHORS` (default
`Aksh1,ES332`), `SCRNA_DRIFT_PP` (default `5`), `RSCRIPT` (override in `build_report.sh`).

## Gotchas

- With no `SCRNA_SAMPLE*` set, `config.R` silently falls back to hardcoded `H1`/`H2` + human.
- `CLUSTER_CELLTYPE_MAP` must stay `NULL` for a new dataset — cluster numbers are not stable across
  runs. The archived maps at the bottom of `config.R` are commented "do NOT activate". The one live
  map (bat wing T1/T2/T6) is wrapped in a `setequal(SAMPLE_NAMES, ...)` guard so it cannot leak into
  a run with different cluster numbering; copy that pattern rather than assigning the map at top level.
- Unmapped clusters in a partial `CLUSTER_CELLTYPE_MAP` fall back to **per-cell** SingleR labels, not
  the cluster majority (`05_annotate.R:379`) — so a map that lists only the clusters you want to fix
  will shatter every cluster you left out. Map all of them or none.
- Manual annotation loop: `logs/05_annotate.log` prints a paste-ready `CLUSTER_CELLTYPE_MAP`; fill it
  in `config.R` after reading `annotation/canonical_markers_dotplot.pdf`, then re-run `05 06 07`.
- `validate_config.R` hard-errors on only two things: a `CLUSTER_CELLTYPE_MAP` type missing from
  `CELLTYPE_COLORS`, and map keys that aren't bare integers.
- The `docs/` tree was reconciled against the code on 2026-07-28 (identifiers, output paths, metadata
  columns, cache contract, `bat_wing` QC/marker/palette overrides). It drifts easily — when you rename a config object, move an output, or add
  a metadata column, grep `docs/`, `DOCUMENTATION.md`, and `ReportGuide.md` for the old name.
- There is no `final_cell_type` column on the Seurat object; the cell metadata label is `cell_type`
  (plus `cell_type_ref` / `cell_type_denovo` after 05r). `final_cell_type` is only the per-cluster
  majority column inside `annotation/cluster_annotation_table.csv`.
- `DOUBLET$doublet_rate` is hardcoded per example sample (`list(H1=0.031, H2=0.077)`).
- `.rds`, `results*/`, `*.log`, `*.pptx` are gitignored; `Samples/H1` + `Samples/H2` matrices are
  force-tracked as the quick-start fixture.

## Doc map

`README.md` quick start · `DOCUMENTATION.md` per-step reference · `ReportGuide.md` how to interpret
each PDF · `VERSION.md` changelog · `TODO.md` backlog · `docs/howto-*.md` task recipes ·
`docs/explanation-*.md` design rationale · `docs/reference-*.md` config/output/step tables.

Study/audit notes that don't fit those four buckets: `docs/bat_wing_readiness.md` (what wing tissue
needs from the pipeline, plus the T1–T6 sample identities and the T1/T2/T6 run) ·
`docs/frozen_reference_scope.md` (design of the frozen-reference labels) ·
`docs/bat_neutrophil_literature.md` (validation of the bat neutrophil fraction) ·
`docs/eng-audit-2026-06-11.md` (cache/path-resolver findings) ·
`docs/design-interactive-report.md` + `docs/report_design_review.md` (HTML report).

---

# context-mode — MANDATORY routing rules

You have context-mode MCP tools available. These rules are NOT optional — they protect your context window from flooding. A single unrouted command can dump 56 KB into context and waste the entire session.

## BLOCKED commands — do NOT attempt these

### curl / wget — BLOCKED
Any Bash command containing `curl` or `wget` is intercepted and replaced with an error message. Do NOT retry.
Instead use:
- `ctx_fetch_and_index(url, source)` to fetch and index web pages
- `ctx_execute(language: "javascript", code: "const r = await fetch(...)")` to run HTTP calls in sandbox

### Inline HTTP — BLOCKED
Any Bash command containing `fetch('http`, `requests.get(`, `requests.post(`, `http.get(`, or `http.request(` is intercepted and replaced with an error message. Do NOT retry with Bash.
Instead use:
- `ctx_execute(language, code)` to run HTTP calls in sandbox — only stdout enters context

### WebFetch — BLOCKED
WebFetch calls are denied entirely. The URL is extracted and you are told to use `ctx_fetch_and_index` instead.
Instead use:
- `ctx_fetch_and_index(url, source)` then `ctx_search(queries)` to query the indexed content

## REDIRECTED tools — use sandbox equivalents

### Bash (>20 lines output)
Bash is ONLY for: `git`, `mkdir`, `rm`, `mv`, `cd`, `ls`, `npm install`, `pip install`, and other short-output commands.
For everything else, use:
- `ctx_batch_execute(commands, queries)` — run multiple commands + search in ONE call
- `ctx_execute(language: "shell", code: "...")` — run in sandbox, only stdout enters context

### Read (for analysis)
If you are reading a file to **Edit** it → Read is correct (Edit needs content in context).
If you are reading to **analyze, explore, or summarize** → use `ctx_execute_file(path, language, code)` instead. Only your printed summary enters context. The raw file content stays in the sandbox.

### Grep (large results)
Grep results can flood context. Use `ctx_execute(language: "shell", code: "grep ...")` to run searches in sandbox. Only your printed summary enters context.

## Tool selection hierarchy

1. **GATHER**: `ctx_batch_execute(commands, queries)` — Primary tool. Runs all commands, auto-indexes output, returns search results. ONE call replaces 30+ individual calls.
2. **FOLLOW-UP**: `ctx_search(queries: ["q1", "q2", ...])` — Query indexed content. Pass ALL questions as array in ONE call.
3. **PROCESSING**: `ctx_execute(language, code)` | `ctx_execute_file(path, language, code)` — Sandbox execution. Only stdout enters context.
4. **WEB**: `ctx_fetch_and_index(url, source)` then `ctx_search(queries)` — Fetch, chunk, index, query. Raw HTML never enters context.
5. **INDEX**: `ctx_index(content, source)` — Store content in FTS5 knowledge base for later search.

## Subagent routing

When spawning subagents (Agent/Task tool), the routing block is automatically injected into their prompt. Bash-type subagents are upgraded to general-purpose so they have access to MCP tools. You do NOT need to manually instruct subagents about context-mode.

## Output constraints

- Keep responses under 500 words.
- Write artifacts (code, configs, PRDs) to FILES — never return them as inline text. Return only: file path + 1-line description.
- When indexing content, use descriptive source labels so others can `ctx_search(source: "label")` later.

## ctx commands

| Command | Action |
|---------|--------|
| `ctx stats` | Call the `ctx_stats` MCP tool and display the full output verbatim |
| `ctx doctor` | Call the `ctx_doctor` MCP tool, run the returned shell command, display as checklist |
| `ctx upgrade` | Call the `ctx_upgrade` MCP tool, run the returned shell command, display as checklist |
