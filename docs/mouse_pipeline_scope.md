# Mouse cardiomyocyte pipeline — scope

Two DNBelab flowcells of mouse cardiomyocytes on the LynnYap mount. Mouse is
**not** one of the validated sample types (see `report_product_design.md` §0), so
this is new work, not a new keyword.

Status: scope only, nothing implemented. Written 2026-10-24. Every fact below
was checked against the data, not assumed.

## 1. The data

`/home/alvin/storeData/DNBelab_Data/LynnYap/`

| flowcell | samples |
|---|---|
| `E250164106_POS_A` | 0hCell_1, 0hCell_2, 0hMed_1, 0hMed_2, 6hCell_1, 6hCell_2, 6hMed_1, 6hMed_2, **MS_Sham1, MS_Sham2** |
| `E250164093_POS_B` | 24hCell_1, 24hCell_2, 24hMed_1, 24hMed_2, 72hCell_1, 72hCell_2, **72hMed_a, 72hMed_b** |

**18 samples**, not the 12 a "timecourse × 2 arms × 2 replicates" reading
suggests. Layout per sample is `outs/filter_matrix` + `outs/raw_matrix`,
identical to the cardiac samples, so copying and SoupX both work unchanged.

Design, as far as the names reveal it: a timecourse (0h, 6h, 24h, 72h) with two
arms (`Cell` and `Med`) in duplicate, plus two `MS_Sham` samples on PosA only.

**Design as clarified by Alvin, 2026-10-25:** `Med` is medium, the vehicle /
negative control for its timepoint; `Cell` is the treated arm. `MS_Sham` is most
likely a control set from an earlier experiment sharing the flowcell, and is
excluded from the timecourse runs. The analysis is therefore two-factor
(arm × time), and the primary contrast is Cell − Med *within* a timepoint, which
is unconfounded by flowcell because both arms of a timepoint were sequenced
together. Time alone remains confounded with flowcell (§1 item 1).

**Step 01 on all 18 (2026-10-25, `Results/results_0hCell1_18samples_filtered`):**
median mito 0.4–1.6 % in every sample, 95th percentile under 10 %, so
`max_percent_mt` 30 is inert and can stay. `72hMed_a` is a failed capture (1 cell
passes QC of 1,434 called barcodes); `72hMed_b` has 60 cells; `72hCell_2` 241. The
72 h timepoint has one usable library (`72hCell_1`, 426 cells) and no Med control.
First end-to-end run: PosA, 0 h + 6 h, both arms, without the shams.

### Three things about the design that need answers before analysis

1. **Flowcell is confounded with timepoint.** 0h and 6h are on PosA, 24h and 72h
   on PosB. Any 0h-vs-72h difference is inseparable from a flowcell batch
   effect. Harmony will integrate them, but integration cannot distinguish a
   batch effect from a biological one when the two are perfectly aligned — it
   will remove both. This is the single most important caveat for this dataset
   and it is a property of how it was sequenced, not something analysis fixes.
2. **`MS_Sham1`/`MS_Sham2` are unexplained** and exist only on PosA. Sham
   surgery controls? A different experiment sharing the flowcell? They should
   not be pooled with the timecourse until that is known.
3. **`72hMed_a` / `72hMed_b` break the replicate naming.** Every other sample
   ends `_1`/`_2`. The report's sample-grouping derives condition by stripping
   a trailing `_[0-9]+`, so these two would each become their own group and the
   change-over-time view would silently show six groups instead of four.

## 2. The blocker: mitochondrial and ribosomal patterns are hardcoded to human

`pipeline/01_load_qc.R:55-56`

```r
seu[["percent.mt"]]   <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
```

Measured in the mouse features file (`0hCell_1/outs/filter_matrix`):

| pattern | genes matched |
|---|---|
| `^mt-` | **13** |
| `^MT-` | **0** |

So on mouse data `percent.mt` is **0.00 for every cell**. Consequences, all
silent — no error, no warning:

* `QC$max_percent_mt` filtering becomes a no-op; dying and stressed cells are
  never removed.
* The QC violin renders a flat line at zero and looks like a clean dataset.
* `percent.ribo` is likewise empty (mouse is `Rps2`, `Rpl3`).

**This must be fixed before a single mouse sample is analysed**, and it is the
kind of fault that produces a confident, publishable-looking, wrong result. The
fix is to make the pattern species-aware in config rather than literal in the
step: a `QC$mt_pattern` / `QC$ribo_pattern` pair defaulting to the human form
and overridden in the mouse overlay. That also makes the assumption visible
instead of buried in step 01.

Worth re-checking for bat at the same time: bat gene symbols come from the
*Eonycteris spelaea* annotation, and whether they are uppercase human-style has
never been verified here. If they are not, every bat run to date has a
`percent.mt` of zero too.

## 3. Species profile: what a mouse overlay needs

A new `pipeline/config_species_mouse.R`, sourced the same way as the bat and cm
overlays, under `if (.species == "mouse_cm")`.

### Reference — this is the good news

`celldex::MouseRNAseqData()` was fetched and inspected. Its `label.main` values:

```
Adipocytes, Astrocytes, B cells, Cardiomyocytes, Dendritic cells,
Endothelial cells, Epithelial cells, Erythrocytes, Fibroblasts, Granulocytes,
Hepatocytes, Macrophages, Microglia, Monocytes, Neurons, NK cells,
Oligodendrocytes, T cells
```

**It contains `Cardiomyocytes`.** Gene symbols are title-case (`Xkr4`, `Sox17`),
matching the data.

This is the opposite of the human cardiac situation, where
`HumanPrimaryCellAtlas` has no cardiomyocyte class and SingleR therefore calls
cardiac cells "Neurons". Mouse can get a correct cardiomyocyte call from the
reference out of the box. `ImmGenData` is immune-only and not useful here.

### Marker panels

Mouse symbols are title-case, so **not one gene** from the human `cm` panels
matches. `Tnnt2` is present in the data; `TNNT2` returns zero. The panel must be
rewritten, not translated by uppercasing:

* cardiomyocyte: `Tnnt2, Myh6, Myh7, Actc1, Myl7, Myl4, Ttn, Actn2, Nkx2-5, Ryr2`
* cardiac progenitor: `Nkx2-5, Isl1, Tbx5, Gata4, Mesp1, Kdr`
* fibroblast: `Col1a1, Col3a1, Dcn, Lum, Postn, Pdgfra, Tcf21`
* endothelial: `Pecam1, Cdh5, Kdr, Tie1`
* immune: `Ptprc, Cd68, Adgre1, Lyz2, Cd3e`
* smooth muscle / pericyte: `Acta2, Tagln, Myh11, Pdgfrb, Rgs5`
* proliferating: `Mki67, Top2a`

Note `Nkx2-5` keeps its hyphen in mouse as in human, but `Lyz2` replaces `LYZ`
and `Adgre1` replaces `ADGRE1`/F4/80 — these are not mechanical case changes and
each needs checking against the actual features file rather than assumed.

### QC thresholds

Cardiomyocytes are mitochondria-dense by function, so a mito cutoff tuned for
PBMCs will discard healthy cells. The `cm` overlay already raises
`max_features` to 9000 and `max_counts` to 70000; mouse CMs need the same
treatment plus a deliberately higher `max_percent_mt`. **The right value cannot
be chosen until `percent.mt` actually computes** — picking one now would be
guessing. Run step 01 with the pattern fixed, read the distribution, then set it.

## 4. Proposed sequence

1. **Make mito/ribo patterns species-aware** in `config.R` + `01_load_qc.R`, with
   the human default unchanged. Verify the existing human and bat runs are
   unaffected (a regression guard in `tests/test_regressions.R` belongs here).
2. **Answer the three design questions** in §1 — especially what `MS_Sham` is.
3. **Copy the 18 samples** into `Samples/Mouse/` with the existing verifier
   (per-file `stat -c %s`, not `du`).
4. **Write `config_species_mouse.R`** with panels checked gene-by-gene against
   the features file, `SINGLER_REF <- "MouseRNAseqData"`, and QC bounds left at
   cm-like values pending step 01.
5. **Run step 01 only**, on two samples, and read the `percent.mt` distribution.
   Set `max_percent_mt` from the data. This is a checkpoint, not a formality.
6. **Run one flowcell end to end** (PosA, 10 samples) before attempting both.
   18 samples exceeds the 4-sample run-dir collapse threshold and the
   `umap_static` switch at 8, so the reports will take the static path.
7. **Only then** consider the combined 18-sample cohort, and state the
   flowcell/timepoint confound prominently wherever its results appear.

## 5. What this does not need

* No new report code. `08b --customer`, `08e` and `08f` all take their panel,
  palette and grouping from config, so they work once the overlay exists — with
  one exception: the grouping regex needs to handle `_a`/`_b` as well as
  `_[0-9]+`, or `72hMed_a` and `72hMed_b` become separate groups.
* No SoupX changes. `soupx_adjust.R` is species-agnostic and the mouse samples
  have `raw_matrix`, so ambient correction works as-is.
* No new palette. `CELLTYPE_COLORS` already covers Cardiomyocyte, Fibroblast,
  Endothelial, Smooth Muscle, Proliferating and Unknown; only mouse-specific
  labels absent from it would need adding, and `validate_config.R` hard-errors
  on any that are missing — so that check is automatic.
