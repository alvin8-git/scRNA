# HTML run report — design / QA / UX review

Rigorous, no-stone-unturned critique of `pipeline/report_template.Rmd` (+ driver
`pipeline/08b_html_report.R`). Severity: **P0** breaks something, **P1** real
defect a user will hit, **P2** polish / correctness debt, **P3** nice-to-have.
Line numbers are against the current `report_template.Rmd`.

---

## 0. What is already good (don't regress these)

- Clean design-token system (`:root` vars, one accent, consistent shadow/line).
- Explicit WCAG-AA contrast pass on small on-white labels (L202-203).
- Real responsive drawer at <820px, not just a squeeze (L206-222).
- Deliberate file-size engineering: violin downsampling to 600/sample (L25-30),
  static UMAP PNGs past 8 samples (driver L176-208), `toWebGL()` for live UMAPs.
- Honest empty-states everywhere (no doublet scores L422, no delta L485, no
  counts L462, no markers L512, no versions L602).
- Self-contained, single file, emailable. Title/date in a fixed bar. Good.

---

## 1. Navigation & information architecture

**P1 — Scrollspy mis-highlights once you reach Per-sample.** `onScroll`
(L618-631) picks the *last* nav link whose target `offsetTop <= y`. The
per-sample sections are shown one at a time; the hidden ones are `display:none`
and report `offsetTop === 0`, which is `<= y` for any scroll position. So every
hidden `#sample-*` link qualifies and the highlight jumps to the **last**
per-sample link in the rail regardless of where you actually are. Fix: skip
links whose target is hidden — `if(secs[i] && secs[i].offsetParent !== null && secs[i].offsetTop<=y) idx=i;`.

**P2 — Clicking a Per-sample nav link doesn't land on that sample.** The
galleries JS (L590) routes every `#sample-*` activation through
`goPS → hdr.scrollIntoView()`, i.e. it scrolls to the `#per-sample` *header*,
not the chosen sample's `<h2>`. Works because only one sample is visible, but a
deep-link from an email (`...#sample-es459`) lands you at the section top, one
chip-row above the content. Consider `document.getElementById('sample-'+sl).scrollIntoView()`.

**P2 — Long "Per sample" list has no affordance to collapse.** 14 samples = 14
stacked links in a 208px rail; the group labels are not collapsible. Fine now,
but it's the section most likely to grow. A `<details>` per group would scale.

**P3 — Tab state isn't reflected in nav.** Sample comparison / Cell proportions
are `{.tabset}`; the rail only knows the section, not which pill (UMAP vs QC) is
open. Acceptable, but a returning reader loses their place.

---

## 2. The bold-nav change you asked for (and a leftover)

**Done:** `.grouplabel` is now `font-weight:700` (L93). RUN / COMPARE / PER
SAMPLE now read as headers.

**P2 — your color tweak is dead.** L93 sets `color:#5b6675`, but L203 re-sets
the same selector to `#6b7686` later in the cascade (equal specificity, later
wins). The bold survives; the color does not. Net effect is fine (AA-compliant
grey, bold), but note the labels (#6b7686) and the links (`--muted` #667085) are
now nearly identical in hue — **bold + uppercase is the only separator.** If you
want them to truly pop, push grouplabels to `--ink`/#3a4250 *in the L203 rule*
(the one that actually wins), not L93.

---

## 3. Content accuracy (the report tells small lies on big runs)

**P1 — The intro over-promises interactivity for >8-sample reports.** L249-251
("drag to zoom, and click a cell type in a legend to toggle it across every
panel — that toggle *is* the depletion view") prints on *every* report. On
14-sample runs the UMAPs are static PNGs and the legend is inert
(`legitem.static`, `js_leg=""`). The L372-373 note corrects it locally, but the
headline claim above the fold is false for exactly the runs this review targets.
Make the intro conditional on `STATIC`.

**P3 — "Drag to zoom" also implies the violins/bars zoom by drag.** They do
(plotly), so this is fine — just flagging the intro sentence does double duty.

---

## 4. Dead / misleading code (correctness debt)

**P1 — The "static violins" optimization is described but never wired up.**
L42-50 define `STATIC` and a local `to_uri()` with a comment saying many-sample
runs "render the violins as static PNGs instead." But `qc_violin` (L398-413),
`doublets` (L421-433), `prop_bar`, `count_bar`, and `marker_dot` all call
`ggplotly`/`plot_ly` **unconditionally** — `STATIC` and this `to_uri` are never
read. So on 14-sample runs the QC + doublet violins are still live widgets (the
exact thing the comment claims to avoid). It's survivable only because `vcells`
caps points at 600/sample, but the code is dead and the comment is misleading.
Either implement the static path under `STATIC` or delete L45-50 and fix the
comment.

**P2 — `fig_h`, `fig_w`, `fac_ncol`, `fac_nrow` are dead** (L51-54). `qc_violin`
hardcodes `ncol = 3` and `height = 430`. Remove or use them.

---

## 5. Accessibility

**P1 — Static UMAP and gallery images have no `alt` text.** L341 (`tags$img`)
and `emit_fig` (L565) emit `<img src=...>` with no `alt`. For the largest reports
the entire UMAP comparison is images — a screen reader gets nothing. Add
`alt = sprintf("UMAP of %s", s)` and pass the caption as alt in `emit_fig`.

**P2 — Cell type is encoded by color alone in the bars.** Composition / Cell
numbers / proportions rely on a 15-color qualitative palette with several close
hues (#4e79a7 vs #76b7b2 vs #86bcb6). Colorblind readers can't separate
adjacent stacks; the only recovery is hover, which doesn't exist on touch. The
static UMAPs are fine (ggrepel text labels). Consider direct labels on the
largest stacks or a pattern fallback.

**P2 — Legend items aren't keyboard operable.** `.legitem` is a `<div>` with a
click handler (L363-366, toggle JS L371), no `role`/`tabindex`. The chips and
controls are real `<button>`s (good) but lack `aria-pressed`.

**P3 — `details`/`summary` glossary** is keyboard-fine; good.

---

## 6. Layout / visual robustness

**P1 — Long run name can overflow the fixed top bar.** `.report-top .title`
(L75) has no `max-width`/ellipsis. These runs have huge concatenated names
(`Aksh1-Aksh2-ES03-ES03_500umi-...`). On a narrow window the title will collide
with or shove the right-aligned meta chips. Add
`white-space:nowrap;overflow:hidden;text-overflow:ellipsis;max-width:48vw` and a
`title=` tooltip.

**P2 — "Unselect all" leaves a blank void.** Deselecting every chip hides all
panels; `umSolo` only special-cases `n===1`, so `n===0` shows an empty
`.umapwrap` with a floating legend and no "nothing selected" hint. Add a
placeholder line when the grid is empty.

**P3 — Sticky legend `top:64px` vs `scroll-padding-top:68px` vs 52px bar.** Minor
4px inconsistency; harmless but worth unifying on one offset constant.

---

## 7. Performance / file size

- 14-sample report: **18 MB** after the resolution bump (was 17 MB) — the
  higher-res `ragg` PNGs cost ~1 MB total, acceptable.
- The bounded violins (600/sample) keep live widgets sane; the real lever if
  size ever bites is wiring up §4's static-violin path.
- One plotly dependency is de-duped by htmlwidgets across all widgets. Good.

---

## Prioritized fix list

1. **P1** Scrollspy: ignore hidden sections (§1).
2. **P1** Conditional intro for static runs (§3).
3. **P1** `alt` text on UMAP + gallery images (§5).
4. **P1** Top-bar title ellipsis (§6).
5. **P1** Resolve the dead static-violin code: implement or delete (§4).
6. **P2** Make bold grouplabels actually darker via the winning L203 rule (§2).
7. **P2** Empty-grid placeholder for Unselect all (§6).
8. **P2** Keyboard/aria on legend + chips (§5).
9. **P2** Remove dead `fig_*`/`fac_*` (§4).
10. **P3** Deep-link to the exact sample; collapsible nav groups; offset constant.
