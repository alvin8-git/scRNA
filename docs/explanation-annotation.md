# Annotation Strategy

This document explains how step 05 assigns `final_cell_type` to every cell, why the pipeline uses two annotation routes simultaneously, and how to interpret annotation quality plots.

---

## The problem

Automated cell type annotation is hard for three reasons:

1. **Reference mismatch.** SingleR references are built from human PBMC datasets. Bat whole blood, unusual tissues, or poorly-represented cell types produce low-confidence scores.
2. **Cluster purity.** Majority-vote labelling at the cluster level works when clusters are pure, but a cluster containing 10% contaminating neutrophils is mislabelled NK if NK is the majority.
3. **Sub-type granularity.** References label cells as `"T cell"` but biology requires `"CD4 T (naive)"` vs `"CD4 T (memory)"` vs `"CD4 T (effector)"`.

---

## The two-route approach

```
Per-cell SingleR scores
        │
        ├──► Contamination-type override
        │    (Neutrophil, RBC, HSPC, Platelet, ...)
        │    Per-cell label applied regardless of cluster majority
        │
        └──► Cluster majority vote
             │
             ├──► CLUSTER_CELLTYPE_MAP override (optional, manual)
             │
             └──► Sub-type refinement (REFINEMENT_MARKERS scoring)
                         │
                         ▼
                  final_cell_type
```

### Route A: Per-cell contamination override

Cells whose SingleR label matches a type in `CONTAMINATION_TYPES` keep that per-cell label. This handles two scenarios:

- **Whole blood:** Neutrophils are the majority in several clusters. Without this override, a cluster that is 60% Neutrophil and 40% NK would be labelled NK if cluster 2 happens to have majority NK.
- **Rare populations:** A single platelet-contaminated cluster of 15 cells would be swallowed by the majority-vote label of the nearest real cluster.

Adjust `CONTAMINATION_TYPES` for your tissue: for sorted PBMC, Neutrophils are genuinely absent and should be removed from the list so they do not appear as contaminants.

### Route B: Cluster majority vote + manual override

For all non-contamination cells, the cluster's majority SingleR label wins. The `CLUSTER_CELLTYPE_MAP` lets you correct errors:

```r
CLUSTER_CELLTYPE_MAP <- c(
  "4" = "CD8 T"  # SingleR said NK, but CD3E=85% and NCAM1=0.5% says CD8 T
)
```

Any cluster not listed falls back to SingleR automatically. Partial maps are safe to use.

### Sub-type refinement

After cluster labelling, coarse labels like `"CD4 T"` or `"B cell"` are refined using marker scoring. For each cluster labelled `"CD4 T"`, the pipeline computes average expression of naive, memory, and effector markers and assigns the sub-type with the highest score:

```
Cluster 2 labelled "CD4 T":
  naive score    (CCR7, TCF7, LEF1, SELL)     → 1.8
  memory score   (IL7R, AQP3, GPR183)          → 0.9
  effector score (GZMK, TNFRSF4, CCL5)         → 0.4
  → assigned "CD4 T (naive)"
```

Set `REFINEMENT_MARKERS = NULL` to skip and keep coarse labels.

---

## Label normalisation (`SINGLER_NORM`)

Raw SingleR labels from `HumanPrimaryCellAtlasData` include entries like `"T_cells:CD4+"`, `"GMP"`, `"BFU-E"`. The 30-entry `SINGLER_NORM` map translates these to canonical pipeline names:

| Raw label | Canonical name |
|-----------|---------------|
| `"T_cells:CD4+"` | `"CD4 T"` |
| `"T_cells:CD8+"` | `"CD8 T"` |
| `"NK_cell"` | `"NK"` |
| `"BFU-E"` | `"RBC"` |
| `"GMP"` | `"Neutrophil"` |
| `"Pre-B_cell_CD34-"` | `"B cell"` |

This normalisation runs before majority vote so that `CLUSTER_CELLTYPE_MAP` keys always use canonical names.

For `MonacoImmune`, labels are already closer to canonical form — `SINGLER_NORM` applies fewer remappings.

---

## How to read annotation quality plots

### SingleR score heatmap (`singler_scores_heatmap.pdf`)

Rows = reference cell types; columns = cells (sampled subset). Bright = high score.

- **Good:** Each cell has one bright column entry and dark everywhere else — confident, unambiguous annotation.
- **Bad:** Many cells have diffuse moderate scores across multiple types — the reference cannot distinguish between types. Usually caused by reference mismatch (e.g. using `HumanPrimaryCellAtlas` for bat tissue).

### Delta score UMAP (`singler_delta_umap.pdf`)

Delta score = top SingleR score minus second-best score. High delta = confident.

- **Good:** Delta > 0.1 uniformly across cells; cluster centres bright.
- **Bad:** All cells near zero — annotation is a guess. Consider switching `SINGLER_REF` to `MonacoImmune` or adding manual overrides.

### Canonical markers dot plot (`canonical_markers_dotplot.pdf`)

Dot size = % of cells in cluster expressing the gene; colour = average expression. This is the primary tool for filling `CLUSTER_CELLTYPE_MAP`.

Look for:
- Clusters where the dot plot clearly shows CD3D/CD3E but SingleR says "NK" → override to CD8 T or CD4 T
- Clusters with high S100A8/S100A9 but SingleR says B cell → override to CD14+ Mono
- Clusters with PPBP/PF4 → Platelet (if not already caught by contamination override)

---

## When SingleR fails

SingleR accuracy degrades when:

1. **Reference is wrong species.** Use `MonacoImmune` for blood; consider a bat-specific reference for bat tissue.
2. **Library quality is poor.** Low-gene cells (<200 features) have noisy profiles. Step 01 QC removes most of these, but check the delta score UMAP for low-confidence zones.
3. **Tissue is non-immune.** Brain, liver, and skin have cell types with no reference equivalent. Expect high `"Unassigned"` rates and use marker-based annotation (`CLUSTER_CELLTYPE_MAP`) instead.

When >20% of cells are `"Unassigned"`, review the delta score UMAP and the canonical markers dot plot. Manual `CLUSTER_CELLTYPE_MAP` entries can cover gaps without rerunning SingleR.

---

## Related

- [How to Override Annotations](howto-override-annotations.md) — step-by-step workflow
- [Configuration Reference](reference-config.md) — `SINGLER_REF`, `REFINEMENT_MARKERS`, `CONTAMINATION_TYPES`
- [Pipeline Architecture](explanation-architecture.md) — why two routes exist
