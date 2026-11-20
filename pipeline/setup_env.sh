#!/usr/bin/env bash
# =============================================================================
# setup_env.sh — Create conda environment for scRNA-seq Seurat pipeline
# Usage: bash setup_env.sh
# =============================================================================
set -euo pipefail

ENV_NAME="${SCRNA_ENV_NAME:-scrna_seurat}"   # run_pipeline.sh reads the same variable

echo "Creating conda environment: ${ENV_NAME}"
echo "This may take 10–20 minutes..."

# R 4.4 is a hard floor, not a preference: config.R:17 calls `%||%` inside its
# source() of pdf_helpers.R, 471 lines before config.R defines the operator itself.
# base R only gained `%||%` in 4.4.0, so on 4.3.x that line errors and nothing in the
# pipeline loads. Do not lower r-base below 4.4 without first hoisting a `%||%`
# definition above config.R:17. Pins below match the verified working env.
# harmony is pinned: 2.x renamed max.iter.harmony to max_iter, which fails step 04.
mamba create -n ${ENV_NAME} \
    -c conda-forge -c bioconda \
    r-base=4.4.3 \
    r-seurat=5.4.0 \
    r-seuratobject=5.3.0 \
    r-harmony=1.2.4 \
    bioconductor-singler \
    bioconductor-celldex \
    bioconductor-scdblfinder \
    bioconductor-biocparallel \
    r-ggplot2 \
    r-patchwork \
    r-cowplot \
    r-pheatmap \
    r-viridis \
    r-dplyr \
    r-remotes \
    r-biocmanager \
    r-matrix \
    r-hdf5r \
    r-ggrepel \
    r-pdftools \
    r-magick \
    bioconductor-deseq2 \
    bioconductor-clusterprofiler \
    bioconductor-enrichplot \
    bioconductor-org.hs.eg.db \
    bioconductor-fgsea \
    r-igraph \
    r-nnls \
    r-circlize \
    r-soupx \
    bioconductor-edger \
    bioconductor-mast \
    r-leidenbase \
    r-rmarkdown \
    r-knitr \
    r-dt \
    r-plotly \
    r-crosstalk \
    r-ragg \
    pandoc \
    -y

# Locate the env through conda itself (envs_dirs may not be under the base install), call its
# Rscript directly (`conda run ... Rscript - <<EOF` drops stdin, so a block fed that way silently
# does nothing), and put its bin/ on PATH as an activated env would (pandoc, compilers).
ENV_BIN="$(dirname "$(conda run -n "${ENV_NAME}" which Rscript)")"
export PATH="${ENV_BIN}:${PATH}"
RS="${ENV_BIN}/Rscript"

echo ""
echo "Verifying all key packages..."
"$RS" - <<'EOF'
pkgs <- c("Seurat", "harmony", "SingleR", "celldex", "scDblFinder",
          "DESeq2", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "fgsea",
          "SoupX", "edgeR", "MAST", "leidenbase", "rmarkdown", "knitr", "DT", "plotly", "ragg",
          "ggplot2", "patchwork", "cowplot", "pheatmap", "viridis", "dplyr",
          "pdftools", "magick")
for (p in pkgs) {
  library(p, character.only = TRUE, quietly = TRUE)
  cat(sprintf("  OK: %s %s\n", p, as.character(packageVersion(p))))
}
stopifnot(rmarkdown::pandoc_available())
cat("\nAll packages loaded successfully.\n")
EOF

echo ""
echo "Environment '${ENV_NAME}' is ready."
echo "Activate with: conda activate ${ENV_NAME}"
echo "Run the smoke test with: bash pipeline/run_pipeline.sh Samples/H1 Samples/H2"
echo "Bat-wing steps 11-14 and the cm trajectory also need: bash pipeline/setup_extras.sh"
