#!/usr/bin/env bash
# =============================================================================
# setup_env.sh — Create conda environment for scRNA-seq Seurat pipeline
# Usage: bash setup_env.sh
# =============================================================================
set -euo pipefail

ENV_NAME="scrna_seurat"

echo "Creating conda environment: ${ENV_NAME}"
echo "This may take 10–20 minutes..."

mamba create -n ${ENV_NAME} \
    -c conda-forge -c bioconda \
    r-base=4.3.3 \
    r-seurat=5.1.0 \
    r-harmony \
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
    -y

echo ""
echo "Installing CellChat v2 and Monocle3 via R (not on conda)..."
conda run -n ${ENV_NAME} Rscript - <<'REOF'
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("CellChat", quietly = TRUE)) {
  remotes::install_github("jinworks/CellChat", upgrade = "never", quiet = TRUE)
}
if (!requireNamespace("monocle3", quietly = TRUE)) {
  BiocManager::install("cole-trapnell-lab/monocle3", update = FALSE, ask = FALSE)
}
if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
  remotes::install_github("satijalab/seurat-wrappers", upgrade = "never", quiet = TRUE)
}
cat("CellChat:", as.character(packageVersion("CellChat")), "\n")
cat("monocle3:", as.character(packageVersion("monocle3")), "\n")
cat("SeuratWrappers:", as.character(packageVersion("SeuratWrappers")), "\n")
REOF

echo ""
echo "Verifying all key packages..."
conda run -n ${ENV_NAME} Rscript - <<'EOF'
pkgs <- c("Seurat", "harmony", "SingleR", "celldex", "scDblFinder",
          "DESeq2", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "fgsea",
          "CellChat", "monocle3", "SeuratWrappers",
          "ggplot2", "patchwork", "cowplot", "pheatmap", "viridis", "dplyr",
          "pdftools", "magick")
for (p in pkgs) {
  library(p, character.only = TRUE, quietly = TRUE)
  cat(sprintf("  OK: %s %s\n", p, as.character(packageVersion(p))))
}
cat("\nAll packages loaded successfully.\n")
EOF

echo ""
echo "Environment '${ENV_NAME}' is ready."
echo "Activate with: conda activate ${ENV_NAME}"
echo "Run pipeline with: bash run_pipeline.sh"
