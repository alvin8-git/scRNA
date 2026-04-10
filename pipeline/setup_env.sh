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
    -y

echo ""
echo "Verifying all key packages..."
conda run -n ${ENV_NAME} Rscript - <<'EOF'
pkgs <- c("Seurat", "harmony", "SingleR", "celldex", "scDblFinder",
          "ggplot2", "patchwork", "cowplot", "pheatmap", "viridis", "dplyr")
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
