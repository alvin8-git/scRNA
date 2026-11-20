#!/usr/bin/env bash
# =============================================================================
# setup_extras.sh — optional packages: CellChat and SeuratWrappers for bat-wing steps 11-14,
# monocle3 for the cm trajectory, presto for faster marker tests. Not needed for the core pipeline or the
# H1/H2 smoke test. Run after setup_env.sh.
# Usage: bash pipeline/setup_extras.sh        (SCRNA_ENV_NAME picks another env)
# These build from GitHub source, so the env gets compilers and the sf/geo stack from conda.
# =============================================================================
set -euo pipefail
ENV_NAME="${SCRNA_ENV_NAME:-scrna_seurat}"

# The core pins are repeated so the solver cannot upgrade them (harmony 2.x breaks step 04).
mamba install -n "${ENV_NAME}" -c conda-forge -c bioconda -y \
    r-base=4.4.3 r-seurat=5.4.0 r-seuratobject=5.3.0 r-harmony=1.2.4 \
    compilers make cmake pkg-config hdf5 \
    r-sf r-spdep r-units r-lme4 r-nloptr r-minqa r-e1071 r-proxy r-quadprog r-ggdist \
    r-slam r-pscl r-pbmcapply r-rsample r-readr r-biglm bioconductor-batchelor

# Locate the env through conda itself (envs_dirs may not be under the base install), call its
# Rscript directly (`conda run ... Rscript - <<EOF` drops stdin, so a block fed that way silently
# does nothing), and put its bin/ on PATH as an activated env would (pandoc, compilers).
ENV_BIN="$(dirname "$(conda run -n "${ENV_NAME}" which Rscript)")"
export PATH="${ENV_BIN}:${PATH}"
"${ENV_BIN}/Rscript" - <<'REOF'
options(repos = c(CRAN = "https://cloud.r-project.org"))
gh <- c(CellChat = "jinworks/CellChat", SeuratWrappers = "satijalab/seurat-wrappers",
        presto = "immunogenomics/presto")
for (p in names(gh)) if (!requireNamespace(p, quietly = TRUE))
  remotes::install_github(gh[[p]], upgrade = "never", quiet = TRUE)
if (!requireNamespace("monocle3", quietly = TRUE))
  BiocManager::install("cole-trapnell-lab/monocle3", update = FALSE, ask = FALSE)
for (p in c(names(gh), "monocle3")) {
  library(p, character.only = TRUE, quietly = TRUE)
  cat(sprintf("  OK: %s %s\n", p, as.character(packageVersion(p))))
}
REOF
echo "Extras installed in '${ENV_NAME}'."
