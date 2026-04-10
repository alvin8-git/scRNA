#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — Master runner for the PBMC scRNA-seq Seurat pipeline.
#
# Usage:
#   # Run all steps on hardcoded defaults in config.R:
#   bash run_pipeline.sh
#
#   # Single sample (auto-detects filter_matrix subfolder):
#   bash run_pipeline.sh /path/to/sampleA
#
#   # Two samples → Harmony integration enabled:
#   bash run_pipeline.sh /path/to/sampleA /path/to/sampleB
#
#   # Specific steps only:
#   bash run_pipeline.sh /path/to/sampleA 01 02 03
#   bash run_pipeline.sh 05 06 07
#
# Outputs:
#   results/report_doublets.pdf
#   results/report_individual.pdf
#   results/report_annotation.pdf
#   results/report_integrated.pdf
# =============================================================================
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="/data/alvin/scRNA"
CONDA_ENV="scrna_seurat"

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'
log()      { echo -e "${GREEN}[$(date '+%H:%M:%S')]${NC} $*"; }
warn()     { echo -e "${YELLOW}[$(date '+%H:%M:%S')] WARN:${NC} $*"; }
err()      { echo -e "${RED}[$(date '+%H:%M:%S')] ERROR:${NC} $*" >&2; }
step_hdr() { echo -e "\n${CYAN}${BOLD}+$(printf '=%.0s' {1..54})+${NC}"
             echo -e "${CYAN}${BOLD}|  $1${NC}"
             echo -e "${CYAN}${BOLD}+$(printf '=%.0s' {1..54})+${NC}"; }

# =============================================================================
# Parse arguments — absolute paths are sample folders; rest are step numbers
# Supports any number of samples: path1 path2 path3 ... [steps]
# =============================================================================
SAMPLE_PATHS=()
STEPS_ARGS=()

for arg in "$@"; do
  if [[ "$arg" == /* ]]; then
    SAMPLE_PATHS+=("$arg")
  else
    STEPS_ARGS+=("$arg")
  fi
done

# Export as SCRNA_SAMPLE1, SCRNA_SAMPLE2, ... so R config.R picks them up
for idx in "${!SAMPLE_PATHS[@]}"; do
  export "SCRNA_SAMPLE$((idx + 1))=${SAMPLE_PATHS[$idx]}"
done

N_SAMPLES=${#SAMPLE_PATHS[@]}
if [[ $N_SAMPLES -eq 1 ]]; then
  SINGLE_SAMPLE=true
else
  SINGLE_SAMPLE=false
fi

# Compute results dir name (mirrors config.R logic) for LOG_DIR
# Derives sample name from path: uses parent of filter_matrix, else basename
_sample_name() {
  local bn; bn=$(basename "$1")
  if [[ "$bn" == "filter_matrix" || "$bn" == "filtered_feature_bc_matrix" ]]; then
    basename "$(dirname "$1")"
  else
    echo "$bn"
  fi
}
if [[ $N_SAMPLES -gt 0 ]]; then
  _parts=()
  for p in "${SAMPLE_PATHS[@]}"; do _parts+=("$(_sample_name "$p")"); done
  _suffix="${_parts[0]}"
  for (( _j=1; _j<${#_parts[@]}; _j++ )); do _suffix="${_suffix}and${_parts[$_j]}"; done
  RESULTS_DIR="${BASE_DIR}/results_${_suffix}"
else
  RESULTS_DIR="${BASE_DIR}/results_H1andH2"  # matches hardcoded config.R defaults
fi
LOG_DIR="${RESULTS_DIR}/logs"
mkdir -p "${LOG_DIR}"

# =============================================================================
# Check conda environment
# =============================================================================
if [[ "${CONDA_DEFAULT_ENV:-}" != "${CONDA_ENV}" ]]; then
  warn "Conda environment '${CONDA_ENV}' is not active."
  warn "Run: conda activate ${CONDA_ENV}"
  echo ""
  read -rp "Continue anyway? [y/N] " yn
  [[ "${yn,,}" == "y" ]] || exit 1
fi

# =============================================================================
# Step definitions
# =============================================================================
declare -A STEP_NAME=(
  ["01"]="Load & QC"
  ["02"]="Doublet Detection"
  ["03"]="Individual Analysis"
  ["04"]="Integration (Harmony)"
  ["05"]="Cell Type Annotation"
  ["06"]="Visualisation"
  ["07"]="Finalise PDF Reports"
)

declare -A STEP_SCRIPT=(
  ["01"]="01_load_qc.R"
  ["02"]="02_doublets.R"
  ["03"]="03_individual.R"
  ["04"]="04_integrate.R"
  ["05"]="05_annotate.R"
  ["06"]="06_visualize.R"
  ["07"]="07_finalize_reports.R"
)

# =============================================================================
# Determine which steps to run
# =============================================================================
if [[ ${#STEPS_ARGS[@]} -gt 0 ]]; then
  STEPS=("${STEPS_ARGS[@]}")
elif [[ "$SINGLE_SAMPLE" == "true" ]]; then
  STEPS=("01" "02" "03" "04" "05" "06" "07")
else
  STEPS=("01" "02" "03" "04" "05" "06" "07")
fi

# =============================================================================
# Run a single step
# =============================================================================
STEP_NUM=0
run_step() {
  local step="$1"
  local script="${STEP_SCRIPT[$step]}"
  local name="${STEP_NAME[$step]}"
  local logfile="${LOG_DIR}/${step}_${script%.R}.log"
  local start_ts
  start_ts=$(date +%s)
  STEP_NUM=$(( STEP_NUM + 1 ))

  step_hdr "Step ${step} [${STEP_NUM}/${#STEPS[@]}]  ${name}"
  echo -e "  ${CYAN}Script :${NC} ${script}"
  echo -e "  ${CYAN}Log    :${NC} ${logfile}"
  echo ""

  if Rscript "${PIPELINE_DIR}/${script}" 2>&1 | tee "${logfile}"; then
    local elapsed=$(( $(date +%s) - start_ts ))
    echo -e "\n  ${GREEN}DONE${NC}  Step ${step} finished in ${elapsed}s"
  else
    echo -e "\n  ${RED}FAILED${NC}  Step ${step} — check log: ${logfile}"
    exit 1
  fi
}

# =============================================================================
# Main
# =============================================================================
echo "=============================================="
echo "  scRNA-seq Pipeline — Seurat + Harmony"
echo "  $(date '+%Y-%m-%d %H:%M:%S')"
if [[ $N_SAMPLES -gt 0 ]]; then
  for idx in "${!SAMPLE_PATHS[@]}"; do
    echo "  Sample $((idx+1)): ${SAMPLE_PATHS[$idx]}"
  done
  [[ "$SINGLE_SAMPLE" == "true" ]]  && echo "  Mode: single-sample (no integration)"
  [[ "$SINGLE_SAMPLE" == "false" ]] && echo "  Mode: ${N_SAMPLES}-sample (Harmony integration)"
else
  echo "  Samples: from config.R defaults"
fi
echo "  Steps: ${STEPS[*]}"
echo "=============================================="
echo ""

TOTAL_START=$(date +%s)

for step in "${STEPS[@]}"; do
  if [[ -z "${STEP_SCRIPT[$step]+_}" ]]; then
    err "Unknown step: ${step}. Valid: 01 02 03 04 05 06 07"
    exit 1
  fi
  run_step "${step}"
done

TOTAL_ELAPSED=$(( $(date +%s) - TOTAL_START ))
step_hdr "Pipeline complete in ${TOTAL_ELAPSED}s"
echo -e "  ${GREEN}Results :${NC} ${RESULTS_DIR}"
echo ""
echo -e "  ${CYAN}Final reports:${NC}"
echo "    ${RESULTS_DIR}/report_qc.pdf"
echo "    ${RESULTS_DIR}/report_doublets.pdf"
echo "    ${RESULTS_DIR}/report_individual.pdf"
echo "    ${RESULTS_DIR}/report_annotation.pdf"
echo "    ${RESULTS_DIR}/report_integrated.pdf"
echo ""
echo -e "  ${YELLOW}Next step (manual annotation, if needed):${NC}"
echo "    1. Open ${RESULTS_DIR}/annotation/canonical_markers_dotplot.pdf"
echo "    2. Fill CLUSTER_CELLTYPE_MAP in pipeline/config.R"
echo "    3. Re-run:  bash run_pipeline.sh 05 06 07"
echo -e "${CYAN}${BOLD}+$(printf '=%.0s' {1..54})+${NC}"
