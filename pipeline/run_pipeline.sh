#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — Master runner for the scRNA-seq Seurat pipeline.
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
#   # Species flag (default: human) — applies bat-specific marker/QC overrides:
#   bash run_pipeline.sh bat /path/to/ES03 /path/to/ES12
#   bash run_pipeline.sh human /path/to/H1 /path/to/H2
#
#   # Specific steps only:
#   bash run_pipeline.sh bat /path/to/ES03 /path/to/ES12 01 02 03
#   bash run_pipeline.sh 05 06 07
#
# Outputs:
#   results/report_doublets.pdf
#   results/report_individual.pdf
#   results/report_annotation.pdf
#   results/report_integrated.pdf
# =============================================================================
set -euo pipefail

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  sed -n '2,/^set -euo pipefail/{/^set -euo pipefail/!p}' "${BASH_SOURCE[0]}"
  exit 0
fi

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
SPECIES="human"   # default; overridden by 'bat' or 'human' keyword in args

for arg in "$@"; do
  if [[ "$arg" == "bat" || "$arg" == "human" || "$arg" == "bat_wing" ]]; then
    SPECIES="$arg"
  elif [[ "$arg" == condition=* ]]; then
    export SCRNA_CONDITION="${arg#condition=}"
  # Treat as sample path if absolute, relative directory, or starts with ./ or ../
  elif [[ "$arg" == /* ]] || [[ -d "$arg" ]] || [[ "$arg" == ./* ]] || [[ "$arg" == ../* ]]; then
    SAMPLE_PATHS+=("$(realpath "$arg")")
  else
    STEPS_ARGS+=("$arg")
  fi
done

export SCRNA_SPECIES="$SPECIES"

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
  if [[ "$bn" == "filter_matrix" || "$bn" == "filtered_feature_bc_matrix" ||
        "$bn" == "raw_matrix"    || "$bn" == "raw_feature_bc_matrix" ]]; then
    basename "$(dirname "$1")"
  else
    echo "$bn"
  fi
}
_matrix_tag() {
  local bn; bn=$(basename "$1")
  if [[ "$bn" == "raw_matrix" || "$bn" == "raw_feature_bc_matrix" ]]; then
    echo "raw"
  else
    echo "filtered"
  fi
}
if [[ $N_SAMPLES -gt 0 ]]; then
  _parts=(); _any_raw=false
  for p in "${SAMPLE_PATHS[@]}"; do
    _parts+=("$(_sample_name "$p")")
    [[ "$(_matrix_tag "$p")" == "raw" ]] && _any_raw=true
  done
  _suffix="${_parts[0]}"
  for (( _j=1; _j<${#_parts[@]}; _j++ )); do _suffix="${_suffix}-${_parts[$_j]}"; done
  $_any_raw && _mtag="raw" || _mtag="filtered"
  RESULTS_DIR="${BASE_DIR}/Results/results_${_suffix}_${_mtag}"
else
  RESULTS_DIR="${BASE_DIR}/Results/results_H1-H2_filtered"  # matches hardcoded config.R defaults
fi
LOG_DIR="${RESULTS_DIR}/logs"
mkdir -p "${LOG_DIR}"

# =============================================================================
# Ensure conda environment is active — auto-activate if needed
# =============================================================================
if [[ "${CONDA_DEFAULT_ENV:-}" != "${CONDA_ENV}" ]]; then
  log "Conda env '${CONDA_ENV}' not active — attempting to activate..."

  # Locate conda's shell init script
  _conda_sh=""
  _conda_base="$(conda info --base 2>/dev/null)"
  [[ -f "${_conda_base}/etc/profile.d/conda.sh" ]] && \
    _conda_sh="${_conda_base}/etc/profile.d/conda.sh"

  # Fallback: common install locations
  if [[ -z "$_conda_sh" ]]; then
    for _c in "${HOME}/miniconda3" "${HOME}/miniforge3" \
               "${HOME}/anaconda3" "/opt/conda" "/usr/local/anaconda3"; do
      if [[ -f "${_c}/etc/profile.d/conda.sh" ]]; then
        _conda_sh="${_c}/etc/profile.d/conda.sh"
        break
      fi
    done
  fi

  if [[ -z "$_conda_sh" ]]; then
    err "Cannot find conda init script. Activate manually:"
    err "  conda activate ${CONDA_ENV}"
    exit 1
  fi

  # shellcheck source=/dev/null
  source "$_conda_sh"

  if ! conda activate "${CONDA_ENV}" 2>/dev/null; then
    err "Environment '${CONDA_ENV}' not found."
    err "Create it with:  bash ${PIPELINE_DIR}/setup_env.sh"
    exit 1
  fi
  log "Activated conda env: ${CONDA_ENV}"
else
  log "Conda env '${CONDA_ENV}' already active."
fi

# =============================================================================
# Pin BLAS/OMP to 1 thread per process — prevents oversubscription when
# future/BiocParallel workers each inherit the full BLAS thread pool.
# =============================================================================
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLAS_NUM_THREADS=1

# =============================================================================
# Validate config before running any steps
# =============================================================================
echo "[run_pipeline] Validating config..."
conda run -n "${CONDA_ENV}" Rscript "${PIPELINE_DIR}/validate_config.R" || exit 1

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
  ["06b"]="Differential Expression"
  ["07"]="Finalise PDF Reports"
  ["08"]="Comparison Report"
  ["11"]="Wing DEGs & Module Scores"
  ["12"]="Pathway Enrichment"
  ["13"]="CellChat Communication"
  ["14"]="Trajectory Analysis"
)

declare -A STEP_SCRIPT=(
  ["01"]="01_load_qc.R"
  ["02"]="02_doublets.R"
  ["03"]="03_individual.R"
  ["04"]="04_integrate.R"
  ["05"]="05_annotate.R"
  ["06"]="06_visualize.R"
  ["06b"]="06b_differential.R"
  ["07"]="07_finalize_reports.R"
  ["08"]="08_comparison_report.R"
  ["11"]="11_wing_degs.R"
  ["12"]="12_pathways.R"
  ["13"]="13_cellchat.R"
  ["14"]="14_trajectory.R"
)

# =============================================================================
# Determine which steps to run
# =============================================================================
if [[ ${#STEPS_ARGS[@]} -gt 0 ]]; then
  STEPS=("${STEPS_ARGS[@]}")
elif [[ "$SPECIES" == "bat_wing" ]]; then
  STEPS=("01" "02" "03" "04" "05" "06" "06b" "07" "11" "12" "13" "14")
elif [[ "$SINGLE_SAMPLE" == "true" ]]; then
  STEPS=("01" "02" "03" "04" "05" "06" "07")
else
  STEPS=("01" "02" "03" "04" "05" "06" "06b" "07")
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
echo "  Species: ${SPECIES}"
echo "  Steps: ${STEPS[*]}"
echo "=============================================="
echo ""

TOTAL_START=$(date +%s)

for step in "${STEPS[@]}"; do
  if [[ -z "${STEP_SCRIPT[$step]+_}" ]]; then
    err "Unknown step: ${step}. Valid: 01 02 03 04 05 06 06b 07 08 11 12 13 14"
    exit 1
  fi
  run_step "${step}"
done

TOTAL_ELAPSED=$(( $(date +%s) - TOTAL_START ))
step_hdr "Pipeline complete in ${TOTAL_ELAPSED}s"
echo -e "  ${GREEN}Results :${NC} ${RESULTS_DIR}"
echo ""
echo -e "  ${CYAN}Final reports:${NC}"
echo "    ${RESULTS_DIR}/01-QC_report.pdf"
echo "    ${RESULTS_DIR}/02-Doublet_report.pdf"
echo "    ${RESULTS_DIR}/03-Individual_report.pdf"
echo "    ${RESULTS_DIR}/04-Annotation_report.pdf"
echo "    ${RESULTS_DIR}/05-Integrated_report.pdf"
echo "    ${RESULTS_DIR}/Overall_report.pdf"
echo ""
if [[ "$SPECIES" == "bat_wing" ]]; then
  echo ""
  echo -e "  ${CYAN}Wing tissue reports:${NC}"
  echo "    ${RESULTS_DIR}/differential/  (DEGs, volcano plots, module scores)"
  echo "    ${RESULTS_DIR}/pathways/      (GO/KEGG/GSEA per cell type)"
  echo "    ${RESULTS_DIR}/cellchat/cellchat_report.pdf"
  echo "    ${RESULTS_DIR}/trajectory/trajectory_report.pdf"
fi
echo -e "  ${YELLOW}Next step (manual annotation, if needed):${NC}"
echo "    1. Open ${RESULTS_DIR}/annotation/canonical_markers_dotplot.pdf"
echo "    2. Fill CLUSTER_CELLTYPE_MAP in pipeline/config.R"
echo "    3. Re-run:  bash run_pipeline.sh 05 06 07"
echo -e "${CYAN}${BOLD}+$(printf '=%.0s' {1..54})+${NC}"
