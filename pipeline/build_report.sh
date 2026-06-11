#!/usr/bin/env bash
# build_report.sh — render the interactive HTML run report for a finished run.
#
# Self-contained so it can run detached:  nohup bash pipeline/build_report.sh <run_dir> &
#
# Usage:
#   bash pipeline/build_report.sh <run_dir> [extra args passed to 08b_html_report.R]
#   bash pipeline/build_report.sh Results/results_ES03_newkit-ES12_newkit_filtered
#
# Output:  <run_dir>/reports/<run>_report.html   (+ build_report.log alongside it)
set -euo pipefail

RSCRIPT="${RSCRIPT:-/home/alvin/miniconda3/envs/scrna_seurat/bin/Rscript}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ $# -lt 1 ]; then
  echo "usage: bash build_report.sh <run_dir> [--samples=A,B] [--max-cells=N]" >&2
  exit 2
fi
RUN_DIR="$1"; shift
[ -d "$RUN_DIR" ] || { echo "no such run dir: $RUN_DIR" >&2; exit 2; }

LOG="$RUN_DIR/reports/build_report.log"
mkdir -p "$RUN_DIR/reports"
{
  echo "=== build_report $(date '+%F %T') ==="
  echo "run_dir: $RUN_DIR"
  echo "Rscript: $RSCRIPT"
  "$RSCRIPT" "$HERE/08b_html_report.R" "$RUN_DIR" "$@"
  echo "=== done $(date '+%F %T') ==="
} >"$LOG" 2>&1
echo "report log: $LOG"
