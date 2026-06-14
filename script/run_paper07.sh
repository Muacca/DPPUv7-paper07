#!/usr/bin/env bash
# run_paper07.sh
# Run all paper07 proof witness scripts (scripts/paper07/).
# Each script writes its own timestamped log to the paper-level data/ directory.
#
# Usage (from the script/ directory):
#   bash run_paper07.sh                          # run all -> ../data/
#   bash run_paper07.sh --output-dir /path/data  # custom log directory
#
# Output: <output-dir>/<script_name>_<timestamp>.log  (written by each script)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${SCRIPT_DIR}/scripts"

# ── argument parsing ──────────────────────────────────────────────────────────

TARGET="all"
DATA_DIR="${SCRIPT_DIR}/../data"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir|-o)
            DATA_DIR="$2"
            shift 2
            ;;
        all|paper07)
            TARGET="$1"
            shift
            ;;
        *)
            echo "Usage: bash run_paper07.sh [--output-dir DIR] [all|paper07]" >&2
            exit 1
            ;;
    esac
done

mkdir -p "${DATA_DIR}"

export PYTHONIOENCODING=utf-8
export PYTHONUTF8=1
export DPPU_LOG_DIR="${DATA_DIR}"
export DPPU_LOG_STDOUT=0

# ── helper ────────────────────────────────────────────────────────────────────

FAILED_SCRIPTS=()

run_script() {
    local rel_path="$1"
    local script_file="${SCRIPTS_DIR}/${rel_path}"
    local base
    base="$(basename "${rel_path}" .py)"
    echo ">>> ${rel_path}"
    local t0
    t0=$(date +%s)

    if python "${script_file}"; then
        local elapsed=$(( $(date +%s) - t0 ))
        echo "    OK  (${elapsed}s)  -> ${DATA_DIR}/${base}_<timestamp>.log"
    else
        local elapsed=$(( $(date +%s) - t0 ))
        echo "    FAILED  (${elapsed}s)  -> ${DATA_DIR}/${base}_<timestamp>.log"
        FAILED_SCRIPTS+=("${rel_path}")
    fi
}

# ── script lists ──────────────────────────────────────────────────────────────

# paper07 proof witness scripts (proposal order; excludes internal QA utilities)
PAPER07_SCRIPTS=(
    paper07/rank_defect_activation_witness.py
    paper07/spectral_completion_witness.py
    paper07/torsion_auxiliary_field_promotion_witness.py
    paper07/finite_fate_branch_witness.py
    paper07/harmonic_kernel_correspondence.py
    paper07/route_comparison_s3_branch_a_witness.py
    paper07/s3_normalization_factor_witness.py
)

# ── dispatch ──────────────────────────────────────────────────────────────────

echo "========================================"
echo "  run_paper07.sh  (target: ${TARGET})"
echo "  output dir: ${DATA_DIR}"
echo "========================================"
echo ""

if [[ "${TARGET}" == "all" || "${TARGET}" == "paper07" ]]; then
    echo "--- paper07/ ---"
    for s in "${PAPER07_SCRIPTS[@]}"; do
        run_script "${s}"
    done
    echo ""
fi

# ── summary ───────────────────────────────────────────────────────────────────

echo "========================================"
if [[ ${#FAILED_SCRIPTS[@]} -eq 0 ]]; then
    echo "  All scripts completed successfully."
else
    echo "  ${#FAILED_SCRIPTS[@]} script(s) failed:"
    for s in "${FAILED_SCRIPTS[@]}"; do
        echo "    - ${s}"
    done
fi
echo "========================================"

if [[ ${#FAILED_SCRIPTS[@]} -gt 0 ]]; then
    exit 1
fi