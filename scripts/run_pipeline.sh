#!/usr/bin/env bash
# run_pipeline.sh
# Execute all pipeline steps in order.
# Skips steps already completed (checks .done markers).
# Usage: bash run_pipeline.sh [--force]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="$SCRIPT_DIR/results"
FORCE=${1:-""}

run_step() {
    local name="$1"
    local script="$2"
    local done_marker="$3"

    if [[ -f "$done_marker" && "$FORCE" != "--force" ]]; then
        echo "⏭  Skipping $name (already done). Use --force to re-run."
        return
    fi

    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "▶  Running: $name"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    if [[ "$script" == *.py ]]; then
        python "$SCRIPT_DIR/scripts/$script"
    elif [[ "$script" == *.R ]]; then
        Rscript "$SCRIPT_DIR/scripts/$script"
    fi
}

mkdir -p "$RESULTS_DIR"

run_step "00 Download Data"         "00_download_data.py"          "$RESULTS_DIR/00_download.done"
run_step "01 Preprocessing"         "01_preprocess.py"             "$RESULTS_DIR/01_preprocess.done"
run_step "02a BBKNN"                "02a_integrate_bbknn.py"       "$RESULTS_DIR/02a_bbknn.done"
run_step "02b Harmony"              "02b_integrate_harmony.py"     "$RESULTS_DIR/02b_harmony.done"
run_step "02c Seurat CCA"           "02c_integrate_seurat_wrap.py" "$RESULTS_DIR/02c_seurat.done"
run_step "02d Scanorama"            "02d_integrate_scanorama.py"   "$RESULTS_DIR/02d_scanorama.done"
run_step "03 Benchmarking (scIB)"   "03_benchmark.py"              "$RESULTS_DIR/03_benchmark.done"
run_step "04 Visualization"         "04_visualize.py"              "$RESULTS_DIR/04_visualize.done"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "✅  Pipeline complete."
echo "   Report: report/summary.html"
echo "   Figures: results/figures/"
echo "   Metrics: results/metrics/scores.csv"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
