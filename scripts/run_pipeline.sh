#!/usr/bin/env bash
# End-to-end BiologicTraj pipeline on one dataset.
#
#   1. processing-stress check            (Python, scanpy)
#   2. summed-count pseudobulk            (Python)
#   3. edgeR/limma-voom count model       (R)  -> logCPM + transition tests
#   4. LongOmicsTraj trajectory shapes    (R)  -> empirical + LMM, several thresholds
#   5. score against ground truth         (Python; synthetic data only)
#
# Usage: scripts/run_pipeline.sh <h5ad> <out_dir> [truth_csv] [estimators]
#   estimators: "empirical" (fast), "lmm" or "empirical,lmm" (default; the LMM
#   takes ~0.1 s per gene per cell type, ~30 min for the synthetic data on one core)
#   e.g. scripts/run_pipeline.sh data/synthetic/prism_like_synthetic.h5ad results/synthetic/drug \
#          data/synthetic/ground_truth_gene_topology.csv
# R must have LongOmicsTraj, edgeR and limma installed. Set RSCRIPT to override the
# R command (e.g. RSCRIPT="micromamba run -n lot Rscript").
set -euo pipefail
H5AD=$1; OUT=$2; TRUTH=${3:-}; EST=${4:-empirical,lmm}
RSCRIPT=${RSCRIPT:-Rscript}
mkdir -p "$OUT"

python3 scripts/stress_check.py --h5ad "$H5AD" --out "$OUT/stress_check"
python3 scripts/pseudobulk.py --h5ad "$H5AD" --mode sum --out "$OUT/pseudobulk_counts" --min-cells 20
$RSCRIPT scripts/normalise_pseudobulk.R "$OUT/pseudobulk_counts" "$OUT/pseudobulk_voom"
$RSCRIPT scripts/run_longomicstraj.R "$OUT/pseudobulk_voom" "$OUT/longomicstraj" 0.35 1 -Inf "$EST"
if [[ -n "$TRUTH" ]]; then
  python3 scripts/score_longomicstraj.py --res "$OUT/longomicstraj" --truth-file "$TRUTH" \
    --transitions "$OUT/pseudobulk_voom/transitions.csv" --query-estimator "${EST##*,}"
fi
