#!/bin/bash
# Run fill seeds (45-49) for 14 new variants, then merge into models 5-9
set -e

BOLTZ_DIR=/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2
PYTHON=/homes/rzhu/miniforge3/bin/python

echo "$(date): === Step 1: Running fill seeds (45-49) ==="
bash $BOLTZ_DIR/run_boltz2_2mvap_newvariants_fill.sh

echo ""
echo "$(date): === Step 2: Merging fill models (seeds 45-49 → models 5-9) ==="
cd $BOLTZ_DIR
$PYTHON merge_fill_models.py

echo ""
echo "$(date): === Step 3: Verifying final model counts ==="
for d in $BOLTZ_DIR/output_2mvap/PMDsc_*_2mvap/; do
    VARIANT=$(basename $d)
    PRED_DIR="$d/boltz_results_${VARIANT}/predictions/${VARIANT}"
    COUNT=$(ls "$PRED_DIR"/${VARIANT}_model_*.pdb 2>/dev/null | wc -l)
    echo "  $VARIANT: $COUNT PDB files"
done

echo ""
echo "$(date): === All fill+merge steps complete ==="
