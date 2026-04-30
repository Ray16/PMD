#!/bin/bash
# Wait for fill+merge pipeline to complete, then rename chain L2 -> S everywhere
set -e

BOLTZ_DIR=/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2
PYTHON=/homes/rzhu/miniforge3/bin/python

# Wait for run_fill_and_merge.sh to finish by monitoring its main PID
echo "$(date): Waiting for fill+merge pipeline to finish..."
while pgrep -f "run_fill_and_merge.sh" > /dev/null 2>&1; do
    sleep 30
done
echo "$(date): Fill+merge pipeline finished."

# Verify merge succeeded
echo ""
echo "$(date): === Verifying merge results ==="
ALL_OK=true
for d in $BOLTZ_DIR/output_2mvap/PMDsc_*_2mvap/; do
    VARIANT=$(basename $d)
    PRED_DIR="$d/boltz_results_${VARIANT}/predictions/${VARIANT}"
    COUNT=$(ls "$PRED_DIR"/${VARIANT}_model_*.pdb 2>/dev/null | wc -l)
    echo "  $VARIANT: $COUNT PDB files"
    if [ "$COUNT" -ne 50 ]; then
        ALL_OK=false
    fi
done

if [ "$ALL_OK" = false ]; then
    echo ""
    echo "ERROR: Not all variants have 50 models. Aborting chain rename."
    exit 1
fi

echo ""
echo "$(date): === All variants verified (50 models each). Running chain rename L2 -> S ==="
cd $BOLTZ_DIR
$PYTHON fix_all_l2_to_s.py

echo ""
echo "$(date): === All done! ==="
