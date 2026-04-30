#!/bin/bash
# Boltz-2 cofolding: PMDsc variants with 2x MVAP + ATP + Mg2+
# Second MVAP (L2) is unconstrained — finds secondary binding site
# Runs 5 diffusion samples per variant for statistical confidence

BOLTZ2_DIR=/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2
INPUT_DIR=$BOLTZ2_DIR/cofold_2mvap
OUTPUT_DIR=$BOLTZ2_DIR/output_2mvap
LOG=$OUTPUT_DIR/run.log

mkdir -p $OUTPUT_DIR

echo "$(date): Starting 2-MVAP cofolding for all 5 variants" | tee $LOG

GPU=6  # Start from GPU 6 (0-4 running MD, 5 running WT MD)
for YAML in $INPUT_DIR/PMDsc_*_2mvap.yaml; do
    NAME=$(basename $YAML .yaml)
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$YAML" \
            --out_dir "$OUTPUT_DIR/$NAME" \
            --use_msa_server \
            --accelerator gpu \
            --recycling_steps 3 \
            --sampling_steps 200 \
            --diffusion_samples 5 \
            --output_format pdb \
            --no_kernels \
            > "$OUTPUT_DIR/${NAME}.log" 2>&1
        echo "$(date): $NAME done (exit code $?)" | tee -a $LOG
    ) &
    GPU=$((GPU + 1))
done

wait
echo "$(date): All 2-MVAP cofolding jobs finished." | tee $LOG
