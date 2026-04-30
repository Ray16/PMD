#!/bin/bash
# Boltz-2 cofolding: 14 new PMDsc variants with 2x MVAP + ATP + Mg2+
# Runs 5 seeds per variant sequentially (T4 can only handle 1 sample at a time)
# Each variant gets its own GPU (GPUs 6-19)

BOLTZ2_DIR=/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2
INPUT_DIR=$BOLTZ2_DIR/cofold_2mvap
OUTPUT_DIR=$BOLTZ2_DIR/output_2mvap
LOG=$OUTPUT_DIR/run_newvariants.log

mkdir -p $OUTPUT_DIR

VARIANTS=(
    PMDsc_K22M_2mvap
    PMDsc_K22Y_2mvap
    PMDsc_R74H_2mvap
    PMDsc_I145A_2mvap
    PMDsc_I145F_2mvap
    PMDsc_R147K_2mvap
    PMDsc_S186C_2mvap
    PMDsc_S208E_2mvap
    PMDsc_T209D_2mvap
    PMDsc_M212Q_2mvap
    PMDsc_I226V_2mvap
    PMDsc_V230E_2mvap
    PMDsc_R74G_R147K_2mvap
    PMDsc_R74G_R147K_Q140L_2mvap
)

echo "$(date): Starting 2-MVAP cofolding for ${#VARIANTS[@]} new variants (5 seeds each)" | tee $LOG

GPU=6
for VARIANT in "${VARIANTS[@]}"; do
    YAML="$INPUT_DIR/${VARIANT}.yaml"
    echo "$(date): Starting $VARIANT on GPU $GPU (5 seeds)" | tee -a $LOG
    (
        for SEED in 0 1 2 3 4; do
            OUT_SUBDIR="$OUTPUT_DIR/${VARIANT}_seed${SEED}"
            rm -rf "$OUT_SUBDIR" 2>/dev/null
            echo "$(date): $VARIANT seed $SEED starting (GPU $GPU)" >> $LOG
            CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
                "$YAML" \
                --out_dir "$OUT_SUBDIR" \
                --use_msa_server \
                --accelerator gpu \
                --recycling_steps 3 \
                --sampling_steps 200 \
                --diffusion_samples 1 \
                --output_format pdb \
                --no_kernels \
                --seed $SEED \
                > "$OUT_SUBDIR.log" 2>&1
            echo "$(date): $VARIANT seed $SEED done (GPU $GPU, exit $?)" >> $LOG
        done
        echo "$(date): $VARIANT complete (all 5 seeds)" | tee -a $LOG
    ) &
    GPU=$((GPU + 1))
done

wait
echo "$(date): All new variant cofolding jobs finished." | tee -a $LOG
