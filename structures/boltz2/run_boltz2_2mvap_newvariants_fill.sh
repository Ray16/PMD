#!/bin/bash
# Boltz-2 cofolding: fill gap at models 5-9 for 14 new variants
# Runs seeds 45-49 (5 per variant), which will be mapped to models 5-9 by merge script
# 14 variants on GPUs 6-19, each running 5 seeds sequentially

BOLTZ2_DIR=/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2
INPUT_DIR=$BOLTZ2_DIR/cofold_2mvap
OUTPUT_DIR=$BOLTZ2_DIR/output_2mvap
LOG=$OUTPUT_DIR/run_newvariants_fill.log

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

echo "$(date): Starting fill run — seeds 45-49 × ${#VARIANTS[@]} variants on GPUs 6-19" | tee $LOG

GPU=6
for VARIANT in "${VARIANTS[@]}"; do
    YAML="$INPUT_DIR/${VARIANT}.yaml"
    echo "$(date): Starting $VARIANT seeds 45-49 on GPU $GPU" | tee -a $LOG
    (
        for SEED in 45 46 47 48 49; do
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
        echo "$(date): $VARIANT fill complete (seeds 45-49)" | tee -a $LOG
    ) &
    GPU=$((GPU + 1))
done

wait
echo "$(date): All fill jobs finished." | tee -a $LOG
