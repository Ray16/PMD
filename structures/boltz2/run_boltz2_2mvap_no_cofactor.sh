#!/bin/bash
# Boltz-2 ablation: 2-MVAP cofolding WITHOUT cofactors
# Protein + 2 MVAP only (no ATP, no Mg2+, no constraints)
# 10 samples per variant (seeds 0-9), 19 variants
# Uses GPUs 14, 15, 16, 17, 18, 19

BOLTZ2_DIR=/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2
INPUT_DIR=$BOLTZ2_DIR/cofold_2mvap_no_cofactor
OUTPUT_DIR=$BOLTZ2_DIR/output_2mvap_no_cofactor
LOG=$OUTPUT_DIR/run.log

mkdir -p $OUTPUT_DIR

VARIANTS=(
    PMDsc_WT_2mvap
    PMDsc_Y19H_2mvap
    PMDsc_K22M_2mvap
    PMDsc_K22Y_2mvap
    PMDsc_R74G_2mvap
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
    PMDsc_R74H_R147K_M212Q_2mvap
    PMDsc_R74G_R147K_M212Q_2mvap
    PMDsc_R74G_R147K_Q140L_2mvap
)

GPUS=(14 15 16 17 18 19)
NUM_GPUS=${#GPUS[@]}

echo "$(date): Starting no_cofactor ablation — 10 seeds × 19 variants on ${NUM_GPUS} GPUs" | tee $LOG

# Distribute variants round-robin across GPUs
for i in "${!VARIANTS[@]}"; do
    VARIANT="${VARIANTS[$i]}"
    GPU_IDX=$((i % NUM_GPUS))
    GPU="${GPUS[$GPU_IDX]}"
    YAML="$INPUT_DIR/${VARIANT}.yaml"

    (
        echo "$(date): $VARIANT on GPU $GPU" | tee -a $LOG
        for SEED in $(seq 0 9); do
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
            echo "$(date): $VARIANT seed $SEED done (GPU $GPU)" >> $LOG
        done
        echo "$(date): $VARIANT complete" | tee -a $LOG
    ) &
done

wait
echo "$(date): All no_cofactor jobs finished." | tee -a $LOG
