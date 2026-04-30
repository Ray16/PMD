#!/bin/bash
# Boltz-2 ablation: 2-MVAP cofolding WITHOUT ATP
# Protein + 2 MVAP + Mg2+ (no ATP)
# L1 pocket-conditioned, Mg2+ contact-constrained to D302
# 10 samples per variant (seeds 0-9), 19 variants
# Uses all 20 GPUs (0-19), one variant per GPU initially

BOLTZ2_DIR=/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2
INPUT_DIR=$BOLTZ2_DIR/cofold_2mvap_no_atp
OUTPUT_DIR=$BOLTZ2_DIR/output_2mvap_no_atp
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

GPUS=(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19)
NUM_GPUS=${#GPUS[@]}

echo "$(date): Starting no_atp ablation — 10 seeds × 19 variants on ${NUM_GPUS} GPUs" | tee $LOG

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
echo "$(date): All no_atp jobs finished." | tee -a $LOG

# Merge results
echo "$(date): Merging seed results..." | tee -a $LOG

for VARIANT in "${VARIANTS[@]}"; do
    MERGED_DIR="$OUTPUT_DIR/${VARIANT}"
    PRED_DIR="$MERGED_DIR/boltz_results_${VARIANT}/predictions/${VARIANT}"
    mkdir -p "$PRED_DIR"
    
    count=0
    for seed in $(seq 0 9); do
        SEED_PRED="$OUTPUT_DIR/${VARIANT}_seed${seed}/boltz_results_${VARIANT}/predictions/${VARIANT}"
        [ ! -d "$SEED_PRED" ] && continue
        
        if [ -f "$SEED_PRED/${VARIANT}_model_0.pdb" ]; then
            cp "$SEED_PRED/${VARIANT}_model_0.pdb" "$PRED_DIR/${VARIANT}_model_${seed}.pdb"
            for prefix in confidence pae pde plddt; do
                if [ "$prefix" = "confidence" ]; then ext=".json"; else ext=".npz"; fi
                [ -f "$SEED_PRED/${prefix}_${VARIANT}_model_0${ext}" ] && \
                    cp "$SEED_PRED/${prefix}_${VARIANT}_model_0${ext}" "$PRED_DIR/${prefix}_${VARIANT}_model_${seed}${ext}"
            done
            count=$((count + 1))
        fi
    done
    echo "  $VARIANT: $count models merged" | tee -a $LOG
done

# Clean up seed dirs
echo "$(date): Cleaning up seed directories..." | tee -a $LOG
find "$OUTPUT_DIR" -maxdepth 1 -type d -name "PMDsc_*_seed*" -exec rm -rf {} +
find "$OUTPUT_DIR" -maxdepth 1 -type f -name "PMDsc_*_seed*.log" -delete

# Verify
echo "$(date): Verification:" | tee -a $LOG
for VARIANT in "${VARIANTS[@]}"; do
    count=$(find "$OUTPUT_DIR/${VARIANT}" -name "*_model_*.pdb" 2>/dev/null | wc -l)
    echo "  $VARIANT: $count models" | tee -a $LOG
done

echo "$(date): All done!" | tee -a $LOG
