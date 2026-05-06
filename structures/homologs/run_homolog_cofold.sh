#!/bin/bash
# Boltz-2 cofolding: PMDsc structural homologs with 2x MVAP + ATP + Mg2+
# Runs all 10 homologs in parallel, each on a separate GPU

CONFIG_DIR=/nfs/lambda_stor_01/homes/rzhu/PMD/structures/homologs/cofold_configs
OUTPUT_DIR=/nfs/lambda_stor_01/homes/rzhu/PMD/structures/homologs/output_2mvap
LOG=$OUTPUT_DIR/run.log

mkdir -p $OUTPUT_DIR

echo "$(date): Starting homolog cofolding (10 jobs in parallel)" | tee $LOG

START_GPU=${1:-10}  # First GPU to use, default 10
GPU=$START_GPU

for YAML in $CONFIG_DIR/*.yaml; do
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
            --diffusion_samples 50 \
            --output_format pdb \
            --no_kernels \
            > "$OUTPUT_DIR/${NAME}.log" 2>&1
        echo "$(date): $NAME done (exit code $?)" | tee -a $LOG
    ) &
    GPU=$((GPU + 1))
done

echo "$(date): All 10 jobs launched on GPUs ${START_GPU}-$((GPU - 1)). Waiting..." | tee -a $LOG
wait
echo "$(date): All homolog cofolding jobs finished." | tee -a $LOG
