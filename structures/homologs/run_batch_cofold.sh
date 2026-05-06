#!/bin/bash
# Boltz-2 cofolding: 90 additional homologs in batches of 10
# Each batch uses GPUs 10-19, waits for completion, then starts next batch

CONFIG_DIR=/nfs/lambda_stor_01/homes/rzhu/PMD/structures/homologs/cofold_configs
OUTPUT_DIR=/nfs/lambda_stor_01/homes/rzhu/PMD/structures/homologs/output_2mvap
LOG=$OUTPUT_DIR/batch_run.log

mkdir -p $OUTPUT_DIR

echo "$(date): Starting batch cofolding — 87 homologs in 9 batches" | tee $LOG


# === Batch 1/9 (10 jobs) ===
echo "$(date): Starting batch 1/9" | tee -a $LOG
GPU=10

NAME=homolog_Q751D8_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q6BY07_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_A0A1D8PC43_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_I1S130_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q4WNV9_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_G9BIY1_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_3F0N_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_F4JCU3_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q0P570_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q62967_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

echo "$(date): Batch 1 launched, waiting..." | tee -a $LOG
wait
echo "$(date): Batch 1 complete" | tee -a $LOG

# === Batch 2/9 (10 jobs) ===
echo "$(date): Starting batch 2/9" | tee -a $LOG
GPU=10

NAME=homolog_Q5U403_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_P53602_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_6N10_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_F8QQQ7_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_2HKE_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q97UL5_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_6N0Z_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_6N0Y_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_3QT6_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_3QT8_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

echo "$(date): Batch 2 launched, waiting..." | tee -a $LOG
wait
echo "$(date): Batch 2 complete" | tee -a $LOG

# === Batch 3/9 (10 jobs) ===
echo "$(date): Starting batch 3/9" | tee -a $LOG
GPU=10

NAME=homolog_6E2T_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_5V2L_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_D4GXZ3_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_4RKS_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_4RKZ_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_3LTO_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q6KZB1_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_7T71_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q97BM8_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q58504_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

echo "$(date): Batch 3 launched, waiting..." | tee -a $LOG
wait
echo "$(date): Batch 3 complete" | tee -a $LOG

# === Batch 4/9 (10 jobs) ===
echo "$(date): Starting batch 4/9" | tee -a $LOG
GPU=10

NAME=homolog_A8GBA5_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_5WAT_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_B2VBV2_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_B4F0A6_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_A5UDQ5_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q65UV5_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q1C960_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q9KRP1_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_A0KQH8_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q9Y946_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

echo "$(date): Batch 4 launched, waiting..." | tee -a $LOG
wait
echo "$(date): Batch 4 complete" | tee -a $LOG

# === Batch 5/9 (10 jobs) ===
echo "$(date): Starting batch 5/9" | tee -a $LOG
GPU=10

NAME=homolog_A7FKP2_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_C0PWW2_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_B7ULN0_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_A6VQK2_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_P57899_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_A5IMZ6_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_B8F7C7_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_C3MU21_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_B0R6D7_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_C4LB24_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

echo "$(date): Batch 5 launched, waiting..." | tee -a $LOG
wait
echo "$(date): Batch 5 complete" | tee -a $LOG

# === Batch 6/9 (10 jobs) ===
echo "$(date): Starting batch 6/9" | tee -a $LOG
GPU=10

NAME=homolog_A7MV01_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_P94169_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_B5ETC9_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_A1JRX5_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_A7NI09_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q6NHU3_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q32IG8_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_B1LCQ5_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_A6URV6_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_P72663_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

echo "$(date): Batch 6 launched, waiting..." | tee -a $LOG
wait
echo "$(date): Batch 6 complete" | tee -a $LOG

# === Batch 7/9 (10 jobs) ===
echo "$(date): Starting batch 7/9" | tee -a $LOG
GPU=10

NAME=homolog_1FWK_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q0W4U0_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q87M60_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_8AOZ_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_P56838_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_B7LK02_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q324G2_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_B2J0A5_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q8U0F3_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q0T6Y7_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

echo "$(date): Batch 7 launched, waiting..." | tee -a $LOG
wait
echo "$(date): Batch 7 complete" | tee -a $LOG

# === Batch 8/9 (10 jobs) ===
echo "$(date): Starting batch 8/9" | tee -a $LOG
GPU=10

NAME=homolog_Q739T5_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_B9K9C8_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_C1AKV1_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_A4W899_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q4JU24_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q7MI80_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_A9WB97_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_A5UZX0_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_B8GCS2_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_A7Z8E0_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

echo "$(date): Batch 8 launched, waiting..." | tee -a $LOG
wait
echo "$(date): Batch 8 complete" | tee -a $LOG

# === Batch 9/9 (7 jobs) ===
echo "$(date): Starting batch 9/9" | tee -a $LOG
GPU=10

NAME=homolog_A6UV93_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q3Z453_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_A8AJ37_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_Q12XK2_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_A7MIX5_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_6ZFH_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

NAME=homolog_C3PFP1_2mvap
if [ ! -d "$OUTPUT_DIR/$NAME/boltz_results_$NAME" ]; then
    echo "$(date): Starting $NAME on GPU $GPU" | tee -a $LOG
    (
        CUDA_VISIBLE_DEVICES=$GPU conda run -n boltz-2 boltz predict \
            "$CONFIG_DIR/${NAME}.yaml" \
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
else
    echo "$(date): $NAME already done, skipping" | tee -a $LOG
fi
GPU=$((GPU + 1))

echo "$(date): Batch 9 launched, waiting..." | tee -a $LOG
wait
echo "$(date): Batch 9 complete" | tee -a $LOG

echo "$(date): All batches finished!" | tee -a $LOG

# Run analysis automatically
echo "$(date): Running entrance fraction analysis..." | tee -a $LOG
python3 /nfs/lambda_stor_01/homes/rzhu/PMD/structures/homologs/analyze_all_homologs.py 2>&1 | tee -a $LOG

echo "$(date): All done!" | tee -a $LOG
