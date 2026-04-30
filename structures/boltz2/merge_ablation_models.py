#!/usr/bin/env python3
"""Merge ablation seed directories into unified per-variant prediction directories.
Creates model_0 through model_9 from seed0 through seed9 outputs.
Works for both no_constraint and no_cofactor conditions."""

import shutil
import sys
from pathlib import Path

BOLTZ2_DIR = Path("/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2")

VARIANTS = [
    "PMDsc_WT_2mvap",
    "PMDsc_Y19H_2mvap",
    "PMDsc_K22M_2mvap",
    "PMDsc_K22Y_2mvap",
    "PMDsc_R74G_2mvap",
    "PMDsc_R74H_2mvap",
    "PMDsc_I145A_2mvap",
    "PMDsc_I145F_2mvap",
    "PMDsc_R147K_2mvap",
    "PMDsc_S186C_2mvap",
    "PMDsc_S208E_2mvap",
    "PMDsc_T209D_2mvap",
    "PMDsc_M212Q_2mvap",
    "PMDsc_I226V_2mvap",
    "PMDsc_V230E_2mvap",
    "PMDsc_R74G_R147K_2mvap",
    "PMDsc_R74H_R147K_M212Q_2mvap",
    "PMDsc_R74G_R147K_M212Q_2mvap",
    "PMDsc_R74G_R147K_Q140L_2mvap",
]

NUM_SEEDS = 10


def merge_condition(condition_name):
    output_dir = BOLTZ2_DIR / f"output_2mvap_{condition_name}"
    if not output_dir.exists():
        print(f"ERROR: {output_dir} does not exist")
        return False

    print(f"\n{'='*60}")
    print(f"Merging condition: {condition_name}")
    print(f"{'='*60}")

    all_ok = True
    for variant in VARIANTS:
        # Create unified prediction directory
        pred_dir = output_dir / variant / f"boltz_results_{variant}" / "predictions" / variant
        pred_dir.mkdir(parents=True, exist_ok=True)

        merged = 0
        for seed in range(NUM_SEEDS):
            seed_dir = output_dir / f"{variant}_seed{seed}" / f"boltz_results_{variant}" / "predictions" / variant

            src_pdb = seed_dir / f"{variant}_model_0.pdb"
            dst_pdb = pred_dir / f"{variant}_model_{seed}.pdb"

            if not src_pdb.exists():
                print(f"  WARNING: missing {variant} seed {seed}")
                continue

            shutil.copy2(src_pdb, dst_pdb)

            # Copy confidence JSON
            src_conf = seed_dir / f"confidence_{variant}_model_0.json"
            dst_conf = pred_dir / f"confidence_{variant}_model_{seed}.json"
            if src_conf.exists():
                shutil.copy2(src_conf, dst_conf)

            # Copy PAE/PDE
            for prefix in ["pae", "pde"]:
                src_npz = seed_dir / f"{prefix}_{variant}_model_0.npz"
                dst_npz = pred_dir / f"{prefix}_{variant}_model_{seed}.npz"
                if src_npz.exists():
                    shutil.copy2(src_npz, dst_npz)

            merged += 1

        total = len(list(pred_dir.glob(f"{variant}_model_*.pdb")))
        status = "OK" if total == NUM_SEEDS else "INCOMPLETE"
        print(f"  {variant}: {merged} merged, {total} total [{status}]")
        if total != NUM_SEEDS:
            all_ok = False

    if all_ok:
        print(f"\nAll {len(VARIANTS)} variants have {NUM_SEEDS} models. Cleaning up seed directories...")
        for variant in VARIANTS:
            removed = 0
            for seed in range(NUM_SEEDS):
                seed_dir = output_dir / f"{variant}_seed{seed}"
                if seed_dir.exists():
                    shutil.rmtree(seed_dir)
                    removed += 1
            # Also clean up log files
            for seed in range(NUM_SEEDS):
                log_file = output_dir / f"{variant}_seed{seed}.log"
                if log_file.exists():
                    log_file.unlink()
        print(f"  Cleaned up seed directories for {condition_name}")
    else:
        print(f"\nWARNING: Some variants incomplete for {condition_name}. Skipping cleanup.")

    return all_ok


if __name__ == "__main__":
    conditions = sys.argv[1:] if len(sys.argv) > 1 else ["no_constraint", "no_cofactor"]

    results = {}
    for cond in conditions:
        results[cond] = merge_condition(cond)

    print(f"\n{'='*60}")
    print("Summary:")
    for cond, ok in results.items():
        print(f"  {cond}: {'ALL OK' if ok else 'INCOMPLETE'}")
