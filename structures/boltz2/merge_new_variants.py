#!/usr/bin/env python3
"""Merge all models (seeds 0-44) for 14 new variants into unified prediction directories.

Seeds 0-4: from run_boltz2_2mvap_newvariants.sh (model_0 through model_4)
Seeds 5-44: from run_boltz2_2mvap_newvariants_batch2.sh (model_10 through model_49)
Total: 50 models per variant, matching the original 5 variants.
"""

import shutil
import sys
from pathlib import Path

OUTPUT_DIR = Path("/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2/output_2mvap")

VARIANTS = [
    "PMDsc_K22M_2mvap",
    "PMDsc_K22Y_2mvap",
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
    "PMDsc_R74G_R147K_Q140L_2mvap",
]


def merge_variant(variant):
    # The seed0 run creates the main prediction directory
    seed0_dir = OUTPUT_DIR / f"{variant}_seed0" / f"boltz_results_{variant}" / "predictions" / variant

    # Target: unified directory under output_2mvap/<variant>/
    pred_dir = OUTPUT_DIR / variant / f"boltz_results_{variant}" / "predictions" / variant
    pred_dir.mkdir(parents=True, exist_ok=True)

    merged = 0

    # Seeds 0-4 -> model_0 through model_4
    for seed in range(5):
        model_idx = seed
        seed_dir = OUTPUT_DIR / f"{variant}_seed{seed}" / f"boltz_results_{variant}" / "predictions" / variant

        src_pdb = seed_dir / f"{variant}_model_0.pdb"
        dst_pdb = pred_dir / f"{variant}_model_{model_idx}.pdb"

        if not src_pdb.exists():
            print(f"  WARNING: missing seed {seed}: {src_pdb}")
            continue

        shutil.copy2(src_pdb, dst_pdb)

        for prefix in ["confidence", "pae", "pde"]:
            ext = "json" if prefix == "confidence" else "npz"
            src = seed_dir / f"{prefix}_{variant}_model_0.{ext}"
            dst = pred_dir / f"{prefix}_{variant}_model_{model_idx}.{ext}"
            if src.exists():
                shutil.copy2(src, dst)

        merged += 1

    # Seeds 5-44 -> model_10 through model_49 (matching original merge scheme)
    for i, seed in enumerate(range(5, 45)):
        model_idx = 10 + i
        seed_dir = OUTPUT_DIR / f"{variant}_seed{seed}" / f"boltz_results_{variant}" / "predictions" / variant

        src_pdb = seed_dir / f"{variant}_model_0.pdb"
        dst_pdb = pred_dir / f"{variant}_model_{model_idx}.pdb"

        if not src_pdb.exists():
            print(f"  WARNING: missing seed {seed}: {src_pdb}")
            continue

        shutil.copy2(src_pdb, dst_pdb)

        for prefix in ["confidence", "pae", "pde"]:
            ext = "json" if prefix == "confidence" else "npz"
            src = seed_dir / f"{prefix}_{variant}_model_0.{ext}"
            dst = pred_dir / f"{prefix}_{variant}_model_{model_idx}.{ext}"
            if src.exists():
                shutil.copy2(src, dst)

        merged += 1

    total = len(list(pred_dir.glob(f"{variant}_model_*.pdb")))
    print(f"  Merged {merged} models. Total PDBs: {total}")
    return total == 50


def cleanup_seed_dirs(variant):
    removed = 0
    for seed in range(45):
        seed_dir = OUTPUT_DIR / f"{variant}_seed{seed}"
        if seed_dir.exists():
            shutil.rmtree(seed_dir)
            removed += 1
        # Also clean up log files
        log_file = OUTPUT_DIR / f"{variant}_seed{seed}.log"
        if log_file.exists():
            log_file.unlink()
    print(f"  Cleaned up {removed} seed directories for {variant}")


if __name__ == "__main__":
    all_ok = True
    for variant in VARIANTS:
        print(f"\n=== {variant} ===")
        ok = merge_variant(variant)
        if not ok:
            all_ok = False

    if all_ok:
        print("\nAll variants have 50 models. Cleaning up seed directories...")
        for variant in VARIANTS:
            cleanup_seed_dirs(variant)
        print(f"\nDone! All {len(VARIANTS) * 50} models merged successfully.")
    else:
        print("\nWARNING: Some variants don't have 50 models. Skipping cleanup.")
        print("Re-run after fixing missing models.")
        sys.exit(1)
