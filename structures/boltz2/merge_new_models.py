#!/usr/bin/env python3
"""Merge batch 2 models (seeds 5-44) into existing prediction directories as model_10 through model_49."""

import shutil
import sys
from pathlib import Path

OUTPUT_DIR = Path("/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2/output_2mvap")

VARIANTS = [
    "PMDsc_WT_2mvap",
    "PMDsc_Y19H_2mvap",
    "PMDsc_R74G_2mvap",
    "PMDsc_R74G_R147K_M212Q_2mvap",
    "PMDsc_R74H_R147K_M212Q_2mvap",
]

def merge_variant(variant):
    pred_dir = OUTPUT_DIR / variant / f"boltz_results_{variant}" / "predictions" / variant
    if not pred_dir.exists():
        print(f"ERROR: prediction dir not found: {pred_dir}")
        return False

    existing = sorted(pred_dir.glob(f"{variant}_model_*.pdb"))
    print(f"{variant}: {len(existing)} existing models")

    merged = 0
    for i, seed in enumerate(range(5, 45)):
        model_idx = 10 + i
        seed_dir = OUTPUT_DIR / f"{variant}_seed{seed}" / f"boltz_results_{variant}" / "predictions" / variant

        src_pdb = seed_dir / f"{variant}_model_0.pdb"
        src_conf = seed_dir / f"confidence_{variant}_model_0.json"

        dst_pdb = pred_dir / f"{variant}_model_{model_idx}.pdb"
        dst_conf = pred_dir / f"confidence_{variant}_model_{model_idx}.json"

        if not src_pdb.exists():
            print(f"  WARNING: missing {src_pdb}")
            continue

        shutil.copy2(src_pdb, dst_pdb)
        if src_conf.exists():
            shutil.copy2(src_conf, dst_conf)

        # Also copy PAE/PDE if present
        for prefix in ["pae", "pde"]:
            src_npz = seed_dir / f"{prefix}_{variant}_model_0.npz"
            dst_npz = pred_dir / f"{prefix}_{variant}_model_{model_idx}.npz"
            if src_npz.exists():
                shutil.copy2(src_npz, dst_npz)

        merged += 1

    total = len(list(pred_dir.glob(f"{variant}_model_*.pdb")))
    print(f"  Merged {merged} new models. Total PDBs: {total}")
    return total == 50


def cleanup_seed_dirs(variant):
    removed = 0
    for seed in range(5, 45):
        seed_dir = OUTPUT_DIR / f"{variant}_seed{seed}"
        if seed_dir.exists():
            shutil.rmtree(seed_dir)
            removed += 1
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
        print("\nDone! All 250 models merged successfully.")
    else:
        print("\nWARNING: Some variants don't have 50 models. Skipping cleanup.")
        print("Re-run after fixing missing models.")
        sys.exit(1)
