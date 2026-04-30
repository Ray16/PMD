#!/usr/bin/env python3
"""Merge fill seeds (45-49) into models 5-9 for the 14 new variants.

After this, each variant should have 50 contiguous models (0-49):
  - models 0-4:   from seeds 0-4   (newvariants run)
  - models 5-9:   from seeds 45-49 (fill run)
  - models 10-49: from seeds 5-44  (newvariants_batch2 run)

Once verified, cleans up ALL seed directories (0-49).
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


def merge_fill(variant):
    pred_dir = OUTPUT_DIR / variant / f"boltz_results_{variant}" / "predictions" / variant

    if not pred_dir.exists():
        print(f"  ERROR: prediction dir not found: {pred_dir}")
        return False

    existing = sorted(pred_dir.glob(f"{variant}_model_*.pdb"))
    print(f"  Existing models: {len(existing)}")

    merged = 0
    for i, seed in enumerate(range(45, 50)):
        model_idx = 5 + i  # seeds 45-49 → models 5-9
        seed_dir = OUTPUT_DIR / f"{variant}_seed{seed}" / f"boltz_results_{variant}" / "predictions" / variant

        src_pdb = seed_dir / f"{variant}_model_0.pdb"
        dst_pdb = pred_dir / f"{variant}_model_{model_idx}.pdb"

        if not src_pdb.exists():
            print(f"  WARNING: missing seed {seed}: {src_pdb}")
            continue

        if dst_pdb.exists():
            print(f"  WARNING: model_{model_idx} already exists, overwriting")

        shutil.copy2(src_pdb, dst_pdb)

        for prefix in ["confidence", "pae", "pde"]:
            ext = "json" if prefix == "confidence" else "npz"
            src = seed_dir / f"{prefix}_{variant}_model_0.{ext}"
            dst = pred_dir / f"{prefix}_{variant}_model_{model_idx}.{ext}"
            if src.exists():
                shutil.copy2(src, dst)

        merged += 1

    total = len(list(pred_dir.glob(f"{variant}_model_*.pdb")))
    print(f"  Merged {merged} fill models. Total PDBs: {total}")

    # Verify contiguous 0-49
    missing = []
    for m in range(50):
        if not (pred_dir / f"{variant}_model_{m}.pdb").exists():
            missing.append(m)
    if missing:
        print(f"  ERROR: missing models: {missing}")
        return False

    return True


def cleanup_all_seeds(variant):
    removed = 0
    for seed in range(50):  # seeds 0-49
        seed_dir = OUTPUT_DIR / f"{variant}_seed{seed}"
        if seed_dir.exists():
            shutil.rmtree(seed_dir)
            removed += 1
        log_file = OUTPUT_DIR / f"{variant}_seed{seed}.log"
        if log_file.exists():
            log_file.unlink()
    print(f"  Cleaned up {removed} seed directories and logs for {variant}")


if __name__ == "__main__":
    all_ok = True
    for variant in VARIANTS:
        print(f"\n=== {variant} ===")
        ok = merge_fill(variant)
        if not ok:
            all_ok = False

    if all_ok:
        print(f"\nAll {len(VARIANTS)} variants have 50 contiguous models (0-49).")
        print("Cleaning up all seed directories...")
        for variant in VARIANTS:
            cleanup_all_seeds(variant)
        print(f"\nDone! All {len(VARIANTS) * 50} models verified and seed dirs removed.")
    else:
        print("\nWARNING: Some variants incomplete. Skipping cleanup.")
        print("Check errors above and re-run after fixing.")
        sys.exit(1)
