#!/usr/bin/env python3
"""Generate YAML inputs for two ablation conditions:
1. no_constraint: protein + 2 MVAP + ATP + Mg2+, but NO pocket/contact constraints
2. no_cofactor:   protein + 2 MVAP only (no ATP, no Mg2+, no constraints)
"""

import os
from pathlib import Path

COFOLD_DIR = Path("/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2/cofold_2mvap")
OUT_NO_CONSTRAINT = Path("/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2/cofold_2mvap_no_constraint")
OUT_NO_COFACTOR = Path("/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2/cofold_2mvap_no_cofactor")

OUT_NO_CONSTRAINT.mkdir(exist_ok=True)
OUT_NO_COFACTOR.mkdir(exist_ok=True)

for yaml_path in sorted(COFOLD_DIR.glob("*.yaml")):
    with open(yaml_path) as f:
        text = f.read()

    name = yaml_path.stem  # e.g. PMDsc_WT_2mvap

    # Extract protein sequence (line after "sequence:")
    lines = text.splitlines()
    seq = None
    for i, line in enumerate(lines):
        if line.strip().startswith("sequence:"):
            seq = line.strip().split("sequence:")[-1].strip()
            break

    if not seq:
        print(f"WARNING: could not parse sequence from {yaml_path.name}")
        continue

    # --- Condition 1: no_constraint (keep all molecules, remove constraints) ---
    no_constraint_yaml = f"""version: 1
sequences:
- protein:
    id: A
    sequence: {seq}
- ligand:
    id: L
    smiles: C[C@](O)(CCOP(=O)(O)O)CC(=O)O
- ligand:
    id: S
    smiles: C[C@](O)(CCOP(=O)(O)O)CC(=O)O
- ligand:
    id: T
    ccd: ATP
- ligand:
    id: M
    ccd: MG
"""

    out_path = OUT_NO_CONSTRAINT / yaml_path.name
    with open(out_path, "w") as f:
        f.write(no_constraint_yaml)

    # --- Condition 2: no_cofactor (protein + 2 MVAP only) ---
    no_cofactor_yaml = f"""version: 1
sequences:
- protein:
    id: A
    sequence: {seq}
- ligand:
    id: L
    smiles: C[C@](O)(CCOP(=O)(O)O)CC(=O)O
- ligand:
    id: S
    smiles: C[C@](O)(CCOP(=O)(O)O)CC(=O)O
"""

    out_path = OUT_NO_COFACTOR / yaml_path.name
    with open(out_path, "w") as f:
        f.write(no_cofactor_yaml)

    print(f"{name}: wrote no_constraint + no_cofactor YAMLs")

print(f"\nTotal: {len(list(OUT_NO_CONSTRAINT.glob('*.yaml')))} no_constraint, "
      f"{len(list(OUT_NO_COFACTOR.glob('*.yaml')))} no_cofactor")
