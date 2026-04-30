#!/usr/bin/env python3
"""Generate YAML inputs for no-ATP ablation:
Protein + 2 MVAP + Mg2+ (no ATP)
L1 pocket-conditioned to active site, Mg2+ contact-constrained to D302.
"""
import os
from pathlib import Path

COFOLD_DIR = Path("/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2/cofold_2mvap")
OUT_DIR = Path("/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2/cofold_2mvap_no_atp")
OUT_DIR.mkdir(exist_ok=True)

for yaml_path in sorted(COFOLD_DIR.glob("*.yaml")):
    with open(yaml_path) as f:
        text = f.read()

    name = yaml_path.stem
    lines = text.splitlines()
    seq = None
    for line in lines:
        if line.strip().startswith("sequence:"):
            seq = line.strip().split("sequence:")[-1].strip()
            break

    if not seq:
        print(f"WARNING: could not parse sequence from {yaml_path.name}")
        continue

    no_atp_yaml = f"""version: 1
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
    id: M
    ccd: MG
constraints:
- pocket:
    binder: L
    contacts:
    - - A
      - 18
    - - A
      - 19
    - - A
      - 20
    - - A
      - 22
    - - A
      - 120
    - - A
      - 153
    - - A
      - 155
    - - A
      - 158
    - - A
      - 208
    - - A
      - 302
    max_distance: 6.0
- contact:
    token1:
    - M
    - MG
    token2:
    - A
    - 302
    max_distance: 4.0
"""

    out_path = OUT_DIR / yaml_path.name
    with open(out_path, "w") as f:
        f.write(no_atp_yaml)
    print(f"{name}: wrote no_atp YAML")

print(f"\nTotal: {len(list(OUT_DIR.glob('*.yaml')))} YAMLs")
