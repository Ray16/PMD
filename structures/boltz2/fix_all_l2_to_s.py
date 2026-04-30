#!/usr/bin/env python3
"""Rename chain 'L2' to 'S' everywhere: PDB files, YAML inputs, analysis scripts.

Boltz2 uses 'L2' as a 2-char chain ID for the second MVAP, which breaks PDB
column alignment (PDB only allows 1-char chain IDs). This script renames L2 -> S
across all files so PDB viewers work correctly and all scripts stay consistent.
"""

import re
import sys
from pathlib import Path

BOLTZ_DIR = Path("/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2")
OUTPUT_DIR = BOLTZ_DIR / "output_2mvap"
COFOLD_DIR = BOLTZ_DIR / "cofold_2mvap"


def fix_pdb_files():
    """Fix chain L2 -> S in all PDB files under output_2mvap/."""
    pdb_files = sorted(OUTPUT_DIR.rglob("*.pdb"))
    print(f"Fixing {len(pdb_files)} PDB files...")

    files_fixed = 0
    lines_fixed = 0
    for pdb in pdb_files:
        text = pdb.read_text()
        if " L2 " not in text:
            continue

        fixed_lines = []
        n = 0
        for line in text.splitlines(keepends=True):
            if line.startswith(("ATOM", "HETATM", "TER")) and " L2 " in line:
                idx = line.index(" L2 ")
                line = line[:idx + 1] + "S" + line[idx + 3:]
                n += 1
            fixed_lines.append(line)

        pdb.write_text("".join(fixed_lines))
        files_fixed += 1
        lines_fixed += n

    print(f"  Fixed {lines_fixed} lines across {files_fixed} PDB files")


def fix_yaml_inputs():
    """Fix chain id: L2 -> id: S in all YAML input files."""
    yamls = sorted(COFOLD_DIR.glob("*_2mvap.yaml"))
    print(f"Fixing {len(yamls)} YAML input files...")

    fixed = 0
    for y in yamls:
        text = y.read_text()
        if "id: L2" in text:
            y.write_text(text.replace("id: L2", "id: S"))
            fixed += 1

    print(f"  Fixed {fixed} YAML files")


def fix_python_scripts():
    """Replace L2 chain references with S in analysis and utility scripts."""
    scripts = [
        OUTPUT_DIR / "comprehensive_l2_analysis.py",
        OUTPUT_DIR / "comprehensive_l2_analysis_all.py",
        BOLTZ_DIR / "fix_l2_chain.py",
        BOLTZ_DIR / "load_boltz2.py",
    ]

    # Also pick up any structural_insights*_all.py scripts
    for p in OUTPUT_DIR.glob("structural_insights*.py"):
        if p not in scripts:
            scripts.append(p)

    # Ordered replacements — longer strings first to avoid partial matches
    replacements = [
        # Dict keys and variable names
        ("'A_L2_iptm_from_L2'", "'A_S_iptm_from_S'"),
        ("'A_L2_iptm_from_A'", "'A_S_iptm_from_A'"),
        ("'mean_A_L2_iptm'", "'mean_A_S_iptm'"),
        ("A_L2_iptm_from_L2", "A_S_iptm_from_S"),
        ("A_L2_iptm_from_A", "A_S_iptm_from_A"),
        ("mean_A_L2_iptm", "mean_A_S_iptm"),
        ("a_l2_iptms", "a_s_iptms"),
        ("l2_coords", "s_coords"),
        ("l2_heavy", "s_heavy"),
        ("l2_com", "s_com"),
        ("l2_dist", "s_dist"),
        ("l2_dists", "s_dists"),
        ("l2_closest", "s_closest"),
        ("l2_plddt", "s_plddt"),
        ("L2_ptm", "S_ptm"),
        # Chain references in quotes
        ("'L2'", "'S'"),
        ('"L2"', '"S"'),
        # Display strings
        ("A-L2", "A-S"),
        ("A→L2", "A→S"),
        ("L2 ", "S "),
        ("(L2)", "(S)"),
        ("L2-L1", "S-L"),
        ("L2-protein", "S-protein"),
        (" L2,", " S,"),
        (" L2.", " S."),
        (" L2\\n", " S\\n"),
        (" L2)", " S)"),
    ]

    print(f"Fixing {len(scripts)} Python scripts...")
    for script in scripts:
        if not script.exists():
            print(f"  Skipping (not found): {script.name}")
            continue

        text = script.read_text()
        original = text
        for old, new in replacements:
            text = text.replace(old, new)

        if text != original:
            script.write_text(text)
            print(f"  Fixed: {script.name}")
        else:
            print(f"  No changes: {script.name}")


def verify():
    """Verify no L2 chain references remain in PDB files."""
    print("\nVerification:")

    # Check a sample PDB
    sample_pdbs = list((OUTPUT_DIR / "PMDsc_WT_2mvap").rglob("*.pdb"))[:1]
    for pdb in sample_pdbs:
        with open(pdb) as f:
            for line in f:
                if line.startswith(("HETATM", "TER")) and " L2 " in line:
                    print(f"  WARNING: L2 still found in {pdb.name}")
                    return False

    # Check chain S exists
    for pdb in sample_pdbs:
        found_s = False
        with open(pdb) as f:
            for line in f:
                if line.startswith("HETATM") and line[21] == "S":
                    found_s = True
                    break
        if found_s:
            print(f"  Chain S found in {pdb.name}")
        else:
            print(f"  WARNING: Chain S not found in {pdb.name}")

    # Verify PDB column alignment
    for pdb in sample_pdbs:
        with open(pdb) as f:
            for line in f:
                if line.startswith("HETATM") and line[21] == "S":
                    chain = line[21]
                    resseq = line[22:26].strip()
                    x = line[30:38].strip()
                    print(f"  Sample S atom: chain={chain}, resSeq={resseq}, x={x}")
                    break

    print("  Verification complete")
    return True


if __name__ == "__main__":
    print("=== Renaming chain L2 -> S across all files ===\n")

    fix_pdb_files()
    print()
    fix_yaml_inputs()
    print()
    fix_python_scripts()

    verify()

    print("\nDone!")
