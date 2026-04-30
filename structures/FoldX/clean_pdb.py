"""
Clean PDB 1FI4 for FoldX input.

Problem: 1FI4 was solved using selenomethionine (MSE) substitution.
All 10 Met residues are MSE with Se instead of S, stored as HETATM records.
FoldX cannot mutate HETATM/MSE residues directly (and M212 is a mutation target).

This script:
  1. Converts MSE → MET (HETATM→ATOM, MSE→MET, SE→SD)
  2. Removes His-tag residues (residue numbers ≤ 0)
  3. Removes water molecules (HOH)
  4. Keeps only chain A ATOM/TER/END records
"""

import os

PDB_DIR = os.path.dirname(os.path.abspath(__file__)) + "/../PDB"
infile = os.path.join(PDB_DIR, "1FI4.pdb")
outfile = os.path.join(PDB_DIR, "1FI4_clean.pdb")

with open(infile) as f:
    lines = f.readlines()

out = []
for line in lines:
    # Skip water
    if "HOH" in line:
        continue

    if line.startswith("ATOM") or line.startswith("HETATM"):
        # Get residue number (PDB columns 23-26)
        try:
            resnum = int(line[22:26].strip())
        except ValueError:
            continue

        # Skip His-tag residues (resnum <= 0)
        if resnum <= 0:
            continue

        # Convert MSE (selenomethionine) to MET (methionine)
        if "MSE" in line:
            line = "ATOM  " + line[6:]       # HETATM → ATOM
            line = line.replace("MSE", "MET")  # MSE → MET
            line = line.replace("SE  ", "SD  ")  # Se atom → S atom
            # Fix element column (cols 77-78)
            if len(line) > 76 and line[76:78].strip() == "SE":
                line = line[:76] + " S" + line[78:]

        out.append(line)
    elif line.startswith("TER"):
        out.append(line)
    elif line.startswith("END"):
        out.append(line)
        break
    elif line.startswith(("HEADER", "TITLE", "COMPND", "REMARK", "CRYST", "SCALE", "ORIG")):
        out.append(line)

with open(outfile, "w") as f:
    f.writelines(out)

# Verify
mse_atoms = sum(1 for l in out if l.startswith("ATOM") and "MSE" in l)
met_atoms = sum(1 for l in out if l.startswith("ATOM") and l[17:20] == "MET")
atom_count = sum(1 for l in out if l.startswith("ATOM"))
print(f"Written: {outfile}")
print(f"  Total ATOM lines: {atom_count}")
print(f"  MET atoms: {met_atoms}")
print(f"  Remaining MSE atoms: {mse_atoms} (should be 0)")
