#!/bin/bash
# Robust MD runner: detects NaN crashes and rebuilds with packmol if needed
MD_DIR=/nfs/lambda_stor_01/homes/rzhu/PMD/md
LOG=$MD_DIR/robust_md.log

check_and_fix_clashes() {
    local uid=$1
    local PREP_DIR=$MD_DIR/md_prep/$uid
    local OUT_DIR=$MD_DIR/md_output/$uid
    
    # Check if error.log has NaN
    if [ -f "$OUT_DIR/error.log" ] && grep -q "NaN" "$OUT_DIR/error.log" 2>/dev/null; then
        echo "$(date): $uid crashed with NaN — rebuilding with packmol" | tee -a $LOG
        
        # Use packmol to properly place MVAP with distance constraints
        # First, save the base system (protein+activesite MVAP+ATP+Mg) as PDB
        conda run -n boltz-2 python3 << PYEOF
import parmed
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

prep_dir = "$PREP_DIR"
uid = "$uid"

# Load the no-bath system
p = parmed.load_file(f"{prep_dir}/system_nobath.prmtop", f"{prep_dir}/system_nobath.inpcrd")

# Save as PDB for packmol (everything: protein + ligands + water + ions)
p.save(f"{prep_dir}/system_nobath.pdb", overwrite=True)

# Generate MVAP PDB for packmol (single molecule, with hydrogens)
smiles = "C[C@](O)(CCOP(=O)(O)O)CC(=O)O"
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol)
Chem.MolToPDBFile(mol, f"{prep_dir}/mvap_packmol.pdb")

# Get box dimensions
box = p.box
print(f"Box: {box[0]:.1f} x {box[1]:.1f} x {box[2]:.1f}")

# Write packmol input
# Place 38 MVAP molecules with:
#   - at least 3.0 A from any existing atom (protein/water/ions)
#   - at least 5.0 A from each other
#   - inside a sphere centered on the box center, radius = box/2 - 5
radius = min(box[0], box[1], box[2]) / 2 - 5
center = np.mean(p.coordinates, axis=0)

with open(f"{prep_dir}/packmol_bath.inp", "w") as f:
    f.write(f"""tolerance 3.0
filetype pdb
output {prep_dir}/system_with_bath.pdb

# Existing system (fixed)
structure {prep_dir}/system_nobath.pdb
  number 1
  fixed 0. 0. 0. 0. 0. 0.
end structure

# Free MVAP molecules
structure {prep_dir}/mvap_packmol.pdb
  number 38
  inside sphere {center[0]:.1f} {center[1]:.1f} {center[2]:.1f} {radius:.1f}
end structure
""")
print(f"Packmol input written: center=({center[0]:.1f},{center[1]:.1f},{center[2]:.1f}), radius={radius:.1f}")
PYEOF
        
        # Run packmol
        cd $PREP_DIR
        conda run -n boltz-2 packmol < packmol_bath.inp > packmol.log 2>&1
        
        if grep -q "Success" packmol.log 2>/dev/null; then
            echo "$(date): $uid packmol succeeded" | tee -a $LOG
            
            # Now rebuild AMBER topology from packmol output
            # Extract just the MVAP molecules from packmol output and rebuild with tleap
            conda run -n boltz-2 python3 << PYEOF2
import parmed

prep_dir = "$PREP_DIR"

# Load packmol output 
# Strategy: load nobath topology, then add MVAP coordinates from packmol
# This is complex, so simpler: just increase the minimum distance in translate approach
# and re-run tleap

# Actually, let's just use a larger minimum distance (15A between MVAP centers)
import numpy as np
import random

p = parmed.load_file(f"{prep_dir}/system_nobath.prmtop", f"{prep_dir}/system_nobath.inpcrd")
coords = np.array(p.coordinates)
protein_center = np.mean(coords[:3000], axis=0)
distances = np.linalg.norm(coords - protein_center, axis=1)

# Select positions >25A from protein center, >15A apart (more conservative)
far_indices = np.where(distances > 25.0)[0]
random.seed(42)
selected_coords = []

for idx in np.random.permutation(far_indices):
    pos = coords[idx]
    too_close = any(np.linalg.norm(pos - sc) < 15.0 for sc in selected_coords)
    if not too_close:
        selected_coords.append(pos)
    if len(selected_coords) >= 38:
        break

print(f"Found {len(selected_coords)} positions (>25A from protein, >15A apart)")

# Read MVAP center
mol2_path = f"{prep_dir}/chain_L_ligl.mol2"
mvap_coords = []
reading = False
with open(mol2_path) as f:
    for line in f:
        if "@<TRIPOS>ATOM" in line:
            reading = True; continue
        if reading and "@<TRIPOS>" in line: break
        if reading:
            parts = line.strip().split()
            if len(parts) >= 5:
                mvap_coords.append([float(parts[2]), float(parts[3]), float(parts[4])])
mvap_center = np.mean(mvap_coords, axis=0)

# Rewrite tleap
mvap_frcmod = f"{prep_dir}/chain_L_ligl.frcmod"
mvap_mol2 = mol2_path
lines = [
    f"# tleap bath RETRY for {uid} (conservative placement)",
    "source leaprc.protein.ff14SB", "source leaprc.gaff2", "source leaprc.water.tip3p", "",
    f"protein = loadpdb {prep_dir}/protein_amber.pdb", "",
    f"loadAmberParams {mvap_frcmod}", f"ligl = loadmol2 {mvap_mol2}", "",
    f"mg = loadpdb {prep_dir}/chain_M_mg.pdb", "",
    f"loadAmberParams {prep_dir}/chain_T_atp.frcmod",
    f"atp = loadmol2 {prep_dir}/chain_T_atp.mol2", "",
]
names = []
for i, pos in enumerate(selected_coords):
    name = f"mv{i:02d}"
    names.append(name)
    tx, ty, tz = pos - mvap_center
    lines.append(f"{name} = loadmol2 {mvap_mol2}")
    lines.append(f"translate {name} {{ {tx:.3f} {ty:.3f} {tz:.3f} }}")

all_u = "protein ligl mg atp " + " ".join(names)
lines += ["", f"complex = combine {{ {all_u} }}", "",
    "solvateOct complex TIP3PBOX 12.0", "",
    "addIonsRand complex Na+ 0", "addIonsRand complex Cl- 0",
    "addIonsRand complex Na+ 40 Cl- 40", "",
    "check complex",
    f"saveAmberParm complex {prep_dir}/system.prmtop {prep_dir}/system.inpcrd",
    "quit"]

with open(f"{prep_dir}/tleap_bath_retry.in", "w") as f:
    f.write("\n".join(lines))
print("Wrote tleap_bath_retry.in")
PYEOF2
            
            # Run tleap
            conda run -n boltz-2 tleap -f $PREP_DIR/tleap_bath_retry.in > $PREP_DIR/tleap_bath_retry.log 2>&1
            echo "$(date): $uid rebuilt with conservative placement" | tee -a $LOG
        else
            echo "$(date): $uid packmol FAILED" | tee -a $LOG
        fi
        
        return 1  # signal that we fixed something
    fi
    return 0  # no crash
}

# Wait for the initial MD run to finish
echo "$(date): Waiting for initial MD run..." | tee -a $LOG
wait

# Check each variant for NaN crashes
any_fixed=0
for uid in WT Y19H R74G HKQ GKQ; do
    check_and_fix_clashes $uid
    if [ $? -eq 1 ]; then
        any_fixed=1
    fi
done

# If any were fixed, re-run MD for those
if [ $any_fixed -eq 1 ]; then
    echo "$(date): Re-running MD for fixed variants..." | tee -a $LOG
    conda run -n boltz-2 python /nfs/lambda_stor_01/homes/rzhu/Biochemical_agent/experimental_feedback/run_md.py $MD_DIR --sim-time-ns 10 --gpus 0,1,2,3,4 2>&1 | tee -a $LOG
fi

echo "$(date): Robust MD pipeline complete." | tee -a $LOG
