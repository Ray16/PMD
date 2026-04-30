#!/bin/bash
# FoldX 5.1 workflow for generating PMD mutant structures
# Starting from PDB 1FI4 (WT PMDsc, S. cerevisiae)

FOLDX=/nfs/lambda_stor_01/homes/rzhu/PMD/foldx5_Linux/foldx_20270131
PDB_DIR=/nfs/lambda_stor_01/homes/rzhu/PMD/structures/PDB
FOLDX_DIR=/nfs/lambda_stor_01/homes/rzhu/PMD/structures/FoldX

# =============================================================================
# STEP 0: Clean the PDB
# =============================================================================
# 1FI4 was solved with selenomethionine (MSE) substitution. FoldX cannot
# handle HETATM/MSE records, so we convert MSE→MET before anything else.
# Also remove the His-tag (residues ≤0) and water molecules (HOH).
#
# This is done by the Python script: clean_pdb.py (see below)
# Output: PDB_DIR/1FI4_clean.pdb

python3 clean_pdb.py

# =============================================================================
# STEP 1: RepairPDB
# =============================================================================
# Fixes bad rotamers, fills missing side-chain atoms, and optimizes the
# structure so FoldX energy calculations are reliable.
# Output: FOLDX_DIR/1FI4_clean_Repair.pdb

$FOLDX \
  --command=RepairPDB \
  --pdb=1FI4_clean.pdb \
  --pdb-dir=$PDB_DIR/ \
  --output-dir=$FOLDX_DIR/

# =============================================================================
# STEP 2: BuildModel (generate mutant structures)
# =============================================================================
# FoldX mutation syntax: WTresidue + Chain + Position + MutantResidue
#   Single mutation:   RA74G;
#   Multiple mutations on one line (applied together): RA74H,RA147K,MA212Q;
#
# --numberOfRuns=5  →  5 independent rotamer sampling runs per mutant
# --out-pdb=true    →  output PDB files (not just energies)
#
# Each mutant gets its own subdirectory and runs in parallel.

# --- Create mutation files ---
mkdir -p $FOLDX_DIR/mut_Y19H
echo "YA19H;" > $FOLDX_DIR/mut_Y19H/individual_list.txt

mkdir -p $FOLDX_DIR/mut_R74G
echo "RA74G;" > $FOLDX_DIR/mut_R74G/individual_list.txt

mkdir -p $FOLDX_DIR/mut_R74H_R147K_M212Q
echo "RA74H,RA147K,MA212Q;" > $FOLDX_DIR/mut_R74H_R147K_M212Q/individual_list.txt

mkdir -p $FOLDX_DIR/mut_R74G_R147K_M212Q
echo "RA74G,RA147K,MA212Q;" > $FOLDX_DIR/mut_R74G_R147K_M212Q/individual_list.txt

# --- Run all 4 in parallel ---
for MUT_DIR in mut_Y19H mut_R74G mut_R74H_R147K_M212Q mut_R74G_R147K_M212Q; do
  (
    cd $FOLDX_DIR/$MUT_DIR
    $FOLDX \
      --command=BuildModel \
      --pdb=1FI4_clean_Repair.pdb \
      --pdb-dir=$FOLDX_DIR/ \
      --output-dir=./ \
      --mutant-file=./individual_list.txt \
      --numberOfRuns=5 \
      --out-pdb=true \
      > foldx.log 2>&1
    echo "$MUT_DIR: done (exit code $?)"
  ) &
done
wait
echo "All BuildModel jobs finished."

# =============================================================================
# OUTPUT FILES EXPLAINED
# =============================================================================
# Each mut_* subdirectory contains:
#
# PDB files:
#   1FI4_clean_Repair_1_0.pdb ... _1_4.pdb   Mutant structures (5 runs)
#   WT_1FI4_clean_Repair_1_0.pdb ... _1_4.pdb WT re-optimized (for fair comparison)
#
# Energy reports:
#   Dif_1FI4_clean_Repair.fxout     ddG (mutant - WT) per run
#   Average_1FI4_clean_Repair.fxout Mean ddG ± SD across runs
#   Raw_1FI4_clean_Repair.fxout     Absolute energies
#   PdbList_1FI4_clean_Repair.fxout List of output PDBs
#
# ddG interpretation:
#   Negative ddG → mutation is stabilizing
#   Positive ddG → mutation is destabilizing
#   |ddG| < 0.5 kcal/mol → roughly neutral
