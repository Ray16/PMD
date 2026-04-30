# PMDsc Boltz-2 Cofolding Results — Best Structures

## Structures

| File | Variant | Mutations | Confidence | Source model |
|---|---|---|---|---|
| WT.pdb | Wild-type PMDsc | None | 0.9391 | model_0 of 5 |
| Y19H.pdb | Y19H | Y19H | 0.9384 | model_0 of 5 |
| R74G.pdb | R74G | R74G | 0.9400 | model_0 of 5 |
| R74H_R147K_M212Q.pdb | HKQ triple | R74H, R147K, M212Q | 0.9395 | model_0 of 5 |
| R74G_R147K_M212Q.pdb | GKQ triple | R74G, R147K, M212Q | 0.9424 | model_0 of 5 |

## Confidence Metrics

| Variant | Confidence | pTM | ipTM | ligand_ipTM | pLDDT | ipLDDT |
|---|---|---|---|---|---|---|
| WT | 0.9391 | 0.9750 | 0.9831 | 0.9831 | 0.9282 | 0.8995 |
| Y19H | 0.9384 | 0.9797 | 0.9870 | 0.9870 | 0.9263 | 0.8955 |
| R74G | 0.9400 | 0.9745 | 0.9845 | 0.9845 | 0.9289 | 0.9069 |
| R74H_R147K_M212Q | 0.9395 | 0.9776 | 0.9854 | 0.9854 | 0.9280 | 0.8990 |
| R74G_R147K_M212Q | 0.9424 | 0.9741 | 0.9853 | 0.9853 | 0.9317 | 0.9138 |

Metric definitions:
- **Confidence**: overall model confidence (0-1)
- **pTM**: predicted TM-score for individual chains
- **ipTM**: interface predicted TM-score across chains
- **ligand_ipTM**: ipTM focusing on protein-ligand interface
- **pLDDT**: predicted local distance difference test (per-residue confidence)
- **ipLDDT**: interface pLDDT at chain-chain contacts

Per-chain pTM (chain 0=Protein, 1=MVAP, 2=ATP, 3=Mg2+):

| Variant | Protein | MVAP | ATP | Mg2+ |
|---|---|---|---|---|
| WT | 0.974 | 0.970 | 0.956 | 0.000 |
| Y19H | 0.979 | 0.977 | 0.964 | 0.000 |
| R74G | 0.973 | 0.970 | 0.963 | 0.000 |
| R74H_R147K_M212Q | 0.977 | 0.968 | 0.961 | 0.000 |
| R74G_R147K_M212Q | 0.973 | 0.970 | 0.965 | 0.000 |

Note: Mg2+ pTM = 0.000 is expected (single atom, no internal structure to score).

## Cofolding Setup

### Protein
- **Enzyme**: PMDsc (phosphomevalonate decarboxylase, S. cerevisiae, UniProt P32377)
- **Sequence**: Full-length 396 aa (residue numbering matches PDB 1FI4)
- **Chain ID**: A

### Ligands

| Chain | Molecule | Representation | Description |
|---|---|---|---|
| L | Mevalonate 5-phosphate (MVAP) | SMILES: `C[C@](O)(CCOP(=O)(O)O)CC(=O)O` | Non-native substrate (R-form, MW 228.1) |
| T | ATP | CCD code: `ATP` | Co-substrate for phosphoryl transfer |
| M | Mg2+ | CCD code: `MG` | Essential catalytic cofactor |

### Pocket Constraints

**MVAP pocket** (chain L, max_distance = 6.0 A):
Constrained to 10 active site residues on chain A:
- K18, Y19, W20, K22 — substrate binding residues
- S120 — substrate anchoring
- S153 — dual substrate/ATP contact
- S155 — ATP gamma-phosphoryl stabilizer
- R158 — substrate carboxylate interaction / decarboxylation
- S208 — nucleotide positioning
- D302 — catalytic base

**ATP pocket** (chain T, max_distance = 6.0 A):
Constrained to 4 ATP binding residues on chain A:
- S121 — ATP gamma-phosphoryl stabilizer
- S153 — dual substrate/ATP contact
- S155 — ATP gamma-phosphoryl stabilizer
- S208 — nucleotide positioning

**Mg2+ contact** (chain M, max_distance = 4.0 A):
- Constrained to D302 (catalytic base) on chain A

### Boltz-2 Parameters

| Parameter | Value |
|---|---|
| Software | Boltz-2 v2.2.1 |
| MSA | MMSeqs2 server (--use_msa_server) |
| Accelerator | GPU (NVIDIA Tesla T4) |
| Recycling steps | 3 |
| Sampling steps | 200 |
| Diffusion samples | 5 (best of 5 selected by confidence) |
| Output format | PDB |

### Residue Sources

Active site residues identified by:
1. Sequence alignment of PMDsc (1FI4, S. cerevisiae) to S. epidermidis MDD (4DPT, 38.8% identity)
2. Verified by amino acid conservation, burial analysis, and spatial clustering
3. Cross-referenced against: Kang et al. 2016/2019, Chen et al. 2020 (Nat Commun), Barta et al. 2012, Michihara et al. 2008
4. See `structures/PDB/1FI4_catalytic_residues.md` for full documentation

## Visualization

```
cd /nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2/best_structures
pymol visualize_best.pml
```
