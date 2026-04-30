# PMD Project Status

**Last updated**: 2026-04-13

## Project Goal

Engineer PMDsc (phosphomevalonate decarboxylase, S. cerevisiae) for the IPP-bypass mevalonate pathway to produce isopentenol. The key challenge is **substrate inhibition**: at high intracellular MVAP concentrations (100-200 mM), the enzyme is inhibited by its own substrate. We need to increase Ki while maintaining or improving kcat/Km.

## Variants Under Study

| Variant | Mutations | kcat/Km (mM-1s-1) | Ki (mM) | Titer (mg/L) |
|---|---|---|---|---|
| WT | None | 0.066 | 18 | 475 |
| R74G | R74G | 0.04 | 110 | 975 |
| Y19H | Y19H | 0.78 | -- | 388 |
| R74H-R147K-M212Q (HKQ) | R74H, R147K, M212Q | 0.4 | 80 | 1079 |
| R74G-R147K-M212Q (GKQ) | R74G, R147K, M212Q | 0.5 | 11 | 8 |

## Completed Work

### 1. Structure Collection (`structures/`)
- **PDB/1FI4.pdb** — WT crystal structure (2.27 A, apo)
- **PDB/1FI4_clean.pdb** — cleaned (MSE→MET, His-tag removed)
- **PDB/4DPT.pdb** — S. epidermidis MDD with substrate analog + ATP (comparative reference)
- **AlphaFold/AF-P32377-F1-model_v6.pdb** — AlphaFold prediction (pLDDT 94.75)
- **PDB/1FI4_catalytic_residues.md** — documented catalytic/active site residues with literature sources

### 2. FoldX Mutant Structures (`structures/FoldX/`)
- Repaired WT: `1FI4_clean_Repair.pdb`
- 4 mutant directories: `mut_Y19H/`, `mut_R74G/`, `mut_R74H_R147K_M212Q/`, `mut_R74G_R147K_M212Q/`
- 5 runs each, ddG predictions:
  - Y19H: +0.44 kcal/mol (neutral)
  - R74G: -1.54 kcal/mol (stabilizing)
  - R74H-R147K-M212Q: +2.90 kcal/mol (destabilizing)
  - R74G-R147K-M212Q: +0.07 kcal/mol (neutral)
- **Key finding**: ddG does NOT correlate with titer — confirms Ki, not stability, drives in vivo performance
- Scripts: `run_foldx.sh`, `clean_pdb.py`

### 3. Boltz-2 Cofolding (`structures/boltz2/`)
- Cofolded all 5 variants with MVAP + ATP + Mg2+
- Pocket constraints applied to 11 verified catalytic residues
- 5 diffusion samples per variant, best models in `best_structures/`
- All models high confidence (>0.93 ipTM)
- **Key finding**: confidence metrics don't correlate with titer — static structures can't capture substrate inhibition
- Visualization: `best_structures/visualize_best.pml`
- Full metrics: `best_structures/README.md`

### 4. Catalytic Residue Identification
11 verified residues (cross-referenced with 5 papers):
- **Catalytic core**: D302 (base), R158 (decarboxylation)
- **ATP binding**: S121, S155, S208
- **Substrate binding**: K18, Y19, W20, K22, S120, S153
- Mutation sites R74, R147, M212 are NOT catalytic — they affect Ki indirectly

## In Progress

### 5. Cosolvent MD — Substrate Bath Simulation (`md/`)
**Status: System preparation running (background job)**

**Why**: Static structures can't explain Ki differences. We need MD with ~20 free MVAP molecules (~100 mM) to observe where a second MVAP binds on the enzyme surface and how mutations affect that binding.

**Pipeline**:
1. [RUNNING] `prepare_md_systems.py` — building AMBER topology for all 5 variants
   - Uses ff14SB (protein) + GAFF2 (ligands) + TIP3P (water)
   - Known issue: ATP parameterization may timeout (needs antechamber fix)
   - Check progress: `cat md/pipeline.log`
2. [TODO] Add ~20 MVAP molecules to solvent using packmol
3. [TODO] Run 10 ns MD on GPUs 0-4 (OpenMM, ~30-50 min each)
4. [TODO] Analyze: MVAP density maps, secondary binding site identification

**Key question being addressed**: Where is the secondary MVAP binding site that causes non-competitive substrate inhibition? Why does R74G increase Ki to 110 mM while R74G-R147K-M212Q drops it to 11 mM?

## Directory Layout

```
PMD/
├── STATUS.md                  ← you are here
├── goal.md                    ← project objectives
├── paper/                     ← 5 reference papers
├── discovered_paper/          ← 3 additional papers
├── foldx5_Linux/              ← FoldX 5.1 binary
├── structures/
│   ├── PDB/                   ← experimental + cleaned structures
│   ├── AlphaFold/             ← predicted structure
│   ├── FoldX/                 ← FoldX mutant structures + ddG
│   ├── boltz2/                ← Boltz-2 cofolding results
│   │   ├── best_structures/   ← best model per variant + .pml + README
│   │   ├── output/            ← 5-sample results
│   │   └── output_1sample/    ← 1-sample backup
│   └── sequences/             ← FASTA files for all variants
└── md/                        ← MD simulations (in progress)
    ├── cofold_input/          ← YAML inputs (copied from boltz2)
    ├── cofold_output/         ← PDB inputs (copied from best_structures)
    ├── md_prep/               ← AMBER topology (being generated)
    └── pipeline.log           ← preparation log
```

## Next Steps (after MD completes)

1. Analyze MVAP density maps — identify secondary binding hot spots
2. Compare hot spots across variants — correlate with experimental Ki
3. If secondary site is identified: design mutations to disrupt it → higher Ki
4. Consider longer MD (50-100 ns) if 10 ns doesn't show clear binding events
