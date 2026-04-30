# PMD Project Status

**Last updated**: 2026-04-30
**Repository**: https://github.com/Ray16/PMD

## Project Goal

Use molecular modeling to understand the mechanism and improve catalytic efficiency of PMDsc (phosphomevalonate decarboxylase, S. cerevisiae) for the IPP-bypass mevalonate pathway to produce isopentenol. The key challenge is **substrate inhibition**: at high intracellular MVAP concentrations (100-200 mM), the enzyme is inhibited by its own substrate. We aim to elucidate the structural basis of substrate inhibition and engineer variants with higher Ki while maintaining or improving kcat/Km.

## Variants Under Study

| Variant | Mutations | kcat/Km (mM-1s-1) | Ki (mM) | Titer (mg/L) |
|---|---|---|---|---|
| WT | None | 0.066 | 18 | 475 |
| R74G | R74G | 0.04 | 110 | 975 |
| Y19H | Y19H | 0.78 | -- | 388 |
| R74H-R147K-M212Q (HKQ) | R74H, R147K, M212Q | 0.4 | 80 | 1079 |
| R74G-R147K-M212Q (GKQ) | R74G, R147K, M212Q | 0.5 | 11 | 8 |

**Additional variants** (14 new, cofolded with 2-MVAP): I145A, I145F, I226V, K22M, K22Y, M212Q, R147K, R74G_R147K, R74G_R147K_Q140L, R74H, S186C, S208E, T209D, V230E

## Completed Work

### 1. Structure Collection (`structures/`)
- **PDB/1FI4.pdb** — WT crystal structure (2.27 A, apo)
- **PDB/1FI4_clean.pdb** — cleaned (MSE→MET, His-tag removed)
- **PDB/4DPT.pdb** — S. epidermidis MDD with substrate analog + ATP (comparative reference)
- **AlphaFold/AF-P32377-F1-model_v6.pdb** — AlphaFold prediction (pLDDT 94.75)
- **PDB/1FI4_catalytic_residues.md** — documented catalytic/active site residues with literature sources

### 2. Catalytic Residue Identification
11 verified residues (cross-referenced with 5 papers):
- **Catalytic core**: D302 (base), R158 (decarboxylation)
- **ATP binding**: S121, S155, S208
- **Substrate binding**: K18, Y19, W20, K22, S120, S153
- Mutation sites R74, R147, M212 are NOT catalytic — they affect Ki indirectly

### 3. FoldX Mutant Structures (`structures/FoldX/`)
- Repaired WT: `1FI4_clean_Repair.pdb`
- 4 mutant directories: `mut_Y19H/`, `mut_R74G/`, `mut_R74H_R147K_M212Q/`, `mut_R74G_R147K_M212Q/`
- 5 runs each, ddG predictions:
  - Y19H: +0.44 kcal/mol (neutral)
  - R74G: -1.54 kcal/mol (stabilizing)
  - R74H-R147K-M212Q: +2.90 kcal/mol (destabilizing)
  - R74G-R147K-M212Q: +0.07 kcal/mol (neutral)
- **Key finding**: ddG does NOT correlate with titer — confirms Ki, not stability, drives in vivo performance
- Scripts: `run_foldx.sh`, `clean_pdb.py`

### 4. Boltz-2 Cofolding (`structures/boltz2/`)
- Cofolded all 19 variants with 2 MVAP molecules + ATP + Mg2+ (50 seeds each)
- Pocket constraints applied to 11 verified catalytic residues
- Best models in `best_structures/`
- All models high confidence (>0.93 ipTM)
- **Key finding**: confidence metrics don't correlate with titer — static structures can't capture substrate inhibition
- Ablation studies completed (4 conditions):
  - Full constraints (baseline)
  - No constraints
  - No cofactors (no ATP/Mg2+)
  - No ATP
- Ablation analysis: `compare_four_conditions.py`, results in `four_condition_comparison.json`
- Orchestration: `run_boltz2_2mvap.sh`, `run_fill_and_merge.sh`
- Merge scripts: `merge_new_models.py`, `merge_new_variants.py`, `merge_fill_models.py`, `merge_ablation_models.py`

### 5. Cosolvent MD — Substrate Bath Simulation (`md/`)
- Ran MD for 5 variants (WT, R74G, Y19H, HKQ, GKQ) with ~20 free MVAP molecules (~100 mM)
- Uses ff14SB (protein) + GAFF2 (ligands) + TIP3P (water)
- Robust pipeline with NaN crash recovery: `run_md_robust.sh`
- MD outputs in `md_output/` (5 variants, ~20GB total)
- Analysis completed: MVAP density maps, contact frequencies, closest residue analysis
- Scripts: `md/analysis/analyze_substrate_bath.py`
- Outputs: density maps (`.dx`, `.npz`), contact frequency plots, distance heatmaps

### 6. Analysis & Figures (`analysis/`)
- **7 main figures** (fig1–fig7): L2 distance distributions, close binding frequency, consensus contacts, R74G competitive binding, R74 loop mechanism, GKQ paradox, confidence correlation
- **2 supplementary figures**: Ki vs close fraction, full contact heatmap
- **Structural panels**: enzyme overview, entrance analysis, Q212 bridge, R74 gateway, variant comparison
- **Ablation comparison**: four-condition pairwise entrance analysis
- **Reports**: `comprehensive_guide.md/pdf`, `substrate_inhibition_report.md`
- **3D viewers**: interactive HTML viewers for consensus sites, R74 loop, competitive binding
- Scripts: `generate_report_figures.py`, `generate_struct_figures.py`, `generate_struct_panels.py`, `run_ablation_analysis.py`, `md_to_pdf.py`

## Directory Layout

```
PMD/
├── STATUS.md                  ← you are here
├── goal.md                    ← project objectives
├── .gitignore                 ← excludes large outputs from git
├── paper/                     ← 5 reference papers
├── structures/
│   ├── PDB/                   ← experimental + cleaned structures
│   ├── AlphaFold/             ← predicted structure
│   ├── FoldX/                 ← FoldX mutant structures + ddG
│   ├── sequences/             ← FASTA files for all variants
│   └── boltz2/                ← Boltz-2 cofolding
│       ├── best_structures/   ← best model per variant
│       ├── cofold_2mvap/      ← 19 variant YAML configs
│       ├── cofold_2mvap_no_*/  ← ablation YAML configs
│       └── output_*/          ← prediction outputs (gitignored)
├── md/                        ← MD simulations
│   ├── cofold_input/          ← YAML inputs
│   ├── cofold_output/         ← PDB inputs
│   ├── analysis/              ← density maps, contact plots
│   ├── md_prep/               ← AMBER topology (gitignored)
│   └── md_output/             ← trajectories (gitignored)
├── analysis/                  ← figures, reports, HTML viewers
├── download_papers.py         ← paper download utility
└── fetch_papers.py            ← paper search utility
```

## Next Steps

1. Design mutations to disrupt secondary MVAP binding site → higher Ki
2. Consider longer MD (50-100 ns) if needed for clearer binding events
3. Validate proposed mutations with Boltz-2 cofolding + MD
