# Catalytic and Active Site Residues of PMDsc (1FI4, S. cerevisiae)

Residues identified by sequence alignment to S. epidermidis MDD (38.8% identity),
verified by: (1) amino acid identity conservation, (2) burial analysis in 1FI4,
(3) spatial clustering within the active site (all within 11.2 A of centroid),
and (4) cross-referenced against five project papers (see below).

Paper-level validation:
- S155E and S208E completely abolish activity (Kang et al. 2016) — confirms essential roles
- S120, S121 provide "essential hydrogen bonds" at C-terminal end of H2 (Kang et al. 2019)
- Y19 "interacts with the phosphate moieties" of MVAP; Y19H increases kcat/KM 10-fold (Kang et al. 2019)
- K22 is "adjacent to the beta-phosphate of MVAPP" (Kang et al. 2016)
- R74, I145 are "located near the active site" but "unlikely to interact directly with substrates" (Kang et al. 2016)
- Picrophilus paper confirms conserved positions 17-21, 72, 139-141, 192-196, 283 (S.epi numbering)

Two candidate residues were excluded after verification:
- K203: surface-exposed outlier (19.0 A from centroid), not conserved (S187 in S. epidermidis)
- T209: not conserved (R193 in S. epidermidis), functional role likely shifted

## Catalytic Residues

| Yeast Res | AA  | Role                                      | S. epidermidis equivalent | Conservation | Source                                             |
|-----------|-----|-------------------------------------------|--------------------------|--------------|-----------------------------------------------------|
| D302      | Asp | Catalytic base, deprotonates 3'-OH of MVAPP | D283                     | Conserved    | Michihara et al. 2008 (PMC2279928); Chen et al. 2020 (Nat Commun) |
| R158      | Arg | Substrate carboxylate interaction, decarboxylation | R144                | Conserved    | Michihara et al. 2008 (PMC2279928)                  |

## ATP Binding Site Residues

| Yeast Res | AA  | Role                                      | S. epidermidis equivalent | Conservation | Source                                             |
|-----------|-----|-------------------------------------------|--------------------------|--------------|-----------------------------------------------------|
| S121      | Ser | Stabilizes gamma-phosphoryl group of ATP  | S107                     | Conserved    | Barta et al. 2012 (PMC4227304)                      |
| S155      | Ser | Stabilizes gamma-phosphoryl group of ATP  | S141                     | Conserved    | Barta et al. 2012 (PMC4227304)                      |
| S208      | Ser | Nucleotide positioning                    | S192                     | Conserved    | Barta et al. 2012 (PMC4227304)                      |

## Substrate (MVAP/MVAPP) Binding Site Residues

| Yeast Res | AA  | Role                                      | S. epidermidis equivalent | Conservation | Source                                             |
|-----------|-----|-------------------------------------------|--------------------------|--------------|-----------------------------------------------------|
| K18       | Lys | Substrate binding                         | K17                      | Conserved    | Chen et al. 2020 (Nat Commun, 10.1038/s41467-020-17733-0) |
| Y19       | Tyr | Substrate binding                         | Y18                      | Conserved    | Chen et al. 2020 (Nat Commun); also mutation target (Y19H) |
| W20       | Trp | Substrate binding                         | W19                      | Conserved    | Chen et al. 2020 (Nat Commun)                       |
| K22       | Lys | Substrate binding                         | K21                      | Conserved    | Chen et al. 2020 (Nat Commun)                       |
| S120      | Ser | Substrate anchoring                       | S106                     | Conserved    | Barta et al. 2012 (PMC4227304)                      |
| S153      | Ser | Dual ligand contact (substrate + ATP)     | S139                     | Conserved    | Barta et al. 2012 (PMC4227304)                      |

## Engineering Mutation Sites (Non-Catalytic)

These residues are NOT part of the catalytic machinery but affect substrate inhibition
kinetics (Ki for MVAP). They were identified by Kang et al. 2017.

| Yeast Res | AA  | Role                                         | Source                          |
|-----------|-----|----------------------------------------------|---------------------------------|
| Y19       | Tyr | Substrate binding AND substrate inhibition   | Kang et al. 2017 (ACS Chem Biol, 10.1021/acschembio.9b00322) |
| R74       | Arg | Substrate inhibition site (Ki sensitivity)   | Kang et al. 2017                |
| R147      | Arg | Substrate inhibition / catalytic region      | Kang et al. 2017                |
| M212      | Met | Substrate inhibition / active site vicinity  | Kang et al. 2017                |

## References

1. Michihara A, et al. "Identification of active site residues in mevalonate diphosphate decarboxylase." PMC2279928.
2. Chen CL, et al. "Visualizing the enzyme mechanism of mevalonate diphosphate decarboxylase." Nat Commun (2020). DOI: 10.1038/s41467-020-17733-0.
3. Barta ML, et al. "Structural basis for nucleotide binding and reaction catalysis in mevalonate diphosphate decarboxylase." Biochemistry (2012). PMC4227304.
4. Vinokur JM, et al. "A single amino acid mutation converts (R)-5-diphosphomevalonate decarboxylase into a kinase." JBC (2017). PMC5313113.
5. Kang A, et al. "Isopentenyl diphosphate (IPP)-bypass mevalonate pathways for isopentenol production." / "Substrate specificity and engineering of mevalonate 5-phosphate decarboxylase." ACS Chem Biol (2017/2019). DOI: 10.1021/acschembio.9b00322.
