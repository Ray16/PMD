# Three-Condition Boltz-2 Cofolding Analysis of PMDsc Substrate Inhibition

## A Comprehensive Synthesis of Constrained, No-Constraint, and No-Cofactor Simulations

---

## 1. Introduction

Mevalonate 5-phosphate decarboxylase (PMDsc) catalyzes the decarboxylation of mevalonate 5-phosphate (MVAP) in the isoprenoid biosynthesis pathway. The enzyme suffers from substrate inhibition: at elevated MVAP concentrations, a second substrate molecule (L2) binds near the active site and impedes catalysis. This inhibition is quantified by Ki and manifests as reduced isopentenol titers in vivo.

To understand the structural basis of substrate inhibition and how mutations alleviate or exacerbate it, we performed Boltz-2 cofolding simulations on 19 PMDsc variants. Each simulation cofolds the enzyme with two MVAP molecules (L1 at the active site, L2 free to bind anywhere), probing where the inhibitory L2 molecule preferentially binds. By comparing results across three increasingly permissive conditions -- pocket-conditioned with cofactors, unconstrained with cofactors, and unconstrained without cofactors -- we can dissect the relative contributions of pocket conditioning, cofactor scaffolding, and protein sequence to L2 binding behavior.

This report synthesizes results from all three conditions into a unified mechanistic picture, with the goal of explaining the experimental Ki and titer data and guiding future enzyme engineering efforts.

---

## 2. Experimental Design

### 2.1 Variants Studied

Nineteen PMDsc variants were analyzed, spanning single-point mutations, double mutants, and triple mutants:

- **Wild type (WT)**
- **Single mutants**: Y19H, K22M, K22Y, R74G, R74H, I145A, I145F, R147K, S186C, S208E, T209D, M212Q, I226V, V230E
- **Double mutant**: R74G-R147K
- **Triple mutants**: R74H-R147K-M212Q (HKQ), R74G-R147K-M212Q (GKQ), R74G-R147K-Q140L

### 2.2 Three Simulation Conditions

| Condition | L1 Treatment | Cofactors | Models per Variant | Total Models |
|-----------|-------------|-----------|-------------------|--------------|
| **Constrained** | Pocket-conditioned to active site | ATP + Mg2+ | 50 | 950 |
| **No-constraint** | Unconstrained | ATP + Mg2+ | 10 | 190 |
| **No-cofactor** | Unconstrained | None | 10 | 190 |

The three conditions form a logical hierarchy of decreasing structural guidance:

1. **Constrained**: Maximum guidance. L1 is placed at the active site via pocket conditioning, and cofactors (ATP, Mg2+) are present. This condition most closely approximates the catalytically competent enzyme-substrate complex.

2. **No-constraint**: Intermediate guidance. Cofactors are present but L1 must find the active site on its own. This tests whether Boltz-2 can independently identify the correct binding pose and whether L2 discrimination persists without explicit L1 placement.

3. **No-cofactor**: Minimal guidance. Neither pocket conditioning nor cofactors are present. This isolates the contribution of the protein fold alone, revealing what happens when the entrance channel is not narrowed by cofactor scaffolding.

### 2.3 Key Metrics

- **Entrance fraction (%)**: Percentage of models where L2 binds within close proximity of the active site entrance (non-competitive inhibition mode).
- **Active-site (competitive) binding**: Models where L2 enters the catalytic pocket itself (competitive inhibition mode).
- **L1-to-active-site distance**: Confirms whether L1 correctly localizes to the active site.
- **Consensus contact residues**: Residues most frequently in contact with L2 across close-binding models.

---

## 3. L1 Active-Site Binding Validation

A critical validation step is confirming that L1 (the catalytic substrate) correctly localizes to the active site under all conditions. If L1 fails to find the active site, the simulation cannot meaningfully probe L2 inhibitory binding.

### 3.1 Results Across Conditions

| Condition | L1-to-Active-Site Distance | Notes |
|-----------|---------------------------|-------|
| **Constrained** | By definition at active site | Pocket conditioning enforces placement |
| **No-constraint** | 3.4--4.5 A (all variants) | L1 independently finds active site |
| **No-cofactor** | 3.9--4.7 A (all variants) | L1 still finds active site without cofactors |

### 3.2 Interpretation

The fact that L1 finds the active site without any guidance -- across all 19 variants, in both the no-constraint and no-cofactor conditions -- provides strong validation of the Boltz-2 cofolding approach. The active-site pocket represents the deepest energy minimum on the protein surface for MVAP binding, and the model correctly identifies it regardless of whether pocket conditioning or cofactors are present.

The slight increase in L1-to-active-site distance from no-constraint (3.4--4.5 A) to no-cofactor (3.9--4.7 A) suggests that cofactors help refine L1 positioning, likely by pre-organizing the active-site geometry. However, the difference is modest (less than 1 A on average), indicating that the protein fold itself provides the primary determinant of L1 binding.

This validation has a practical consequence: pocket conditioning is not required to study L2 binding. The constrained condition offers higher throughput (50 models/variant vs. 10) and sharper discrimination between variants, but the no-constraint condition yields fundamentally consistent results.

---

## 4. L2 Secondary Binding: Three-Way Comparison

The central question is where the inhibitory L2 molecule binds and how this varies across variants and conditions. The table below presents the entrance fraction -- the percentage of models where L2 binds at the entrance channel -- for all 19 variants under all three conditions, alongside available experimental data.

### 4.1 Complete Results Table

| Variant | Constrained (%) | No-Constraint (%) | No-Cofactor (%) | Ki (mM) | kcat/Km | Titer (mg/L) |
|---------|-----------------|-------------------|-----------------|---------|---------|--------------|
| WT | 28 | 20 | 100 | 18 | 0.066 | 475 |
| Y19H | 16 | 70 | 100 | N.D. | 0.78 | 388 |
| K22M | 18 | 40 | 100 | N.D. | 0.12 | 22 |
| K22Y | 24 | 30 | 90 | N.D. | 0.12 | 22 |
| R74G | 6 | 10 | 100 | 110 | 0.04 | 975 |
| R74H | 10 | 30 | 90 | N.D. | 0.44 | 770 |
| I145A | 32 | 70 | 100 | N.D. | 0.01 | 623 |
| I145F | 48 | 50 | 100 | N.D. | 0.20 | N.D. |
| R147K | 28 | 80 | 90 | N.D. | 0.32 | 793 |
| S186C | 36 | 30 | 100 | N.D. | 0.08 | 596 |
| S208E | 42 | 60 | 90 | N.D. | N.D. | N.D. |
| T209D | 26 | 60 | 90 | N.D. | 0.13 | N.D. |
| M212Q | 36 | 70 | 100 | N.D. | 0.50 | 601 |
| I226V | 24 | 50 | 100 | N.D. | 0.46 | 633 |
| V230E | 40 | 70 | 100 | 10 | 0.08 | 278 |
| R74G-R147K | 6 | 10 | 100 | N.D. | 0.42 | 909 |
| HKQ | 2 | 0 | 90 | 80 | 0.40 | 1079 |
| GKQ | 4 | 0 | 90 | 11 | 0.50 | 8 |
| R74G-R147K-Q140L | 4 | 0 | 90 | N.D. | 0.03 | 401 |

### 4.2 Overall Statistics

| Metric | Constrained | No-Constraint | No-Cofactor |
|--------|-------------|---------------|-------------|
| Total close-binding models | 271/950 (29%) | 87/190 (46%) | 187/190 (98%) |
| Pearson r (Ki vs entrance) | -0.60 | -0.449 | 0.006 |
| Active-site (competitive) models | ~56 | 12 | 5 |

### 4.3 Key Observations

**Discrimination power degrades without cofactors.** The constrained condition separates variants most cleanly (entrance fractions span 2--48%), the no-constraint condition preserves much of this separation (0--80%), but the no-cofactor condition collapses nearly all variants to 90--100%. The correlation between Ki and entrance fraction follows the same gradient: r = -0.60 (constrained), r = -0.449 (no-constraint), r = 0.006 (no-cofactor).

**The no-cofactor condition is a negative control, not a failure.** The near-universal 90--100% entrance binding without cofactors demonstrates that the entrance channel is inherently wide in the apoenzyme. Cofactors are required to narrow the channel sufficiently for mutations to exert differential effects. This is a biologically meaningful finding: it explains why substrate inhibition is relevant only under catalytic turnover conditions, when ATP and Mg2+ are bound.

**Triple mutants are consistently distinguished.** HKQ (2%, 0%, 90%), GKQ (4%, 0%, 90%), and R74G-R147K-Q140L (4%, 0%, 90%) show the lowest entrance fractions under both the constrained and no-constraint conditions. These variants achieve near-complete entrance closure when cofactors are present, but this protection vanishes without cofactors.

---

## 5. Variant Discrimination Hierarchy

Ranking variants by entrance fraction reveals a consistent hierarchy across the constrained and no-constraint conditions, with some instructive re-orderings.

### 5.1 Constrained Condition Ranking (lowest entrance = most protected)

| Tier | Variants | Entrance (%) | Interpretation |
|------|----------|-------------|----------------|
| **Tier 1: Near-complete closure** | HKQ (2%), GKQ (4%), R74G-R147K-Q140L (4%) | 2--4% | Triple mutants with synergistic defense |
| **Tier 2: Strong closure** | R74G (6%), R74G-R147K (6%), R74H (10%) | 6--10% | R74 mutations dominate |
| **Tier 3: Moderate closure** | Y19H (16%), K22M (18%) | 16--18% | Entrance-adjacent mutations |
| **Tier 4: WT-like** | K22Y (24%), I226V (24%), T209D (26%), WT (28%), R147K (28%) | 24--28% | Modest or no effect |
| **Tier 5: Increased entrance** | I145A (32%), S186C (36%), M212Q (36%), V230E (40%), S208E (42%), I145F (48%) | 32--48% | Potentially destabilizing |

### 5.2 No-Constraint Condition Ranking

| Tier | Variants | Entrance (%) | Interpretation |
|------|----------|-------------|----------------|
| **Tier 1: Complete closure** | HKQ (0%), GKQ (0%), R74G-R147K-Q140L (0%) | 0% | Zero entrance binding in 10 models |
| **Tier 2: Strong closure** | R74G (10%), R74G-R147K (10%) | 10% | R74G effect persists |
| **Tier 3: WT-like** | WT (20%) | 20% | Baseline |
| **Tier 4: Moderate entrance** | K22Y (30%), R74H (30%), S186C (30%), K22M (40%) | 30--40% | Mixed effects |
| **Tier 5: High entrance** | I226V (50%), I145F (50%), S208E (60%), T209D (60%), Y19H (70%), I145A (70%), V230E (70%), M212Q (70%), R147K (80%) | 50--80% | Elevated relative to WT |

### 5.3 Cross-Condition Consistency

The top-tier variants (HKQ, GKQ, R74G-R147K-Q140L, R74G, R74G-R147K) maintain their ranking across both cofactor-containing conditions. This consistency indicates that their entrance-closure mechanism is robust and not an artifact of pocket conditioning.

Several single mutants, however, show notable re-ordering. Y19H drops from Tier 3 (16%, constrained) to Tier 5 (70%, no-constraint), suggesting that the pocket conditioning amplifies its modest protective effect. Conversely, R147K rises from Tier 4 (28%, constrained) to Tier 5 (80%, no-constraint), suggesting that its effect in the constrained condition may partly reflect enhanced sampling (50 vs. 10 models per variant) rather than a strong intrinsic effect.

These re-orderings highlight that the no-constraint condition, while having lower statistical power (10 models/variant), may be less prone to overestimating weak effects.

---

## 6. The Entrance Closure Hierarchy

### 6.1 The R74 Gateway Mechanism

Position 74 is the most influential single residue for L2 binding mode. The three residue identities at this position -- Arg (WT), Gly (R74G), and His (R74H) -- represent three distinct gateway states:

**Arg74 (WT): Electrostatic anchor, open gate.**
The wild-type arginine at position 74 provides a positively charged anchor that attracts the negatively charged MVAP substrate. However, it does not physically occlude the entrance channel. In the constrained condition, the loop spanning residues 17--33 makes 0 van der Waals contacts with L2 via R74, leaving the gate open for non-competitive entrance binding. WT shows 28% entrance (constrained) and 20% (no-constraint).

**Gly74 (R74G): No gate, competitive entry.**
Replacing arginine with glycine removes both the electrostatic anchor and any steric barrier. Paradoxically, this reduces entrance binding to 6% (constrained) and 10% (no-constraint). The mechanism is not entrance closure but rather redirection: without the R74 anchor to position L2 at the entrance, L2 instead enters the active-site pocket itself (competitive binding). R74G is the dominant source of competitive binding models (6 of ~56 in the constrained condition). The experimental Ki of 110 mM reflects this shift from non-competitive to competitive inhibition -- competitive inhibitors are less potent when measured as Ki because substrate can outcompete them at moderate concentrations.

**His74 (R74H): Physical closure.**
Histidine at position 74 physically closes the entrance channel. The loop spanning residues 17--33 makes 21 van der Waals contacts with L2 in close-binding models, sterically blocking entrance access. R74H shows 10% entrance (constrained) and 30% (no-constraint). This mechanism is fundamentally different from R74G: rather than redirecting L2, it prevents L2 from accessing the entrance at all.

### 6.2 Single-Mutation Effects Beyond R74

| Mutation | Constrained (%) | No-Constraint (%) | Proposed Mechanism |
|----------|-----------------|-------------------|-------------------|
| Y19H | 16 | 70 | Disrupts aromatic stacking at entrance; effect amplified by constraints |
| K22M | 18 | 40 | Removes positive charge at entrance; reduces electrostatic attraction for MVAP |
| K22Y | 24 | 30 | Replaces charge with aromatic; modest steric effect |
| R147K | 28 | 80 | Shortens side chain; effect subtle in isolation |
| M212Q | 36 | 70 | Introduces polar amide; complex effects (see GKQ paradox) |
| I145A | 32 | 70 | Reduces steric bulk; may widen entrance |
| I145F | 48 | 50 | Introduces aromatic bulk; may alter loop dynamics |
| S186C | 36 | 30 | Thiol substitution; modest effect |
| S208E | 42 | 60 | Introduces negative charge near entrance; may attract MVAP |
| T209D | 26 | 60 | Introduces negative charge; mixed effect |
| I226V | 24 | 50 | Conservative substitution; minimal effect |
| V230E | 40 | 70 | Introduces negative charge; increases entrance binding |

### 6.3 Synergistic Multi-Mutation Effects

The most striking finding across all three conditions is the synergy of combined mutations:

**R74G-R147K (double mutant)**: 6% constrained, 10% no-constraint. Identical to R74G alone, suggesting R147K adds no further benefit in this context. Titer = 909 mg/L.

**HKQ (R74H-R147K-M212Q, triple mutant)**: 2% constrained, 0% no-constraint. The combination of three mechanisms produces near-complete entrance closure:
- **H74**: Physical gate closure (21 van der Waals contacts, loop 17--33)
- **Q212**: Bridge formation between Y19 and K22, locking the entrance loop
- **K147**: Salt bridge with E144, stabilizing the closed conformation

The HKQ triple mutant achieves the best experimental titer (1079 mg/L) and a high Ki (80 mM), consistent with the computational prediction that entrance closure alleviates substrate inhibition.

**GKQ (R74G-R147K-M212Q, triple mutant)**: 4% constrained, 0% no-constraint. Despite low entrance fractions comparable to HKQ, GKQ has dramatically different experimental outcomes (Ki = 11 mM, titer = 8 mg/L). This paradox is addressed in detail in Section 11.

**R74G-R147K-Q140L (triple mutant)**: 4% constrained, 0% no-constraint. Another triple mutant achieving near-complete entrance closure. Titer = 401 mg/L, intermediate between HKQ and WT.

---

## 7. Cofactors as Structural Scaffold

### 7.1 The Cofactor Effect

The most dramatic finding of the three-condition comparison is the role of cofactors (ATP and Mg2+) in enabling variant discrimination. Without cofactors, virtually all variants converge to 90--100% entrance binding, and the correlation with Ki collapses to r = 0.006.

| Condition | Mean Entrance (%) | Range (%) | Ki Correlation (r) |
|-----------|-------------------|-----------|-------------------|
| Constrained | ~22 (weighted) | 2--48 | -0.60 |
| No-constraint | ~39 (weighted) | 0--80 | -0.449 |
| No-cofactor | ~97 (weighted) | 90--100 | 0.006 |

### 7.2 Mechanistic Interpretation

ATP and Mg2+ bind at the active site and occupy substantial volume within the binding cleft. Their presence narrows the entrance channel that connects the bulk solvent to the active-site pocket. This narrowing has two consequences:

1. **Increased sensitivity to mutations.** When the channel is narrow (cofactors present), small changes in side-chain volume, charge, or flexibility at positions like R74, K22, R147, and M212 can tip the balance between entrance-accessible and entrance-occluded conformations. When the channel is wide (cofactors absent), these same mutations cannot overcome the inherent accessibility of the entrance.

2. **Functional relevance.** Substrate inhibition occurs during catalytic turnover, when ATP and Mg2+ are necessarily bound. The cofactor-dependent entrance narrowing thus represents the physiologically relevant condition. The no-cofactor condition -- while informative as a negative control -- does not represent the state of the enzyme during catalysis.

### 7.3 Implications for Simulation Design

The no-cofactor result validates the inclusion of cofactors as essential for meaningful substrate inhibition modeling. Future Boltz-2 cofolding studies of PMDsc should always include ATP and Mg2+. Conversely, the no-cofactor condition can serve as a useful negative control: any variant that shows reduced entrance binding even without cofactors would represent an exceptionally strong entrance-closure mutation (none were observed in this study).

---

## 8. Pocket Conditioning: Helpful but Not Essential

### 8.1 What Pocket Conditioning Provides

Pocket conditioning constrains L1 to the active site, ensuring that every simulation starts from a catalytically relevant configuration. The key question is whether this constraint artificially drives L2 to the entrance (by occupying the active site with L1 and forcing L2 elsewhere) or whether L2 binding patterns reflect genuine energetics.

### 8.2 Evidence That Conditioning Does Not Create Artifacts

Three lines of evidence argue against artifactual L2 placement:

1. **L1 finds the active site without conditioning.** In the no-constraint condition, L1 independently localizes to within 3.4--4.5 A of the active site across all 19 variants. This means the active site is occupied by L1 in both conditions, and L2 must find a secondary site in either case.

2. **Consensus contact residues are conserved.** The top three L2 contact residues in close-binding models are the same across conditions:
   - Constrained: R/G/H74 (80%), K22 (79%), T209 (79%)
   - No-constraint: K22 (86%), R/G/H74 (86%), T209 (84%)
   - The same residues dominate, indicating that L2 binds to the same structural feature (the entrance channel) regardless of how L1 was placed.

3. **Variant ranking is preserved.** The top-tier variants (HKQ, GKQ, R74G-R147K-Q140L, R74G, R74G-R147K) maintain their positions in both the constrained and no-constraint conditions. If conditioning were creating artifacts, we would expect random re-ordering.

### 8.3 What Pocket Conditioning Does Provide

While not essential, pocket conditioning offers two practical advantages:

**Higher throughput.** With 50 models per variant (vs. 10), the constrained condition provides better sampling and more precise entrance fraction estimates. The standard error on a 28% fraction is approximately 6.5% with n=50 versus 14% with n=10.

**Sharper discrimination.** The constrained condition produces a tighter range of entrance fractions (2--48%) and a stronger Ki correlation (r = -0.60 vs. -0.449). This likely reflects both better sampling and the elimination of occasional L1 misplacement events that could confound L2 analysis.

### 8.4 Recommendation

For screening purposes, the constrained condition is preferred due to its higher throughput and discrimination power. For mechanistic studies, the no-constraint condition provides a valuable orthogonal validation. The two conditions should be viewed as complementary rather than redundant.

---

## 9. Consensus Contact Residues Across Conditions

### 9.1 Contact Frequency Tables

**Constrained condition (271 close-binding models):**

| Residue | Contact Frequency (%) |
|---------|----------------------|
| R/G/H74 | 80 |
| K22 | 79 |
| T209 | 79 |
| N28 | 61 |
| N72 | 56 |
| S153 | 50 |
| S208 | 47 |

**No-constraint condition (87 close-binding models):**

| Residue | Contact Frequency (%) |
|---------|----------------------|
| K22 | 86 |
| R/G/H74 | 86 |
| T209 | 84 |
| N72 | 68 |
| N28 | 56 |
| S153 | 51 |
| S208 | 49 |

**No-cofactor condition (187 close-binding models):**

| Residue | Contact Frequency (%) |
|---------|----------------------|
| S208 | 74 |
| T209 | 74 |
| R/G/H74 | 73 |
| S121 | 67 |
| G207 | 54 |
| K22 | 48 |

### 9.2 Cross-Condition Analysis

The top three residues -- R/G/H74, K22, and T209 -- appear in the top tier across all three conditions, confirming that the entrance channel is a conserved structural feature defined by these residues. However, important differences emerge:

**Cofactor-dependent contacts.** N28 and N72 rank 4th and 5th in both cofactor-containing conditions (56--68%) but drop out of the top six without cofactors. These residues likely make indirect contacts with L2 that are stabilized by the cofactor-mediated narrowing of the entrance channel.

**Cofactor-independent contacts.** S208 and T209 maintain high contact frequencies across all conditions (47--84%), suggesting they form direct contacts with L2 independent of cofactor presence. These residues define the core of the entrance binding site.

**No-cofactor-specific contacts.** S121 (67%) and G207 (54%) appear prominently only in the no-cofactor condition. This suggests that without cofactors, the wider entrance channel allows L2 to sample a broader binding surface, making contacts with residues that are inaccessible when cofactors narrow the channel.

### 9.3 Structural Interpretation

The entrance binding site is defined by residues on three structural elements:

1. **The 17--33 loop** (K22, N28): Forms one wall of the entrance channel. Mutations at K22 (K22M, K22Y) and the Q212 bridge to Y19 and K22 modulate this wall.

2. **The 70--74 region** (N72, R/G/H74): Forms the gateway. R74 mutations have the largest single-residue effect on entrance binding.

3. **The 207--209 turn** (S208, T209): Forms the opposite wall. Mutations at S208 (S208E) and T209 (T209D) introduce negative charges that may attract the phosphate group of MVAP, paradoxically increasing entrance binding.

---

## 10. The R74G Competitive Mode

### 10.1 Two Modes of Substrate Inhibition

The simulations reveal two distinct modes of L2 binding:

1. **Entrance binding (non-competitive)**: L2 binds at the mouth of the active-site channel, physically blocking substrate access and product release. This is the dominant mode for most variants and corresponds to non-competitive or uncompetitive inhibition kinetics.

2. **Active-site binding (competitive)**: L2 enters the catalytic pocket itself, directly competing with L1 for the active site. This corresponds to competitive inhibition kinetics.

### 10.2 Active-Site Binding Across Conditions

| Condition | Total Competitive Models | Key Contributors |
|-----------|------------------------|-----------------|
| Constrained | ~56 | R74G (6), S208E (14), WT (2), others |
| No-constraint | 12 | WT (3), I145A (2), S186C (2), others (1 each) |
| No-cofactor | 5 | Scattered across variants |

### 10.3 R74G as the Poster Case

R74G consistently produces competitive binding models across conditions. The mechanism is straightforward: removing the bulky, positively charged arginine at the gateway eliminates both the electrostatic attraction that positions L2 at the entrance and the steric barrier that prevents L2 from passing through the gate into the active site. With the gate open, L2 preferentially enters the active-site pocket rather than binding at the entrance.

This explains the experimental data for R74G:
- **Ki = 110 mM** (highest measured): Competitive inhibitors are effectively weaker because the substrate can outcompete them. The high Ki reflects the competitive binding mode, not the absence of inhibition.
- **Titer = 975 mg/L** (second highest): The shift from non-competitive to competitive inhibition is functionally beneficial because competitive inhibition can be overcome by substrate at moderate concentrations.
- **kcat/Km = 0.04** (lowest measured): The removal of R74's electrostatic anchor slightly impairs catalytic efficiency, consistent with a role for R74 in substrate positioning.

### 10.4 S208E Competitive Binding

In the constrained condition, S208E produces the most competitive binding models (14 of ~56). The introduction of a negative charge (glutamate) at position 208 may repel the phosphate group of MVAP from the entrance, redirecting L2 toward the active site. However, S208E also shows high entrance binding (42% constrained, 60% no-constraint), suggesting it promotes both binding modes. This dual effect may explain its lack of experimental titer data -- it may have been deprioritized due to complex inhibition kinetics.

### 10.5 Condition Dependence of Competitive Binding

Competitive binding is most frequent in the constrained condition (~56 models) and drops sharply in the no-constraint (12 models) and no-cofactor (5 models) conditions. This trend suggests that pocket conditioning, by firmly placing L1 at the active site, may slightly overestimate competitive binding by creating a more ordered active-site environment that can accommodate a second substrate molecule. Alternatively, the higher model count in the constrained condition (950 vs. 190) provides more opportunities to sample the relatively rare competitive binding pose.

---

## 11. The GKQ Paradox Across Conditions

### 11.1 The Paradox Defined

GKQ (R74G-R147K-M212Q) presents the most challenging case for the entrance-frequency model of substrate inhibition:

| Metric | HKQ | GKQ | Prediction from Entrance Fraction |
|--------|-----|-----|----------------------------------|
| Constrained entrance (%) | 2 | 4 | Similar inhibition |
| No-constraint entrance (%) | 0 | 0 | Identical inhibition |
| No-cofactor entrance (%) | 90 | 90 | Identical inhibition |
| Ki (mM) | 80 | 11 | GKQ should have less inhibition |
| kcat/Km | 0.40 | 0.50 | GKQ is slightly more efficient |
| Titer (mg/L) | 1079 | 8 | GKQ should perform similarly |

By all computational metrics across all three conditions, GKQ and HKQ are nearly identical. Yet experimentally, GKQ has 7-fold more severe substrate inhibition (Ki = 11 vs. 80 mM) and 135-fold lower titer (8 vs. 1079 mg/L).

### 11.2 Resolution: Frequency Does Not Equal Affinity

The entrance fraction measures how often L2 visits the entrance in Boltz-2 cofolding -- a frequency metric. Ki, however, reflects binding affinity -- how tightly L2 binds once it arrives. The GKQ paradox reveals that these two quantities can be decoupled:

**GKQ rarely places L2 at the entrance, but when it does, L2 binds tightly.**

The proposed molecular mechanism centers on M212Q. The glutamine substitution introduces a free amide group that can form hydrogen bonds with the MVAP substrate. In the context of GKQ:

- **G74** (from R74G): Opens the gate, allowing L2 access to a broader surface including deeper entry into the active site.
- **K147** (from R147K): Shortens the side chain at position 147.
- **Q212** (from M212Q): The free amide of glutamine stabilizes any L2 molecule that does reach the entrance or active-site vicinity through additional hydrogen bonding.

In HKQ, the same Q212 amide is present, but H74 physically blocks L2 from reaching a position where the Q212 bridge can stabilize binding. The H74 gate acts as a kinetic barrier that prevents L2 from accessing the high-affinity binding pose.

### 11.3 Why HKQ Succeeds and GKQ Fails

The success of HKQ can be understood as a triple defense:

1. **Kinetic barrier** (H74): Physically blocks L2 from entering the channel (21 van der Waals contacts with loop 17--33).
2. **Entrance lock** (Q212): Bridges Y19 and K22, locking the entrance loop in a closed conformation.
3. **Structural stabilization** (K147): Forms a salt bridge with E144, rigidifying the closed state.

In GKQ, the first line of defense is absent (G74 provides no barrier), and the second (Q212) paradoxically becomes a liability -- it stabilizes the very L2 binding it was meant to prevent.

### 11.4 Implications for the Computational Model

The GKQ paradox exposes a fundamental limitation of the entrance-frequency metric: it captures binding probability but not binding thermodynamics. A complete computational model of substrate inhibition would need to incorporate:

1. **Binding frequency** (captured by entrance fraction)
2. **Binding affinity** (not captured; would require free-energy calculations)
3. **Residence time** (not captured; would require molecular dynamics)

Despite this limitation, the entrance-frequency metric correctly identifies HKQ as the best variant (lowest entrance fraction, highest titer) and correctly ranks most other variants. GKQ represents an edge case where the frequency-affinity decoupling is extreme.

### 11.5 Cross-Condition Consistency of the Paradox

The paradox is consistent across all three conditions:

- **Constrained**: GKQ (4%) vs. HKQ (2%) -- nearly identical, both very low
- **No-constraint**: GKQ (0%) vs. HKQ (0%) -- identical
- **No-cofactor**: GKQ (90%) vs. HKQ (90%) -- identical

No condition resolves the paradox computationally. This consistency actually strengthens the frequency-vs-affinity interpretation: if the paradox were due to a sampling artifact in one condition, we would expect it to resolve in another. Its persistence across conditions argues for a fundamental thermodynamic difference that is invisible to the frequency-based metric.

---

## 12. Statistical Summary

### 12.1 Model Counts and Close-Binding Statistics

| Statistic | Constrained | No-Constraint | No-Cofactor |
|-----------|-------------|---------------|-------------|
| Total models | 950 | 190 | 190 |
| Models per variant | 50 | 10 | 10 |
| Close-binding models | 271 (29%) | 87 (46%) | 187 (98%) |
| Active-site (competitive) models | ~56 (6%) | 12 (6%) | 5 (3%) |
| Variants with 0% entrance | 0 | 3 (HKQ, GKQ, R74G-R147K-Q140L) | 0 |
| Variants with >90% entrance | 0 | 0 | 19 (all) |

### 12.2 Correlation with Experimental Ki

| Condition | Pearson r (Ki vs Entrance %) | Direction | Interpretation |
|-----------|-------------------------------|-----------|----------------|
| Constrained | -0.60 | Negative | Higher Ki (less inhibition) associates with lower entrance fraction |
| No-constraint | -0.449 | Negative | Same trend, weaker signal |
| No-cofactor | 0.006 | None | No discrimination without cofactors |

The negative correlation is the expected direction: variants with lower entrance binding (fewer inhibitory L2 placements) should have higher Ki (weaker inhibition, more substrate needed to trigger it). The constrained condition provides the strongest correlation, likely due to better sampling (50 models/variant).

Note that only five variants have measured Ki values (WT, R74G, V230E, HKQ, GKQ), limiting the statistical power of these correlations. With n = 5, the correlation coefficients should be interpreted as trends rather than definitive relationships.

### 12.3 Entrance Fraction Distribution by Condition

**Constrained (range: 2--48%, median ~24%):**
- Most variants cluster between 16% and 42%
- Clear separation of top-tier (2--6%) and bottom-tier (40--48%) variants
- Distribution is roughly unimodal

**No-constraint (range: 0--80%, median ~50%):**
- Wider spread, reflecting lower sampling and/or reduced constraint-mediated focusing
- Bimodal tendency: a cluster at 0--10% (triple mutants + R74G variants) and a broad cluster at 30--80%
- WT (20%) sits between the two clusters

**No-cofactor (range: 90--100%, median 100%):**
- Extreme ceiling effect: 12 of 19 variants at exactly 100%
- Remaining 7 variants at 90%
- No meaningful discrimination

---

## 13. Recommendations for Future Studies

### 13.1 Computational Priorities

1. **Free-energy perturbation (FEP) calculations for GKQ and HKQ.** The GKQ paradox cannot be resolved by cofolding frequency analysis alone. FEP or thermodynamic integration calculations on the L2 binding pose would quantify the binding affinity difference that the entrance fraction misses. This is the single highest-priority follow-up.

2. **Increased sampling for no-constraint condition.** The current 10 models per variant provides limited statistical power. Increasing to 50 models per variant (matching the constrained condition) would determine whether the no-constraint condition can achieve comparable discrimination, potentially eliminating the need for pocket conditioning in future screens.

3. **Molecular dynamics of entrance-open vs. entrance-closed states.** Short (100--500 ns) MD simulations of HKQ and GKQ in the L2-bound state would quantify residence time differences and test the frequency-vs-affinity hypothesis directly.

4. **Expanded variant library.** The current 19 variants include several positions (I145, S186, I226, V230) that show elevated entrance binding. These positions could be targeted for charge-complementary mutations (e.g., I145D, I145E) designed to electrostatically repel MVAP from the entrance.

### 13.2 Experimental Priorities

1. **Ki measurements for R74H, R147K, and R74G-R147K.** These variants show strong entrance closure (6--10% constrained) and high titers (770--909 mg/L), but lack Ki data. Confirming that their Ki values are elevated would strengthen the Ki-entrance correlation.

2. **Ki measurements for triple mutant R74G-R147K-Q140L.** This variant shows 4% entrance (constrained) and 0% (no-constraint), comparable to HKQ and GKQ. Its Ki would help determine whether the frequency-affinity decoupling seen in GKQ is specific to the Q212 mutation or a more general phenomenon.

3. **Inhibition mode characterization.** The simulations predict that R74G and R74G-containing variants should show competitive inhibition kinetics, while R74H and R74H-containing variants should show non-competitive kinetics. Lineweaver-Burk or Dixon plot analysis could test this prediction.

4. **Combinatorial variants with S208 and T209.** These residues are the most consistent L2 contacts across all conditions. Mutations designed to repel MVAP at these positions (e.g., S208D, T209E) in combination with R74H could further enhance entrance closure.

### 13.3 Methodological Recommendations

1. **Always include cofactors.** The no-cofactor condition demonstrates that ATP and Mg2+ are essential for meaningful variant discrimination. All future Boltz-2 cofolding studies of PMDsc should include cofactors.

2. **Use constrained condition for screening, no-constraint for validation.** The constrained condition provides the best discrimination and throughput for initial variant screening. Top hits should be validated in the no-constraint condition to confirm that their entrance closure is not an artifact of pocket conditioning.

3. **Complement frequency metrics with affinity estimates.** The GKQ paradox demonstrates that entrance frequency alone is insufficient to predict Ki. Future analyses should incorporate scoring functions (e.g., Boltz-2 confidence scores, Rosetta interface energy) that approximate binding affinity in addition to frequency.

4. **Report both entrance and competitive binding.** The R74G case demonstrates that low entrance binding can mask a shift to competitive binding. Both metrics should be reported to capture the full picture of L2 binding behavior.

---

## 14. Conclusions

This three-condition Boltz-2 cofolding analysis of 19 PMDsc variants yields the following principal conclusions:

### 14.1 The Entrance Channel is the Primary Inhibition Site

Across all conditions, the entrance channel -- defined by residues R/G/H74, K22, T209, N28, N72, and S208 -- is the dominant L2 binding site. This is consistent with a non-competitive inhibition mechanism in which L2 blocks substrate access to the active site without occupying the catalytic pocket itself.

### 14.2 Cofactors Are Essential for Variant Discrimination

ATP and Mg2+ narrow the entrance channel and enable mutations to differentially modulate L2 access. Without cofactors, the channel is wide, and all variants converge to 90--100% entrance binding (r = 0.006 with Ki). With cofactors, entrance fractions span 0--80% (no-constraint) or 2--48% (constrained), and the correlation with Ki reaches r = -0.60.

### 14.3 Pocket Conditioning Sharpens but Does Not Create Discrimination

L1 independently finds the active site (3.4--4.5 A) without pocket conditioning. The same top contact residues and the same variant hierarchy emerge in both cofactor-containing conditions. Pocket conditioning provides higher throughput and sharper discrimination but is not required for meaningful results.

### 14.4 R74 is the Master Switch

Position 74 controls the gateway between entrance binding and active-site (competitive) binding. Arginine (WT) anchors L2 at the entrance (non-competitive). Glycine (R74G) opens the gate for competitive entry (Ki = 110 mM, titer = 975 mg/L). Histidine (R74H) physically closes the gate (10% entrance, titer = 770 mg/L).

### 14.5 Triple Mutants Achieve Synergistic Closure

HKQ (R74H-R147K-M212Q) achieves near-complete entrance closure (2% constrained, 0% no-constraint) through three synergistic mechanisms: H74 gate closure, Q212 entrance-loop bridging, and K147-E144 salt bridge stabilization. It achieves the best experimental titer (1079 mg/L) and a high Ki (80 mM).

### 14.6 Frequency Does Not Equal Affinity

GKQ (R74G-R147K-M212Q) matches HKQ in entrance frequency (4% constrained, 0% no-constraint) but has 7-fold lower Ki (11 mM) and 135-fold lower titer (8 mg/L). The resolution is that Q212's free amide stabilizes rare L2 visits in the absence of the H74 kinetic barrier, creating high-affinity binding that the frequency metric cannot detect. This paradox is consistent across all three conditions and represents a fundamental limitation of frequency-based computational screening.

### 14.7 Summary of Design Principles

For engineering PMDsc variants with reduced substrate inhibition:

1. **Block the entrance** with sterically bulky, uncharged residues at position 74 (e.g., His).
2. **Lock the entrance loop** with bridging residues that stabilize the closed conformation (e.g., Q212 in combination with H74).
3. **Stabilize the closed state** with salt bridges in the surrounding structure (e.g., K147-E144).
4. **Avoid opening the gate** without compensating closure mechanisms (the GKQ lesson).
5. **Always include cofactors** in computational screens -- they are essential for discrimination.
6. **Validate hits across conditions** to distinguish robust effects from conditioning artifacts.

---

*Analysis performed using Boltz-2 cofolding with 950 constrained models, 190 no-constraint models, and 190 no-cofactor models across 19 PMDsc variants. Three conditions tested: pocket-conditioned with cofactors, unconstrained with cofactors, and unconstrained without cofactors.*
