# PMD specific proton transfer evidence

Research date: 2026-11-13. Scope: PMD/MPD with mevalonate 5-monophosphate (MVAP/MVA5P), especially Anaerolinea thermophila and Herpetosiphon aurantiacus; distinguish this from MDD/DMD with MVAPP/MVA5PP. This is a focused evidence search, not a claim of exhaustive systematic coverage. Searches used primary publisher/PubMed/PMC sources, the user-provided Thomas2019 and VanNice2014 PDFs, and primary full-text XML where PMC browser access was blocked.

## What does the direct At and Ha evidence establish?

### Takeaway

The central experimental study of the user's two PMD enzymes establishes substrate preferences, substrate-binding interactions and enzymatic activity. It does not establish the initial C3-OH proton acceptor. There is no experimentally justified At/Ha-specific choice between Asp-mediated and ATP-gamma-oxygen-mediated proton transfer in this paper.

### Cited Findings

- Thomas et al. 2019 characterized both enzymes with MVAP and MVAPP. At is substrate-promiscuous; Ha showed no measurable MVAPP turnover under the reported conditions. Thus an MVAPP-specific MDD mechanism is an analogy for these PMDs, not a demonstrated identity. — [Thomas et al., ACS Chemical Biology 2019, DOI 10.1021/acschembio.9b00322](https://pubs.acs.org/doi/10.1021/acschembio.9b00322), [PubMed](https://pubmed.ncbi.nlm.nih.gov/31268677/).
- At structures 6N0X and 6N0Y contain MVAP and MVAPP respectively. Attempts to cocrystallize with ADP gave no additional nucleotide-pocket density. These structures therefore do not directly resolve an ATP-containing proton-transfer arrangement. — [Thomas2019](https://pubs.acs.org/doi/10.1021/acschembio.9b00322), [6N0X](https://www.rcsb.org/structure/6N0X), [6N0Y](https://www.rcsb.org/structure/6N0Y).
- The paper identifies conserved Arg143/Asp280 in At as catalytic residues, citing earlier decarboxylase studies; this is not a new measurement isolating proton abstraction during At phosphorylation. — [Thomas2019](https://pubs.acs.org/doi/10.1021/acschembio.9b00322).
- At H194M changes kcat insignificantly but increases KM for MVAP 24-fold. Ha H204M increases KM for MVAP approximately 50-fold. The reported histidine role is primarily substrate-phosphate recognition. — [Thomas2019](https://pubs.acs.org/doi/10.1021/acschembio.9b00322).

Full-text verification: local user-supplied PDF `/nfs/lambda_stor_01/homes/rzhu/PMD/paper/substrate-specificity-and-engineering-of-mevalonate-5-phosphate-decarboxylase.pdf`, extracted with PyMuPDF. Relevant printed pages: 1768–1770 (enzyme characterization and structures), 1772–1774 (specificity and mutations). Whole paper was searched for deprotonation/proton/acid/base/Asp/His/Arg language; its focus is specificity, not a resolved initial proton-transfer mechanism.

### Inferences

- A nearby His in a generated pocket should not automatically be labeled a general base. In the best direct PMD study it is a phosphate-binding determinant, and losing it predominantly changes substrate affinity.
- An Asp near MVAP is a plausible catalytic-site feature, but neither its proximity nor its conservation identifies its chemical task. It could orient substrate, participate in proton transfer, affect electrostatics, or contribute to subsequent elimination/decarboxylation.
- The lack of a nucleotide-bound At crystal structure in Thomas2019 matters: a predicted ATP-Mg pose supplies a hypothesis that has not been observed there at atomic resolution.

### Gaps

- No At/Ha-specific positional proton-transfer measurement, pre-steady-state separation of kinase and decarboxylase chemistry, solvent isotope/proton inventory that identifies the acceptor, or experimental assignment of the accepting ATP oxygen was found in the examined direct study.
- No verified later primary At/Ha study resolving this particular acceptor question was found in focused searches of 2023–2026 literature. Absence in this search does not prove that none exists.

## What independent PMD experiments constrain the question?

### Takeaway

There is experimental PMD evidence beyond the one MDD QM/MM paper, but it constrains reaction architecture and ionization rather than identifying the receiving atom. The most relevant PMD result is selective loss of phosphorylation at low pH while decarboxylation of a prephosphorylated intermediate remains active.

### Cited Findings

- Dellas et al. 2013 demonstrated an ATP-dependent MVAP-to-IP activity in Roseiflexus castenholzii, establishing the alternative pathway. It is a pathway/product-discovery result rather than an assignment of the C3-OH proton acceptor. — [Dellas et al., eLife 2013, DOI 10.7554/eLife.00672](https://elifesciences.org/articles/00672), [PMC3857490](https://pmc.ncbi.nlm.nih.gov/articles/PMC3857490/).
- VanNice et al. 2014 characterized Haloferax volcanii PMD and verified its product using coupled enzymatic tests and mass spectrometry. Fluorinated substrate inhibition supported similarity to classical MDD chemistry, without establishing who accepts the initial proton. Assays used MgCl2; an assay metal requirement is not itself a resolved active-site coordination model. — [VanNice et al., Journal of Bacteriology 2014, DOI 10.1128/JB.01230-13](https://pmc.ncbi.nlm.nih.gov/articles/PMC3957691/).
- Vinokur et al. 2016 supplied Roseiflexus PMD with its normal MVAP substrate or the preformed 3,5-bisphosphorylated intermediate. At pH 4.3, normal substrate conversion stopped while decarboxylation of the intermediate remained. Authors suggested altered ionization of an active-site residue involved in C3-OH deprotonation as **one possible explanation**, not a measured proton-acceptor identity. — [Vinokur et al., Scientific Reports 2016, DOI 10.1038/srep39737](https://pmc.ncbi.nlm.nih.gov/articles/PMC5177888/), [publisher](https://www.nature.com/articles/srep39737).
- Aoki et al. 2022 studied ATP-independent Picrophilus torridus mevalonate 3,5-bisphosphate decarboxylase (MBD), which performs the downstream chemistry. Mutation of conserved Asp309 to Thr or Asn abolished MBD activity. Its structure places Asp309 analogously to AtPMD Asp280, supporting an important Asp role in downstream chemistry without establishing an initial PMD proton acceptor. — [Aoki et al., Journal of Biological Chemistry 2022, DOI 10.1016/j.jbc.2022.102111](https://pmc.ncbi.nlm.nih.gov/articles/PMC9254496/).

Vinokur full text was independently read as primary XML from `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5177888/fullTextXML` because browser access to PMC intermittently returned a challenge page. The paper's discussion explicitly frames the active-site-residue interpretation as a possible explanation. Aoki2022 is a distinct decarboxylase and must not be described as direct replication of ATP-mediated proton transfer.

### Inferences

- The PMD pH observation is compatible with a protein-base requirement and should be included as a counterweight to an overconfident gamma-phosphate narrative. It is also compatible with changes in substrate/ATP phosphate protonation, metal coordination, and other coupled ionizations. It does not distinguish them.
- Conserved Asp being required for ATP-independent downstream decarboxylation helps explain why its mutational effects on full turnover cannot by themselves establish a first-step general-base role.
- Several studies supporting the two-stage reaction do not equal several direct observations of ATP-gamma oxygen accepting the proton. Count experiments by the specific claim they test.

### Gaps

- No PMD experiment in these sources directly tracks the hydroxyl proton onto ATP phosphate or Asp.
- Salt-adapted H. volcanii PMD, Roseiflexus PMD, AtPMD and HaPMD are distinct enzymes with different assay environments; their behavior should not be pooled as a single enzyme-specific experiment.
- MBD is ATP-independent and has no initial phosphorylation. Its Asp result cannot exclude a separate Asp role in a dual-function enzyme's first step.

## What can be concluded for the RFD4 designs and what would discriminate routes?

### Takeaway

For At/Ha PMD and the RFD4 structures, the initial proton-acceptor mechanism remains unresolved by direct enzyme-specific evidence. ATP-gamma-assisted transfer is a credible working hypothesis supported in related MDD literature, while an Asp-assisted route should remain a live alternative. This is an evidence assessment, not a quantified probability.

### Cited Findings

- Both acceptor proposals concern the same initial substrate C3-OH phosphorylation, followed by distinct downstream elimination/decarboxylation chemistry. The step-specific distinction is established in the PMD intermediate/pH study. — [Vinokur2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5177888/).
- The direct At/Ha study demonstrates substantial substrate-recognition differences even between homologous decarboxylases; this limits unqualified transfer of a computed MDD/MVAPP mechanism to PMD/MVAP. — [Thomas2019](https://pubmed.ncbi.nlm.nih.gov/31268677/).
- The related MDD structural/mechanistic study discusses Asp as an orienting/downstream residue, but treats the identity of the proton acceptor as unresolved and phosphate acceptance as possible. Its model does not image the transferring hydrogen. — [Chen et al., Nature Communications 2020, DOI 10.1038/s41467-020-17733-0](https://pmc.ncbi.nlm.nih.gov/articles/PMC7414129/).

### Inferences

- Plot the substrate hydroxyl O to both Asp carboxylate oxygens, and to **each** of ATP's three nonbridging gamma oxygens. Showing the minimum gamma distance is useful only if the individual oxygen identities and all three values are retained. In the RFD4 naming scheme these are substrate IMV O16 versus ATP IAT O13/O16/O17, not the attacking phosphorus IAT P41.
- These O···O distances test access, not proton-transfer favorability. A short contact could be a hydrogen bond with the wrong donor/acceptor direction; proton orientation, phosphate/Asp protonation and Mg binding can change its significance. The shorter distance is not evidence of the lower activation barrier.
- Use the common substrate C3-O/Pgamma attack coordinate and leaving-bond geometry for both hypotheses; otherwise an apparent preference might be caused by incomparable reaction definitions.
- A useful experimental discrimination strategy would measure early ADP formation and the 3-phosphorylated intermediate separately from final IP, compare carefully chosen Asp variants with WT, verify retained fold and substrate binding, and examine pH/solvent isotope behavior. No individual mutation or isotope effect alone uniquely identifies an accepting atom; interpretation needs a kinetic scheme and structural controls.
- For mechanistic calculations, compare both pathways from matched protonation, hydration, Mg coordination, enzyme/substrate identity and conformational ensembles. Explicitly follow the O-H bond and candidate acceptor-H bonds together with phosphoryl transfer. Validate endpoint chemistry and, when claiming a transition state, the relevant mode and connectivity. A one-distance scan or a static pose cannot resolve the mechanism.
- Honest colleague-facing wording: "The RFD4 pockets provide candidate geometries compatible with both proposed proton-acceptor routes. Literature supports ATP-phosphate-assisted transfer as a plausible hypothesis in related MDD, but direct PMD-specific identification of the initial proton acceptor is lacking. We therefore compare both routes and do not infer mechanism or enzyme superiority from distance alone."

### Gaps

- No model-specific reaction-path or experimental proton-transfer evidence exists in the supplied RFD4 coordinate files. Their local contacts are structural hypotheses.
- PMD-specific pKa assignments, metal-bound solvent occupancy and free-energy comparisons of the two routes remain unresolved in this source set.
- Results concerning downstream Asp catalysis should not be used to remove Asp from a design simply because a gamma-assisted first step is plausible; the full enzyme still needs downstream chemistry.
