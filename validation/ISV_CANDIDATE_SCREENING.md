# ISV — Initial scientific validation candidate screening (v1)

**Status:** candidate sourcing only; **no candidate has passed final evidence-to-structure acceptance**. This is Dataset 2 of the [two-dataset plan](ACTIVE_DATASETS.md). The technical dataset TV48 is unchanged.

**Version:** 23 November 2026. [Machine-readable candidate list](ISV_candidate_screening_v1.csv).

## Current screening pool

Eight **distinct protein candidates** plus one nonindependent backup PDB structure (4EPC for R2ab). No inclusion count, experimental localization performance, or independence of these cases is being claimed. The first six protein candidates are the priority source-audit queue; two have secondary/indirect evidence and are *not* valid primary residue-contact ground truth without further evidence.

| Priority | Protein | Proposed PDB | Material | Source DOI | Evidence class | Key remaining qualification |
| --- | --- | --- | --- | --- | --- | --- |
| High | Acylphosphatase | [2ACY](https://www.rcsb.org/structure/2ACY) | Silica | [10.1021/la401985d](https://doi.org/10.1021/la401985d) | NMR + proteomics localization | Verify experimental AcP isoform/species vs bovine testis 2ACY; map precise residues |
| High | RNase A | [7RSA](https://www.rcsb.org/structure/7RSA) | Silica | [10.1021/nl2044524](https://doi.org/10.1021/nl2044524) | Chemical modification + proteomics | Original modification site list and matching bovine chain |
| High | Hen egg-white lysozyme | [1HEL](https://www.rcsb.org/structure/1HEL) | Silica | Same study | Chemical modification + proteomics | Original residue evidence; independently count same-study conditions correctly |
| High | Cytochrome c | [1HRC](https://www.rcsb.org/structure/1HRC) | Silica | Same study | Chemical modification + proteomics | Source organism, residue mapping, heme state and ligand policy |
| High | *S. epidermidis* R2ab | [7KWI](https://www.rcsb.org/structure/7KWI) | Polystyrene | [10.1021/acs.langmuir.3c02609](https://doi.org/10.1021/acs.langmuir.3c02609) | Residue-specific chemical modification / proteomics | Confirm variant and NMR ensemble selection; chemical changes may reflect unfolding, not direct contact |
| Medium | GB3, protein G domain | [2OED](https://www.rcsb.org/structure/2OED) | Gold NP | [10.1038/s41467-022-34749-w](https://doi.org/10.1038/s41467-022-34749-w) | NMR, mutations and orientation/affinity | Do not mislabel mutation/affinity changes as observed atomic contact residues |
| Secondary | Myoglobin | [1A6M](https://www.rcsb.org/structure/1A6M) | Silica | [10.1016/S0021-9797(03)00401-6](https://doi.org/10.1016/S0021-9797(03)00401-6) | HDX-MS structural change | Experiment-to-PDB species, ligand state and whether labels describe stability rather than contact |
| Secondary | Proteinase K | [5B1E](https://www.rcsb.org/structure/5B1E) | Gold NP | Same 2022 Nature Communications study | Functional accessibility | Activity readout is not experimental localization ground truth |
| Backup only | R2ab alternative partial repeat structure | [4EPC](https://www.rcsb.org/structure/4EPC) | Polystyrene | Same R2ab study | Structural alternative | Not a ninth independent protein; align experimental construct to 7KWI and 4EPC |

Source studies for the first five candidates and the GB3 orientation approach were verified to exist and to address these types of measurements. This does **not** establish exact residue labels or correct cross-study PDB mapping.

## Independence and leakage controls

1. **TV48 overlap:** none of the PDB IDs in this screening table appears verbatim among the 48 frozen technical PDB identifiers. This is **ID-level nonoverlap only**. Protein-domain homology, repeated structural families and prior exploratory use must be audited before claiming an independent holdout.
2. **Same-study clustering:** RNase A, lysozyme and cytochrome c share experimental publication [10.1021/nl2044524](https://doi.org/10.1021/nl2044524) and must not count as three independent publications. GB3 and proteinase K share a second publication. R2ab/4EPC are structural alternatives for the same protein and must not count twice.
3. **Historical overlap:** check every proposed protein against older nine-condition labels and any other previously inspected literature/localization examples. No outcome-driven channel or parameter selection.
4. **Evidence resolution:** NMR perturbations, footprinting, mutations, fluorescence, proteolysis and orientation restraints are different measures. Retain their methodological uncertainty; a change after adsorption may not be a direct contact.
5. **Matched input:** source species, isoform, variant, pH, buffer, nanoparticle chemistry, construct, relevant cofactors, numbering and experimentally observed regions must match the selected structure.
6. **Fair baseline:** preregister channels and region-size/coverage comparators; compare only supportable, predeclared metrics and cluster uncertainty by independent protein/publication.

## Admission stages

**CANDIDATE → ORIGINAL_SOURCE_AUDITED → STRUCTURE_MATCHED → NUMERIC_QA_PASSED → ACCEPTED.**

The current entries have not passed the second, third or fourth stage. Specific evidence ranges are *deliberately omitted* rather than invented from abstracts.

### Next source-audit outputs

For each of the six priority proteins, create a detailed row that lists DOI, experimental section/figure/SI, exact measured positions (original sequence numbering), evidence type, sequence to structure alignment, missing structural residues and treatment of the material/measurement conditions. Only then choose **ACCEPTED**, **CONDITIONAL**, or **REJECTED**, and freeze an independent `ISV_manifest.csv`.

In contrast with the existing 48-protein technical verification, **no ISV result has been computed or claimed**.
