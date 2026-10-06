# Methodology and provenance

## Scope

The target is an independent experimental benchmark for protein–material and protein–nanoparticle interfaces. The benchmark is designed to test whether InterfaceScout's material-chemistry-conditioned structural descriptors associate with observed adsorption/corona occurrence, quantitative abundance, protected/contact-compatible regions, and adsorption-associated conformational response.

## Source discovery

Searches covered protein-corona databases, biomaterial adsorption databases, nanoparticle-protein proteomics repositories, limited-proteolysis/LiP-MS studies, NMR/CSP mapping, crosslink/covalent-capture MS, solid-state NMR, PRIDE datasets, GitHub repositories, and publication supplementary data. Each candidate source is tracked at the independent experimental-family level to avoid double-counting a database and a later re-analysis of the same records.

## Data ingestion

Raw source tables are copied without modifying source values. Harmonization occurs in separate derived files. Source-row identifiers preserve reversibility to the source. A source that is verified but not technically retrievable is marked deferred; numbers are never reconstructed from article prose or figures.

## PC-DB handling

The peer-reviewed ACS Nano PC-DB publication reports 817 nanoparticle formulations and 2,497 proteins from 83 studies. The public GitHub repository currently exposes a directly analyzable human subset of 597 experiments and 1,513 protein columns. The frozen benchmark uses the accessible subset, not the larger publication-level count.

The wide matrix contains 597 × 1,513 experiment-protein cells. A value >0 is interpreted exactly as the repository's Protein_Frequency.R script interprets it: detected/present. Source quantitative values are retained in the raw matrix. A sparse normalized master index is used to avoid duplicating >900,000 cells in a second derived CSV.

## Protein identifiers and structures

Protein names are retained verbatim from source columns. A 2,497-accession file exists in the same repository, but it is not a safe one-to-one positional crosswalk to the 1,513-column human matrix. Accessions are therefore unresolved rather than guessed.

For LEVEL 3, only experimentally mapped PDB structures are stored as experimental structures. For LEVEL 1, PDB/AlphaFold enrichment is a separate downstream metadata-enrichment step; no predicted structure is mislabeled as experimental.

## Material mapping

Mapping is multi-label because one exposed surface can simultaneously be charged, hydrogen-bonding, hydrophobic, aromatic, etc. Rules use exposed surface modification/charge first and core chemistry only where the core is plausibly exposed. Uncertain cases remain ambiguous/unmapped. The exact rule string is stored per material/modification combination.

## LEVEL 2 semantics

The benchmark never collapses directionality:
- increased proteolytic accessibility -> newly exposed / conformational response;
- decreased proteolytic accessibility -> protected / contact-compatible candidate;
- direct covalent-capture or NMR localization -> direct localization/contact evidence.

Protection is not automatically equated with direct contact when stabilization can produce the same experimental signature.

## LEVEL 3 quality hierarchy

A: residue/peptide localization by direct NMR/CSP, crosslink/covalent capture MS, site-directed fluorescence, or ssNMR with folded-state interpretability.
B: broader domain/region evidence that is experimentally localized but lower resolution.
Predicted structures and inferred contacts are never upgraded to A.

## Freeze and anti-cherry-picking

The dated branch was created before InterfaceScout batch scoring. Existing InterfaceScout case scores were not read or used as eligibility criteria. Existing benchmark candidate citations were re-checked against original literature and kept only when they satisfy the frozen source-based rules. PRS, harmonic ANM, and future structure-based/breakable CG outputs are outside this benchmark-construction logic.

## Deferred-source policy

A verified source may remain deferred if raw tables are inaccessible, binary-only and not parsable in the current ingestion path, or if provenance overlap cannot yet be resolved. Deferred is not equivalent to excluded. Inclusion requires a future versioned release; it cannot be retroactively added only because InterfaceScout performs well on it.
