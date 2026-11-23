# Biofilm protein representation in ISV (Dataset 2)

The ISV initial-scientific-validation candidate pool now includes **12 biofilm-related PDB entries**, alongside 68 earlier protein candidates. [All 80 proposed entries](ISV_protein_candidate_pool_80_v2.csv) are a **screening pool**, not 80 accepted or structurally qualified proteins.

**RCSB source-record screen:** [GitHub Actions run 45268134181](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45268134181), **12/12 metadata checks PASS**. All 12 entries have experimentally determined protein-bearing structures: **10 X-ray structures and 2 cryo-EM structures**. These checks establish source existence, experimental method, basic deposited polymer description and reported resolution—not coordinate completeness, full-protein representation, family independence or InterfaceScout applicability.

## Biofilm panel by biological mechanism

| Subgroup | PDB | Protein | Resolution (Å) | Scientific use | Admission priority |
| --- | --- | --- | ---: | --- | --- |
| Initial cell adhesion | [1R17](https://www.rcsb.org/structure/1R17) | SdrG fibrinogen-binding adhesin | 1.86 | Adhesin-side surface mapping; receptor peptide is a separate entity | High |
| Biofilm surface layer | [4BHU](https://www.rcsb.org/structure/4BHU) | BslA, bacterial hydrophobin | 1.91 | Air/water and biofilm outer surface chemistry; oligomer caution | High |
| EPS-associated protein | [5OF2](https://www.rcsb.org/structure/5OF2) | TasA soluble conformation | 1.86 | Monomeric folded-state chemistry; not a substitute for amyloid fibers | High |
| Matrix cell-cell linkage | [4KKP](https://www.rcsb.org/structure/4KKP) | RbmA | 2.50 | Biofilm matrix adhesion/multimer contexts | High |
| Abiotic surface attachment | [6MLT](https://www.rcsb.org/structure/6MLT) | Bap1 | 1.90 | V. cholerae adhesion-relevant extracellular matrix protein | High |
| Lectin-polysaccharide binding | [1GZT](https://www.rcsb.org/structure/1GZT) | LecB | 1.30 | P. aeruginosa lectin; calcium and sugar-binding context | High |
| EPS synthesis/turnover | [5TCB](https://www.rcsb.org/structure/5TCB) | PelA hydrolase domain | 1.535 | Extracellular polysaccharide processing; not cell-surface adhesin | Contextual |
| EPS synthesis regulation | [4DN0](https://www.rcsb.org/structure/4DN0) | PelD cytosolic receptor domain | 2.30 | Intracellular c-di-GMP recognition, not external interface contact | Contextual |
| EPS degradation | [5BXA](https://www.rcsb.org/structure/5BXA) | PslG | 1.90 | Polysaccharide hydrolase with mannose in deposited structure | Contextual |
| Curli secretion | [6LQH](https://www.rcsb.org/structure/6LQH) | CsgG–CsgF complex | 2.94 cryo-EM | Membrane/oligomer-bound secretion; requires separately specified context | Special challenge |
| Amyloid biofilm fiber | [8C50](https://www.rcsb.org/structure/8C50) | Curli-like CsgA from *Pontibacter korlensis* | 7.60 cryo-EM | Very-low-resolution, non-globular fibril; NOT routine residue-level validation | Special challenge |
| EPS turnover | [9OZL](https://www.rcsb.org/structure/9OZL) | RbmB | 2.10 | VPS polysaccharide lyase; not direct substrate adsorption site | Contextual |

**Scientific inclusion rule:** The target of ISV is **protein-side physicochemical surface characterization**. Some biofilm-related proteins are directly relevant to abiotic-surface adhesion (SdrG, Bap1, BslA); others are components of matrix assembly, enzymes, intracellular regulators or membrane fibrils. These categories must remain separate. A known biofilm role does **not** constitute experimental proof of where a protein binds an abiotic material.

**Quality gates:** verify the correct deposited entity/chain (e.g., exclude the fibrinogen peptide when testing SdrG alone), biological assembly, residues observed in coordinates, missing non-backbone heavy atoms, construct species, ligand/cofactor handling, domain redundancy, 11-map viability and numerical finite outputs. For amyloid/membrane examples, validate applicability and either put in a separate structural stress-test stratum or exclude from the core protein analysis.

**Dataset policy:** Biofilm entries are part of **Dataset 2 (ISV)**, not a third dataset. The final goal is **at least 50 accepted protein structures** overall. A limited, stratified subset of actual biofilm-surface/matrix representatives is preferable to inflating the group with proteins that only indirectly regulate biofilm production.

**Existing limitation:** The separate 68-structure model-1 Cα coordinate coverage workflow encountered 62 external coordinate-file HTTP 404 errors and 6 provisional positive results. This was a **download/access failure, not evidence of bad structures**; it must be repaired and rerun before coordinate completeness is asserted for the pool. The successful 68-entry RCSB metadata screen remains valid as a metadata-only audit.
