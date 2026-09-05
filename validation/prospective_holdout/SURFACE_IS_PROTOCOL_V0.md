# Surface-IS atomic-surface prototype V0

Status: exploratory comparison frozen before Surface-IS scores are computed. This does not modify IS-v1-Core-primary.

## Scientific question
Does a protein-molecular-surface-first representation improve recovery of experimentally mapped protein/material interfaces relative to the current residue-first scRSA formulation?

## Systems
All three first-pass systems use anionic/carboxyl-rich material chemistry so that the representation change can be isolated without introducing a new material ontology.

- GB3 / 2OED chain A / citrate-AuNP / pH 6.4 / experimental K4,K13,K50.
- beta-2-microglobulin / 1JNJ chain A / citrate-AuNP / pH 7.7 / experimental residues 2,3,26,28,29,30,33,55,56,58,59.
- HSA / 2VUF chain A / PAA-Fe3O4 / pH 7.0 modeling convention / experimental peptide regions 373-389,403-410,414-428.

## Required ablation arms
A. `RAW_RESIDUE_CORE`: canonical residue-first Core on the source-matched chain, retaining the PDB's explicit hydrogens exactly as the current backend does.
B. `HSTRIP_RESIDUE_CORE`: identical canonical residue-first Core after removing explicit H/D atoms before Shrake-Rupley SASA. This isolates the SASA-convention effect.
C. `ATOMIC_SURFACE_CORE`: heavy-atom molecular-surface-first prototype on the same H-stripped structure.

No GNM, RIN, APBS, curvature, planarity, learned coefficients or fitted weights are allowed in V0.

## Atomic-surface representation
- Shrake-Rupley probe radius 1.4 A, 200 points.
- Explicit H/D atoms are removed before SASA calculation.
- Every heavy atom with SASA > 0 contributes a molecular-surface point located at its atomic coordinate with exposed area equal to that atom's SASA.
- Candidate patch centers are all exposed heavy atoms; there is no residue-level scRSA gate.

## Anionic-surface compatible functional groups
The atom ontology is a direct atom-level translation of the frozen residue chemistry ontology, without literature-energy weighting:
- LYS: NZ, weighted by protonated-state availability.
- ARG: NE, NH1, NH2, weighted by protonated-state availability.
- HIS: ND1, NE2, weighted by protonated-state availability.
- SER: OG, weight 1.
- THR: OG1, weight 1.

For compatible atom a:
`q_a = SASA_a * f_state(residue,pH)`.
All other atoms have q_a=0.

## Multiscale surface patch score
For every exposed heavy-atom center k:
`D_k(R) = sum_a [distance(k,a)<=R] q_a`, for R=5 and 8 A.
Each radius is normalized by its within-protein maximum.
`M_k = 100 * min(D5_norm, D8_norm)`.

Residue-level prediction is derived from the surface, not used to construct it:
`SurfaceScore_i = max(M_k)` over exposed heavy atoms belonging to residue i.

This V0 score is a deterministic surface-chemistry patch localization score, not a probability or adsorption free energy.

## Evaluation
For residue-level GT: rank all residues present in the structure by the arm-specific residue score and compute AUROC, AP, Recall@5/10/20, and top-10 spatial recovery at 5/8/10 A.
For HSA peptide-level GT: peptide-region spatial recovery is primary; residue AUROC/AP are descriptive only.
Also export atom-level functional-group SASA for each experimental GT residue so that GB3 K4/K13/K50 can be audited directly.

## Interpretation rule
Improvement in arm B versus A is attributed to SASA hydrogen-convention correction, not to surface-first modeling. Additional improvement in arm C versus B is the evidence attributable to the atomic-surface-first representation.
No formula, atom ontology, radii or thresholds may be changed after V0 results are observed.
