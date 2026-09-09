# InterfaceScout external-validation protocol

This protocol is fixed before inspecting external experimental benchmark scores.

## Ground-truth inclusion
A benchmark case must use a folded protein, a real solid/material or nanoparticle interface, and experimental protein-side localization at residue, peptide/segment, domain, or molecular-face level. Adsorption amount, DLS, CD/FTIR, activity-only inference, simulation-only contact maps, and systems dominated by major adsorption-induced unfolding are excluded from primary ground truth.

The analyzed coordinate model must represent the experimentally relevant oligomeric state. If the PDB asymmetric unit is not the relevant biological assembly, the relevant assembly must be reconstructed before scoring and the mapping documented.

## Evidence tiers
- Tier A: exact residue or residue-anchor evidence.
- Tier B: experimentally localized peptide/segment/loop/region.
- Tier C: experimentally supported domain or molecular-face orientation.

## Parameter-selection separation
The external benchmark is not used to choose model parameters. Before benchmark evaluation, two parameters were selected on a separate adsorption-label-free structural panel (PDB 1CRN, 1R69, 1SHG, 2PPN, 4AKE, 1OMP, 1TIM, and 5CSC):
- Shrake-Rupley sampling density: 200 points/atom, selected as the lowest density satisfying predefined convergence criteria against a 1000-points/atom numerical reference.
- Multiscale aggregation pair: 6/9 A, selected by equal consensus comparison across predefined local radius pairs, without designating any candidate pair as the reference.

Other settings are not benchmark-optimized: the solvent probe radius is 1.40 A, the operational scRSA exposure threshold is 0.05, generic side-chain reference pKa values are used for state availability, and the same-face gate is parameter-free (positive dot product between centroid-to-C-alpha outward vectors).

## Common prediction geometry
All model and baseline variants use the same exposed-residue set, the same 8 A non-transitive same-face coarse-patch geometry, the same hotspot redundancy rule, and the same number of reported candidate centers. The 8 A coarse-patch neighbourhood is a separate structural rule and is not equated with the outer 9 A multiscale aggregation radius.

## Baselines
- B0: geometry-matched exposed-surface null.
- B1: scRSA-only center ranking.
- B2: chemistry-membership-only center ranking.
- B3: residue-local compatibility L = I x scRSA x f_state without spatial persistence.
- B4: single-radius 9 A density ranking.
- B5: full InterfaceScout 6/9 A multiscale-persistence ranking.

## Primary background
Enrichment denominators use all exposed surface residues in the prepared structure, not only chemistry-eligible residues.

## Primary metrics
For Tier A cases:
- exact overlap at top k;
- exposed-surface enrichment at top k;
- minimum Euclidean C-alpha distance from each experimental anchor to a predicted patch;
- proximity recall at prespecified coarse thresholds.

For Tier B cases:
- exact patch/segment overlap recall and precision;
- minimum C-alpha distance from the experimental region to a predicted patch;
- near-region recall at prespecified coarse thresholds.

For Tier C cases:
- molecular-face agreement is scored only when the experimental source defines a directional or domain-level orientation that can be mapped to the structural model. The rule and mapping are documented case by case before scores are generated.

## No post-benchmark retuning
External benchmark labels must not be used to alter scRSA = 0.05, SASA sampling = 200 points/atom, the 6/9 A multiscale pair, the minimum persistence operator, chemistry-channel definitions, the same-face gate, the 8 A coarse-patch geometry, hotspot suppression, or Pareto descriptors. Any post-benchmark change constitutes a new model version and requires a new held-out evaluation.
