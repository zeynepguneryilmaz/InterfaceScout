# InterfaceScout external-evaluation protocol

This protocol separates parameter selection from experimental evaluation and prevents experimental interface labels from entering prediction.

## Ground-truth inclusion
A benchmark case must use a folded protein, a real solid/material or nanoparticle interface, and experimental protein-side localization at residue, peptide/segment, domain, or molecular-face level. Adsorption amount alone, DLS alone, CD/FTIR alone, activity-only inference, simulation-only contact maps, and systems dominated by major adsorption-induced unfolding are not treated as primary localization ground truth.

The analyzed coordinate model must represent the experimentally relevant oligomeric state as closely as possible. Structural mismatches or unresolved experimentally implicated regions must be recorded rather than silently treated as prediction failures or successes.

## Evidence resolution
- **Tier A:** exact residue or residue-anchor evidence.
- **Tier B:** experimentally localized peptide, segment, loop, helix, or broader region.
- **Tier C:** experimentally supported domain or molecular-face orientation.

Comparisons must match the resolution of the experimental evidence. Exact-residue metrics are not used to overinterpret region- or orientation-level evidence.

## Parameter-selection separation
The external benchmark is not used to choose model parameters. Before external experimental evaluation, the numerical parameters requiring empirical selection were examined on a separate adsorption-label-free structural panel: PDB **1CRN, 1R69, 1SHG, 2PPN, 4AKE, 1OMP, 1TIM, and 5CSC**.

- Shrake–Rupley sampling density: **200 points/atom**, selected against a 1000-points/atom numerical reference.
- Multiscale aggregation pair: **6/9 Å**, selected by equal consensus comparison across predefined candidate radius pairs without experimental adsorption labels.

The solvent probe radius is **1.40 Å**, the operational exposure threshold is **scRSA ≥ 0.05**, and coarse patches use a separate **8 Å non-transitive same-face neighborhood**.

## Prediction-first experimental evaluation
Prediction input is restricted to the protein/structure identifier, analyzed chain(s), surface-chemistry mapping required for the experimental system, and solution pH. Experimental localization labels, evidence tiers, DOI metadata, and comparison outcomes are not passed to the prediction function.

Predictions are saved before the comparison stage loads experimental ground truth. The candidate-panel manifest and the ground-truth manifest are retained separately so that this ordering is auditable.

## Evaluation
For residue- or anchor-resolved evidence, evaluation can report direct overlap, first-hit rank, and prespecified Cα spatial proximity (for example ≤5 Å and ≤8 Å). For peptide/segment/region evidence, interpretation is performed at regional resolution. Broad regions are not presented as equivalent to exact-residue recovery. Orientation/domain evidence is kept separate unless an explicit mapping rule is available.

InterfaceScout reports multiple candidate patches because a protein surface can contain more than one chemically plausible contact face. Results therefore should not be collapsed into a single universal "accuracy" value across heterogeneous evidence types.

## No post-evaluation retuning
External experimental labels must not be used to alter the publication settings or chemistry definitions. Any post-evaluation change to the prediction-defining model requires a new held-out evaluation.
