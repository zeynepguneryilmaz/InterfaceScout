# InterfaceScout external-validation protocol (publication freeze)

This protocol is fixed before inspecting benchmark scores.

## Ground-truth inclusion
A benchmark case must use a folded protein, a real solid/material or nanoparticle interface, and experimental protein-side localization at residue, peptide/segment, domain, or molecular-face level. Adsorption amount, DLS, CD/FTIR, activity-only inference, simulation-only contact maps, and systems dominated by major adsorption-induced unfolding are excluded from primary ground truth.

The analyzed coordinate model must represent the experimentally relevant oligomeric state. If the PDB asymmetric unit is not the relevant biological assembly, the relevant assembly must be reconstructed before scoring and the mapping documented.

## Evidence tiers
- Tier A: exact residue or residue-anchor evidence.
- Tier B: experimentally localized peptide/segment/loop/region.
- Tier C: experimentally supported domain or molecular-face orientation.

## Common prediction geometry
All model and baseline variants use the same exposed-residue set, the same 8 A non-transitive same-face coarse-patch geometry, the same hotspot redundancy rule, and the same number of reported candidate centers. Only the center-ranking signal changes.

## Baselines
- B0: geometry-matched exposed-surface null.
- B1: scRSA-only center ranking.
- B2: chemistry-membership-only center ranking.
- B3: residue-local compatibility L = I x scRSA x f_state without spatial persistence.
- B4: single-radius density ranking using a frozen comparator radius.
- B5: full InterfaceScout 5/8 A multiscale-persistence ranking.

## Primary background
Enrichment denominators use all exposed surface residues in the prepared structure, not only chemistry-eligible residues.

## Primary metrics
For Tier A cases:
- exact overlap at top k;
- exposed-surface enrichment at top k;
- minimum Euclidean C-alpha distance from each experimental anchor to a predicted patch;
- sensitivity analysis at 5 A and 8 A proximity thresholds.

For Tier B cases:
- exact patch/segment overlap recall and precision;
- minimum C-alpha distance from the experimental region to a predicted patch;
- near-region recall at 5 A and 8 A.

For Tier C cases:
- molecular-face agreement is scored only when the experimental source defines a directional or domain-level orientation that can be mapped to the structural model. The rule and mapping are documented case by case before scores are generated.

## Model-freeze rule
Benchmark labels must not be used to alter scRSA = 0.05, the 5/8 A pair, the minimum operator, chemistry-channel definitions, the same-face gate, hotspot suppression, or Pareto descriptors. Any post-benchmark change constitutes a new model version and requires a new held-out evaluation.
