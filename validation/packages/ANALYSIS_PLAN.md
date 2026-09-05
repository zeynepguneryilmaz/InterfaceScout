# InterfaceScout literature-validation analysis plan

## Packages
1. `all15`: full validation package, 15 unique proteins.
2. `strict_unseen`: sensitivity/generalization package, excluding development-seen BSA and lysozyme.

## Per-condition quantitative metrics (Tier 1)
- AUROC
- Average precision (AP)
- Recall@5, Recall@10, Recall@20
- Exact anchor rank and percentile
- Spatial recovery at 5, 8, and 10 Å using top-10 predicted patch centers
- Nearest predicted patch distance
- Exposure-matched permutation p-value (10,000 permutations)
- Spearman correlation / NDCG when continuous residue-contact frequencies are available

## Tier 2 metrics
- Orientation/face recovery
- Distance from predicted patch to reported binding face/domain
- Active-site accessibility agreement where experimentally assessed
- Supported-mode recovery for multimodal/orientation-specific literature

## Aggregation/statistics
- Report condition count and unique-protein count separately.
- Macro-average/median across conditions.
- Protein-clustered bootstrap 95% confidence intervals.
- Leave-one-protein-out sensitivity analysis.
- Report fraction of Tier-1 conditions with exposure-matched permutation p<0.05.
- Do not pool repeated conditions from the same protein as statistically independent replicates.
- Do not tune InterfaceScout parameters based on benchmark outcomes.

## Preprocessing rules
- Use exact literature PDB when specified.
- For NMR ensembles, use MODEL 1 unless the source explicitly validates an ensemble representation.
- Preserve biologically required chains/assemblies (e.g. alpha-chymotrypsin A/B/C; fibrinogen multichain assembly).
- Never collapse orientation-specific contact sets into a single exhaustive truth set for the primary metric unless the source explicitly defines the union as the adsorption ensemble.
- Material chemistry mapping must be frozen before evaluating blocked chemistry conditions (H-bond donor/acceptor and HAp cases).

## Reporting policy
- Keep failures and weak conditions.
- Distinguish `complete`, `ready`, `pending`, `blocked`, and `stress` conditions.
- The full-package headline can summarize All-15 performance; strict-unseen performance must be reported as a sensitivity/generalization analysis.
