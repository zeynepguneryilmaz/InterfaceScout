# InterfaceScout external validation package — All-15

This package is the complete literature-validation set. It contains 15 unique proteins and all accepted protein–surface conditions. Conditions may be Tier 1 (residue-level quantitative), Tier 2 (orientation/spatial), or stress-test/limitation conditions.

## Statistical policy
- Report the number of unique proteins separately from the number of protein–surface conditions.
- Compute per-condition AUROC/AP/Recall@k/spatial recovery only when residue-level ground truth is sufficiently explicit.
- Use orientation/spatial metrics for Tier 2 conditions.
- Cluster summary uncertainty by protein.
- Keep development-seen proteins in the full package, but label them explicitly.
- Do not tune InterfaceScout parameters to improve this package.

## Development-seen proteins
- BSA
- Lysozyme

All other proteins are assigned to the strict-unseen sensitivity set unless later provenance review shows they were inspected during model development.
