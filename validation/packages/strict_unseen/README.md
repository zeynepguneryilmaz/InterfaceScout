# InterfaceScout external validation package — Strict unseen sensitivity set

This package is the sensitivity/generalization subset of the full 15-protein validation package. It excludes proteins explicitly used during InterfaceScout development or previously inspected in development-facing analyses.

## Excluded from strict-unseen sensitivity statistics
- BSA
- Lysozyme

## Included proteins
- Fibronectin III8-10
- Cytochrome c
- Alpha-chymotrypsin
- Myoglobin
- Hemoglobin
- Trypsin
- RNase A
- Acetylcholinesterase
- WNV E protein domain III
- GB3 / Protein G B1
- Proteinase K
- Carbonic anhydrase II
- Fibrinogen

This yields 13 unique proteins before any later provenance-driven exclusion.

## Statistical policy
- Same frozen InterfaceScout parameters as the All-15 package.
- No benchmark-specific tuning.
- Per-condition metrics are identical to the full package where ground truth permits.
- Summary statistics use condition-level macro summaries with protein-clustered bootstrap confidence intervals and leave-one-protein-out sensitivity analysis.
