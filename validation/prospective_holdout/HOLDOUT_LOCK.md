# Prospective holdout lock

Lock date: 2026-09-05
Frozen model: `IS-v1-primary-score-freeze`

This directory is intentionally empty of benchmark outcomes at freeze time.

Any protein/material system added here after the freeze is prospective validation evidence only if its experimental/MD/CG ground truth was not used to select or modify the IS v1 formula, cutoffs, chemistry ontology, pH/state rules, or auxiliary-feature policy.

Rules:
1. Register PDB/structure, material chemistry, source conditions and ground-truth provenance before calculating IS performance whenever practicable.
2. Do not change IS v1 after viewing holdout performance.
3. Report all registered holdout systems, including failures.
4. Repeated conditions on one protein are clustered at the protein level.
5. If a model rule is changed after holdout inspection, increment the model version and establish a new untouched holdout set.

At creation of this file there are zero registered prospective holdout proteins.
