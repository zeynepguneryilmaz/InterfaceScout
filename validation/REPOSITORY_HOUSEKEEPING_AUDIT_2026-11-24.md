# InterfaceScout repository housekeeping and scientific-claim cross-check

Date: 2026-11-24. Scope: connected GitHub file-by-file inspection of the primary README, validation index, authoritative dataset file, supported software requirements, workflow configuration, the optional external tools, and the citation file. **This is not a complete source-tree or end-to-end scientific verification.**

## Retained authoritative public content

- Public application: [README](../README.md); version 1.1.0 and the 11 separate protein-side physicochemical mapping channels.
- Active validation datasets: [TV48 and ISV76](ACTIVE_DATASETS.md). TV48 is a 48-structure technical verification; ISV76 is a frozen *candidate structural panel* containing 76 structures, not 76 independently accepted experimental truths.
- ISV selection: [fixed manifest](ISV_active_candidates_76_v3.csv) and [decision record](ISV_SELECTION_DECISION_76.md); exclude 3GHG, 1QJS, 4DN0, 6LQH; retain 1STP with biotin-state caveat.
- Previous 80-structure source audit, early literature conditions, QA measurements and fixed result archives remain historical provenance and must not be relabeled as ISV76 experimental validation.
- TV48 production results and explicit structural/algorithmic limits remain [documented](FINAL_V110_RESULTS.md).
- Optional real APBS/PDB2PQR/PROPKA is not one of the eleven canonical channel definitions and is not adsorption free energy.

## Housekeeping performed

1. Removed automatic GitHub workflows for discontinued SURFMAP native full-panel, targeted PEP-Patch rescue and 1GCT atom repair.
2. Removed their two newly introduced unvalidated experiment-only Python helpers. Earlier actual benchmark tools and historical logs were *not* deleted.
3. Retained external P2Rank/PEP-Patch/preflight workflows as **manual opt-in only** for reproducibility; no further push-triggered full-panel benchmarking.
4. Retained the Protein-sol 76-PDB full-panel as **manual opt-in only** to avoid repeat submissions to an external academic server on minor edits. The limited 2-PDB diagnostic remains separately available.
5. Corrected contradictory references to an ISV panel 'still to be selected' in both [validation index](README.md) and [dataset index](ACTIVE_DATASETS.md).
6. Corrected historical 80-structure QA strata: numbers from an 80-structure screen cannot serve as subgroup totals for the current 76.
7. Added [local automated repository consistency checks](repository_consistency_check.py) to the cross-platform publication CI; these inspect structural panel counts/exclusions, documentation scope and retired workflow configuration.

## Outstanding verification gates — not claimed complete

- Inspect the latest cross-platform CI run after these edits; committing a file is not passing CI.
- Inspect real FreeSASA and Protein-sol run artifacts and explicitly record numerical counts and parser/provider failures.
- Repeat production-integrated APBS/SAS sampling regression on current head.
- Audit original experimental source figures, residue numbering, construct matching and structural modifications for ISV scientific admission.
- Perform full recursive inventory (broken links, duplicates, unused internal modules, licenses, files sizes, frontend-to-backend API compatibility). The current GitHub connector supports fetching exact paths and searching but did not expose a complete recursively enumerated working tree in this review. No blanket statement that every repository file has been checked should be made.

## Publication-safe interpretation

A successful 11-map computation on TV48 or ISV76 demonstrates software executability and numerical consistency, not protein–material adsorption accuracy. PEP-Patch is a physicochemical-patch comparator, Protein-sol is a surface-descriptor comparator, FreeSASA checks solvent-accessibility computations, while P2Rank and SURFMAP are not required for the current manuscript. No unobserved run, guessed timing, or inferred result is recorded as a successful execution.
