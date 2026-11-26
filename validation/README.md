# Validation and reproducibility

InterfaceScout v1.1.0 uses two principal structural panels.

| Panel | Purpose | Records |
| --- | --- | --- |
| TV48 | Technical verification of the eleven maps, accessibility, region construction, ranking and numerical robustness | [48-structure manifest](strict_apo_48_selected_candidates_v2_final.csv), [production verification](FINAL_V110_RESULTS.md), [extended verification](EXTENDED_V110_SCIENTIFIC_RESULTS.md) |
| ISV76 | Application and structural applicability of protein-surface mapping; not 76 accepted independent experimental cases | [76-structure manifest](ISV_active_candidates_76_v3.csv), [selection record](ISV_SELECTION_DECISION_76.md), [application results](ISV76_FINAL_MODEL_INITIAL_RESULTS.md) |

**Download:** [Research datasets](../datasets/README.md) contains persistent archives of the original calculation outputs and their checksums. [Publication evidence](PUBLICATION_EVIDENCE.md) maps analyses to source runs. [Active datasets](ACTIVE_DATASETS.md) records selection history and evidence criteria.

## Numerical model and methods

- [Dual accessibility](DUAL_ACCESSIBILITY_LITERATURE_AND_METHOD.md)
- [Regional Chemistry Profile](REGIONAL_CHEMISTRY_EVIDENCE_AND_METHOD.md)
- [Parameter selection](parameter_selection.json)
- [FreeSASA comparison](ISV76_FREESASA_76_INDEPENDENT_COMPARISON_FINAL.md)
- [Continuum electrostatics protocol](TV48_APBS_SAS_PROTOCOL.md)

## Historical development and exploratory evidence

Earlier structure screens, parameter-development panels and retrospective literature comparisons remain in this directory for auditability. The historical 80-structure screen is superseded by ISV76. The frozen nine-condition literature comparison is documented in [benchmark_manifest.json](benchmark_manifest.json) and [benchmark_case_results.csv](benchmark_case_results.csv); the check `python validation/verify_publication.py` reproduces its historical aggregate checks and parameter-selection assertions. These records are distinct from the publication application panel.

The analysis scripts and original workflow definitions are retained so source runs can be reproduced using the commit recorded in each data archive.
