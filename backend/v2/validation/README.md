# Validation and parameter-selection source

This directory contains the scripts and frozen definitions used for the manuscript's parameter-selection and experimental comparison analyses. Generated intermediate files are intentionally not stored in the repository.

## Development panel

`select_core_parameters_final.py` defines the eight-protein adsorption-label-free development panel, and `select_core_parameters_consensus.py` performs the numerical and spatial comparison used to select the final settings.

Development structures: 1CRN, 1R69, 1SHG, 2PPN, 4AKE, 1OMP, 1TIM, and 5CSC.

## Experimental comparison

`benchmark_experimental_final.json` contains the literature-derived experimental annotations and source provenance used for the manuscript comparison. Predictions are generated before those annotations are loaded for comparison.

The associated scripts are:

- `generate_experimental_predictions.py` — generates InterfaceScout predictions for the literature cases.
- `compare_experimental_predictions_v2.py` — compares saved predictions with the literature-supported protein-side annotations.
- `summarize_experimental_validation.py` — produces the compact case-level summary used in the manuscript.
- `reproduce_publication.py` — optionally regenerates the complete local analysis package for audit or rechecking.

Local generated files are ignored by git. The public `examples/` directory instead contains only the normal InterfaceScout user-facing outputs (Excel workbook and B-factor PDB) for the PDB structures used in the manuscript evaluation.

See `PROTOCOL.md` for the frozen comparison rules and interpretation.
