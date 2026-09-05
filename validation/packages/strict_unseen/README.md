# InterfaceScout external validation package — Strict unseen prospective holdout

This package is reserved for a true prospective holdout analysis after the InterfaceScout chemistry ontology, pKa source, preprocessing rules and scoring implementation are frozen.

## Critical provenance rule
Any protein or condition whose literature ground truth, contact residues, orientation outcome or validation result has already been inspected during model development or benchmark design is excluded from this strict holdout package, even if that protein was not used to fit numerical parameters.

Because the current 15-protein literature panel has already been inspected to varying degrees while designing the benchmark, it belongs to the All-15 external validation package and must not be relabeled as strict blind validation.

## Current status
No protein is yet admitted to the strict prospective holdout set.

## Admission criteria after freeze
- previously uninspected protein-surface condition;
- exact PDB or defensible sequence-to-structure mapping;
- published residue-level, spatial, orientation, or continuous-contact ground truth defined before InterfaceScout scoring is examined;
- surface chemistry maps unambiguously to the frozen InterfaceScout ontology;
- no benchmark-specific parameter or chemistry-map tuning;
- preprocessing rules fixed in advance, including MODEL selection and biologically required chains.

## Statistical policy
The future holdout will use the same frozen scoring implementation and metrics as the All-15 package. Results will be reported without retuning, including failures, with per-condition metrics plus protein-clustered summary confidence intervals when sample size permits.
