# Archived 282-protein patch audit — 2026-10-16

## Scope and provenance
All 12 original frozen Github Actions artifacts (run 38614780417) were downloaded and read; 282 unique protein JSON files and 3102 chemistry-channel rank-1 patch records were extracted. Frozen V1 was not rerun or changed.

## Initial AF model-confidence screen
- 282 run successes
- 239 protein-level STANDARD by original frozen QC
- 97 provisional tighter candidates: median pLDDT >=90; fraction pLDDT <50 <=0.05; AF reported coverage between 0.95 and 1.05
- No structure-level filter proves every predicted region is locally reliable or belongs to a mature chain.

## Terminal-bias audit using actual JSON rank-1 patch members
- All 282: 660/3102 (21.28%) centers within precursor positions 1-30.
- STANDARD 239: 612/2629 (23.28%) centers within precursor positions 1-30.
- Strict provisional 97: 320/1067 (29.99%) centers within positions 1-30; 319/1067 rank-1 patches have *at least one* member within positions 1-24.
- Under a crude reference assigning each protein's 11 patch centers uniformly over its **reported surface residues**, 97-candidate expected count first 30 is 186.52. Observed/expected 1.716. This is descriptive, NOT spatial randomized significance.
- All-channel terminal-near fraction among provisional 97 (within first or last 30 residues): 520/1067 (48.73%).
- Human albumin LP025 (precursor P02768): 6/11 rank-1 patches overlap residues 1-24. Albumin maturation removes signal/propeptide positions 1-24; the frozen protein-input run cannot be equated with mature circulating HSA without preparation-aware reevaluation.

## Interpretation and required gate
Current 97 is a **provisional structure-global-confidence stratum**, not a biologically/locally validated main set. Some terminal sites are legitimate; some correspond to signal peptides, propeptides or mismatched isoform precursors. Do not apply a universal first-30 deletion or hide problem cases. Do not change V1 ranking based on output quality. Use sequence-specific UniProt features and residue-wise AF-pLDDT, and optionally PAE/domain arrangement and assembly state. The script `scripts/validation/patch_residue_confidence_qc.py` performs source-specific, post-hoc rank-1 inspection when full JSON outputs and source structures can be downloaded. A local extracted QC CSV with per-protein/channel positions exists as separate deliverable. Residue pLDDT and UniProt processing have not yet been executed for all 97 due to unavailable external PDB/UniProt API access in the execution container.

## Data integrity warning
Filtering to 97 proteins after looking at 282 results is post-hoc structural QC. Original source denominator 285 and 282 completed run records remain archived. Differences in chemistry ranking after filtering should be described as sensitivity, not independent validation.
