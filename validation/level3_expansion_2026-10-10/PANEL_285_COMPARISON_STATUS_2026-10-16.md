# InterfaceScout — full 285-panel evidence comparison status (2026-10-16)

## Frozen computational run: COMPLETE; per-protein literature validation: INCOMPLETE

- Frozen source size: 285 (PC-DB identification frequency panel, not residue-contact gold standard).
- Frozen V1 successful: 282 (11 channels each = 3102 landscapes); three without usable AlphaFold models; no runtime errors.
- The present crosswalk reconciles **every one of the 285 identifiers** to a V1 status, UniProt and AF model where available, structure coverage, QC, surface-residue count and chemistry-specific output statistics.
- Files: `PANEL_285_V1_AND_EVIDENCE_AUDIT_2026-10-16.csv` (285 rows) and earlier `PANEL_285_EVIDENCE_CROSSWALK_INITIAL.csv`.
- **NOT** a completed independent literature search over 285 protein × material pairs; unreviewed entries remain transparently marked. No unreferenced accuracy/concordance statistic is calculated.

## Source-side verified overlaps
- LP025 human serum albumin: legacy HSA–PAA/Fe3O4 interface peptide mapping exists (DOI 10.1021/ac201889j). The legacy pH=6.2 is *capture-stage* pH, not documented pre-crosslink adsorption pH; existing frozen score not edited. Also appears among proteins in Subbotina et al. 2026 DOI 10.1039/D5EN00784D; paper's adsorption-energy maps cannot automatically be converted into residue-contact gold.
- LP028 serotransferrin: prior admitted prospective human **apo**-transferrin–iron oxide (DOI 10.1021/acs.analchem.9b00482); requires exact apo/panel AF structural-state and numbering mapping for any panel-specific comparison; do not duplicate old prospective scores.
- LP178 human Lysozyme C: poultry HEWL MD contacts are different species and cannot transfer as exact labels.
- LP283 ubiquitin-ribosomal eS31 fusion: not equivalent to free ubiquitin without cleavage mapping.
- LP001 ApoA-I and LP002 C3: included in the 2026 RSC PS nanoparticle *independent computational* paper; evaluate source orientation/contact data and construct match before Level-C, no auto-score.
- Frozen original experimental benchmark contains 7 primary cases; two prospective ADMIT are separately scored. These are not evidence that 285 distinct panel proteins are independently validated.

## 2026 independent computational paper screening
Subbotina, McElligott and Lobaskin. `10.1039/D5EN00784D`, open access, accepted Jan 22 and first published Jan 28 2026. UA/CoronaKMC assesses HSA, fibrinogen, IgG, C3, ApoA-I on 10-nm PS, PS-NH2 and PS-COOH, with residue-granularity bead potentials and orientation/adsorption-energy maps. **No directly extracted, same-protein-mapped contact-residue gold sets** verified here; their model predicts adsorption energies/corona, not the V1 claim. Candidate Level-C only, condition mapping must be separated by outer chemistry and PDB. Strong risk of using source code-derived orientation estimates rather than physical contact labels without pre-score locked definitions.

## Current valid inference
285-wide test establishes *executability/applicability* and chemically variable coarse interface landscapes, not residue-localization experimental accuracy. The numerical pairwise overlap is confounded by protein size and related channels can be identical. Claims/weights/thresholds unchanged.

## Actual next step for completed literature concordance
Systematically retrieve paper/SI or supplementary deposited MD trajectory contact tables for **all 285**, with a zero-reference outcome for any protein lacking source data; record source DOI, protein construct, pH, external contact residues, surface exposure chemistry, and matching V1 output map **before score**. Freeze admissible labels and score on an independent branch. No shortcut converts PC-DB composition/detection frequency into interface-localization ground truth.
