# 282 successful proteins — evidence gap-closure report (16 Oct 2026)

## Completed and independently checkable
- Exactly 282 frozen V1 SUCCESS records carried from the 285 frozen source panel; 3 source proteins with no AF model excluded from successful-run denominator only.
- 282-row source/status/QC/evidence crosswalk: `PANEL_282_EVIDENCE_CROSSWALK_2026-10-16.csv`.
- 282-row targeted literature screen: `PANEL_282_TARGETED_LITERATURE_SCREEN_2026-10-16.csv`; 10 specific protein entries have newly indexed literature leads. **This was a targeted, non-exhaustive literature search, not 282 comprehensive systematic reviews.** A blank citation does not mean literature does not exist.
- Literature leads include 10.3390/ijms22168722 (MUS-Au MD hemoglobin, albumin, alpha1-antitrypsin, C3); 10.1039/D5EN00784D (2026 UA PS nanoparticle modelling of albumin, fibrinogen, immunoglobulin, ApoA-I, C3); 10.1021/acs.molpharmaceut.4c01533 (lipid bilayers albumin/ApoA-I), 10.1016/j.jmb.2022.167771 (HSA/C3b MD), 10.1021/ac201889j (HSA/PAA Fe3O4 legacy experimental), 10.1021/acs.analchem.9b00482 (apo-transferrin/IONP prospective).
- Crucial structure-mapping caveats: fibrinogen literature models heterotrimer but panel rows are alpha, beta and gamma separate protein chains; hemoglobin MD commonly treats tetramer, panel has separately modelled alpha and beta chains; C3b fragment ≠ full C3; precursor ALBU/TRFE can differ from mature crystallographic protein numbering.
- None of these *newly identified* literature leads was converted to verified exact residue-contact labels for the specific panel AF model and frozen chemistry channel. The new Level-C scoreable count is therefore zero.

## Frozen findings still supported
- V1 run success 282/285 = 98.95%, 3102 channel landscapes.
- Original experimental primary benchmark top3 near8 hits 7/7, top5 direct hits 4/7, 2 secondary structural-response challenge misses.
- Independently frozen prospective Level3: GB3/Au K13 strong rank2 overlap; apoTf/IONP only top3 near8 recall 0.075, top5 near8 recall 0.10, no top5 direct hit.
- Existing experimental crosswalk of panel proteins is not an independent new revalidation against AF V1; experimental PDB versus AF precursor mapping must be checked first.

## Critical non-result
**The requested 282-protein independent experimental/MD residue-concordance analysis is not completed**, because independently published site-level reference annotations with correct matched species, construct, material-facing chemistry, pH and sequence numbering have not been acquired for each protein. Reporting 282 correspondence scores would be fabricated. The honest status is 282 verified V1 run entries and targeted source triage, not 282 independent localization validations.

## Frozen integrity
No V1 source, weights, thresholds, channels, original benchmark, prospective ground truth, or experimental scores changed. Level-C MD simulations are literature concordance only, not experimental Level-3 accuracy.
