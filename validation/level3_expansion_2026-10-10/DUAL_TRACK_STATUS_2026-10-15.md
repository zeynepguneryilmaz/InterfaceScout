# Experimental Level-3 and independent computational Level-C — audit status

No frozen model/threshold/weight, benchmark manifest, 285 panel outputs, prospective ground truth, or score artifacts changed.

## Experimental Level-3
- Prospective independent scoring unchanged: L3E001 GB3/citrate AuNP and L3E025 apo-transferrin/iron-oxide NP. Published score file on `level3-expansion-score-2026-10-15` has 2 conditions, 2 top3 near8 hits, 1 top5 direct hit.
- Original benchmark is 7 primary + 2 structural-response secondary cases; do **not** pool these with prospective results without a clearly defined stratum. HSA/PAA-Fe3O4 and HCAII silica cases are present in the old frozen manifest.
- HCAII source DOI 10.1529/biophysj.104.054809, open text PMC1305500: engineered pseudo-wild-type C206S with additional site-specific cysteine/AEDANS labels, not exact WT 2CBA construct. Original 2CBA is a structural surrogate, so site-specific comparison must be described as anchor-level and construct qualification highlighted. pH 6.3 phosphate and pH 8.3 Tris adsorption conditions are explicitly given, with 15nm silica.
- HSA source original adsorption step in unbuffered deionized water versus MES pH 6.2 chemical crosslinking stage remains methodological caveat in frozen old case. No retrospective rescore.

## Computational Level-C
- New prospect registry `computational_concordance_candidate_registry.csv` has C001 HEWL/silica, C002 three proteins/silica, C003 HEWL/SWCNT.
- C003 previously transcribed candidate MD contact sites from 2025 literature: 61,62,63,68,69,70,73,103,106,107,111,112,113,116. This is an *unverified candidate extraction* pending direct article figure/methods check, including protein PDB and exact CNT/chirality.
- C001/C002 still lack exact source contact residue lists and mappings.
- None yet passes pre-score inclusion; C003 SWCNT external surface does not presently have verified frozen canonical map. Do not score it by silently borrowing hydrophobic/pi mapping from pristine C60 or polystyrene.
- Keep Level-C visibly distinct from experimental Level-3, especially any MD results reused from the same silica experimental research lineage.

## Next admissible decision
Acquire open primary MD contact table and exact condition mappings for C001 or C002 and independently verify canonical chemistry and PDB; only then lock a source-ground-truth snapshot on a new freeze branch and score on another branch. For experimental Level-3 prioritize source-pH-confirmed residue-defined and folded-state novel family. No change of status is justified by current evidence alone.
