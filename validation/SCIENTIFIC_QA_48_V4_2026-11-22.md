# InterfaceScout 48-protein scientific screening — executed v4

**Run (2026-11-22):** https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/44981110375
**Code and data:** validation/final_panel_quality_48.py, validation/final_panel_compat_48.py, validation/development_panel_48_candidates_v4_REVIEW.csv.
**This is screening, not a declaration that all 48 structures are scientifically approved.**

## Actual completed results
- 48/48 RCSB metadata, biological-assembly mmCIF structures and wwPDB validation XML reports retrieved and parsed.
- 48/48 PDB IDs CATH domain annotations found (123 domain assignments).
- 48/48 full biological-assembly inputs produced **11 InterfaceScout maps** without code exceptions under pH 7.4, ionic strength 150 mM, 298 K. The chemistry engine does not model all ligand/metal effects, so this is execution compatibility only.
- CATH **8 homologous-superfamily IDs** shared by **13 PDB records**; do not call 48 independent homologous families.
- 33/48 R-free fields not present in RCSB historical refinement metadata; **missing != failed**. wwPDB DCC fields are exported separately.

## Screening categories
| Category | Number | PDB IDs |
|---|---:|---|
| HIGH_CONCERN | 14 | 2PPN, 1OMP, 1TIM, 5CSC, 1GFL, 1ALD, 1A4Y, 1HVR, 1LDM, 1CLL, 1B0L, 1MUP, 1BEB, 1RWZ |
| REPAIR_REQUIRED | 5 | 1TEN, 1BRS, 1E9I, 1CSP, 1DPS |
| HOMOLOGY_REVIEW | 6 | 4AKE, 1A6M, 1A3N, 1GKY, 1CDG, 1SMD |
| OTHER_REVIEW | 9 | 1BPI, 2CBA, 1J8U, 4M6J, 1PPA, 1BTL, 1C5E, 1PLC, 2TRX |
| SCREENING_OK | 14 | 1CRN, 1R69, 1SHG, 1UBQ, 1PGA, 1HEL, 1HRC, 1CSE, 1TCA, 1STN, 1RNH, 1RBL, 1A1X, 1THW |

The categories are nonoverlapping *priority buckets*. Individual flags can overlap. HIGH_CONCERN uses **review triggers**, not official wwPDB pass/fail thresholds: clashscore >25, Ramachandran outliers >2%, rotamer outliers >10%, RSRZ outliers >10%, or legacy X-ray resolution >2.5 Å. No output-based tuning was performed.

## Highest-priority wwPDB flags
| PDB | Resolution Å | Clashscore | Rama outliers % | Rotamer outliers % | RSRZ outliers % | Flags |
|---|---:|---:|---:|---:|---:|---|
| 2PPN | 0.92 | 6.21 | 0.00 | 0.00 | 47.66 | WWPDB_RSRZ_OUTLIERS_GT10PCT |
| 1OMP | 1.8 | 10.90 | 0.82 | 11.64 | 1.08 | WWPDB_ROTAMER_OUTLIERS_GT10PCT |
| 1TIM | 2.5 | 114.32 | 12.45 | 14.84 | — | WWPDB_CLASHSCORE_GT25; WWPDB_RAMA_OUTLIERS_GT2PCT; WWPDB_ROTAMER_OUTLIERS_GT10PCT |
| 5CSC | 2.8 | 76.62 | 6.26 | 31.40 | — | RESOLUTION_EXCEPTION_OR_FAIL; WWPDB_CLASHSCORE_GT25; WWPDB_RAMA_OUTLIERS_GT2PCT; WWPDB_ROTAMER_OUTLIERS_GT10PCT |
| 1GFL | 1.9 | 9.85 | 0.44 | 11.75 | 1.52 | WWPDB_ROTAMER_OUTLIERS_GT10PCT |
| 1ALD | 2 | 40.75 | 4.43 | 18.56 | — | WWPDB_CLASHSCORE_GT25; WWPDB_RAMA_OUTLIERS_GT2PCT; WWPDB_ROTAMER_OUTLIERS_GT10PCT |
| 1A4Y | 2 | 16.86 | 1.04 | 14.53 | — | WWPDB_ROTAMER_OUTLIERS_GT10PCT |
| 1HVR | 1.8 | 7.80 | 0.52 | 12.96 | — | WWPDB_ROTAMER_OUTLIERS_GT10PCT |
| 1LDM | 2.1 | 34.54 | 0.92 | 15.38 | — | WWPDB_CLASHSCORE_GT25; WWPDB_ROTAMER_OUTLIERS_GT10PCT |
| 1CLL | 1.7 | 29.17 | 0.00 | 17.21 | 2.08 | WWPDB_CLASHSCORE_GT25; WWPDB_ROTAMER_OUTLIERS_GT10PCT |
| 1B0L | 2.2 | 15.94 | 1.16 | 16.06 | — | WWPDB_ROTAMER_OUTLIERS_GT10PCT |
| 1MUP | 2.4 | 51.36 | 0.65 | 19.40 | — | WWPDB_CLASHSCORE_GT25; WWPDB_ROTAMER_OUTLIERS_GT10PCT |
| 1BEB | 1.8 | 15.50 | 0.65 | 10.51 | 7.69 | WWPDB_ROTAMER_OUTLIERS_GT10PCT |
| 1RWZ | 1.8 | 9.32 | 0.00 | 0.95 | 16.39 | WWPDB_RSRZ_OUTLIERS_GT10PCT |

**Important interpretation:** 2PPN (RSRZ 47.66%) and 1RWZ (RSRZ 16.39%) require density-fit inspection. 1TIM (clashscore 114.32; Ramachandran 12.45%) and 5CSC (clashscore 76.62; Ramachandran 6.26%; 2.8 Å) should not be blindly included in high-quality *primary* parameter selection even though they belong to the legacy eight.

## CATH homologous superfamily collisions
| CATH H | PDB entries |
|---|---|
| 3.40.50.300 | 1GKY / 4AKE |
| 3.40.190.10 | 1B0L / 1OMP |
| 3.20.20.70 | 1ALD / 1TIM |
| 1.10.490.10 | 1A3N / 1A6M |
| 2.60.40.10 | 1CDG / 1TEN |
| 2.60.40.1180 | 1CDG / 1SMD |
| 3.20.20.80 | 1CDG / 1SMD |
| 2.40.128.20 | 1BEB / 1MUP |

These comparisons were made across PDB-level CATH domain lists, and selected biological assemblies should be checked against exact CATH chain/domain identifiers. CATH homologous superfamily classification tests remote evolutionary/structural relationships not caught by a simple sequence-identity threshold. This is an argument for deduplication or grouping by superfamily before considering N=48 as independent.

## Preparation and software compatibility
- All 48 entire biological assembly 1 structures were normalized through Gemmi to protein-heavy-atom PDB, with unique single-character chain remapping.
- Actual `interface_engine.analyze_all_maps` succeeded and returned 11 maps for each candidate; prepared PDB SHA256 and detailed map results are in the attached GitHub Actions artifact.
- Missing side-chain heavy atoms and altlocs remain relevant to the scientific quality of surface descriptor estimates. Software execution is not evidence that these input uncertainties are harmless.
- Legacy 5CSC, and all quality-flagged or homologous entries, need exclusion/replacement/explicit sensitivity handling before the panel is labeled fully approved.
- Actual final validation is **not complete**: no adsorption prediction claim or physical parameter validation was made here.

## Traceability
Official wwPDB validation report information: https://www.wwpdb.org/validation/validation-reports
CATH domain hierarchy: https://www.cathdb.info
Run output: https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/44981110375
