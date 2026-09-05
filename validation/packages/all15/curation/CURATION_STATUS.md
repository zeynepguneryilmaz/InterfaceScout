# InterfaceScout 15-protein benchmark curation status

Date: 2026-09-05
Branch: `validation-study`

## Authoritative curation records

Use these files for all subsequent model-development analyses:

1. `assembly_map.csv` — biological-unit/chain policy for all 15 proteins.
2. `structure_filter_rules.csv` — explicit coordinate filters. Critical filters: FNIII8-10 = 1FNF chain A residues 1236-1509; WNV E-domain III = 2HG0 chain A PDB residues 8-109.
3. `primary_quantitative_conditions_v2.csv` — authoritative primary residue-level model-selection conditions and ground-truth lists.
4. `source_condition_audit_v2.csv` — all registry conditions, source/assumption status, and permitted analysis role.
5. `source_condition_verified_updates.csv` and `provenance_corrections.csv` — audit trail showing why records changed.
6. `tier2_verified_conditions.csv` — verified conditions/interpretation for orientation/activity-only systems.

## Primary quantitative development pool

12 positive conditions across 7 unique proteins:
- Fibronectin III8-10: FN_COO, FN_NH3, FN_CH3_HEAD, FN_CH3_SIDE, FN_CH3_BETA
- Cytochrome c: CYTC_CH3, CYTC_COOH
- Lysozyme: LYZ_SWCNT
- Alpha-chymotrypsin: CHT_CNT_HYDRO
- Myoglobin: MB_CIT_AUNP (citrate-anionic proxy caveat)
- WNV E protein domain III: WNV_GRAPHENE
- RNase A: RNASE_SILICA_4NM

Repeated conditions on one protein are not independent statistical units. Protein-level macro summaries and protein-clustered bootstrap are mandatory.

## Separate evidence sets

Curvature stress tests (not pooled with primary model-selection conditions):
- CYTC_SILICA_4NM
- CYTC_SILICA_11NM
- RNASE_SILICA_11NM

Intentional chemistry negative control:
- CHT_CNT_PI_NEG

Tier-2 orientation/activity systems (not residue-AUROC labels):
- RNASE_COOH_SAM, RNASE_NH2_SAM
- ACHE_COOH_SAM, ACHE_NH2_SAM, ACHE_CIT_AUNP
- GB3_CIT_AUNP
- PK_CIT_AUNP
- CA2_CIT_AUNP
- FIB_CIT_AUNP

Excluded/pending from model selection:
- BSA_TIO2, BSA_SILICA, BSA_GRAPHITIC (exact ground-truth/source audit not frozen)
- FN_OH, CYTC_OH (H-bond ontology semantics not frozen)
- CYTC_NH2 (exhaustive residue ground truth not frozen)
- HB_CIT_AUNP, TRP_CIT_AUNP (paper provides contact-probability/binding-region evidence, but an exact binary residue list has not been frozen without digitization)
- FIB_GRAPHENE (ground-truth audit incomplete)

## Critical provenance corrections

- The previous fibronectin DOI `10.3390/ma11122570` was incorrect and unrelated to fibronectin. Correct SAM source: Liamas et al., *Int J Mol Sci* 2018, DOI `10.3390/ijms19113321`.
- Correct fibronectin source construct is 1FNF residues 1236-1509, not the full deposited chain.
- Fibronectin source conditions: pH 7; 1.0 M NaCl for carboxyl/negative trajectory, 0.8 M NaCl for amine/positive trajectory, 0.05 M NaCl for uncharged systems, 300 K.
- WNV graphene paper models E-domain III only (original E residues 299-400), mapping to 2HG0 PDB residues 8-109. The earlier whole-ectodomain WNV run is invalid for benchmark validation.
- WNV frequent-contact binary ground truth uses >35% regions: original E 346-352 and 397-400, mapping to 2HG0 PDB 55-61 and 106-109.
- RNase 4-nm silica exact key adsorption residues are 15,16,17,18,21,50,51,52,53,55 and are now a primary quantitative development condition.

## Historical results that must not be used as final curated validation

`results_full_capacity_fast/` is a useful historical screening snapshot, but it used deposited structures without the present chain/range curation, SASA=120, and mostly screening pH/ionic assumptions. Do not cite its aggregate predictions as final benchmark results.

`results_curated_validation/` and `results_curated_validation_v2/`, if present, are intermediate curation runs. The source-corrected optimized run `results_curated_validation_v3_fast/` is the current model-selection run, provided its workflow completes and QC passes.

## Statistical interpretation

All current benchmark outcomes were inspected during model development. Therefore none of these analyses is a prospective independent validation, even after curation. After the feature/rule set is frozen, a new untouched protein/material holdout set is required for prospective validation.
