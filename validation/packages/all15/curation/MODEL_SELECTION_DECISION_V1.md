# InterfaceScout model-selection decision V1

Date: 2026-09-05
Branch: `validation-study`
Status: development-benchmark decision; not prospective validation.

## Decision

The current primary InterfaceScout score should remain the **Core IS** score for model freeze candidate status. GNM, APBS, radial prominence, and RIN descriptors remain available as auxiliary/diagnostic descriptors and ablation-tested alternatives, but are **not** justified as universal multiplicative factors in the primary score on the present curated benchmark.

This decision is based on the source-corrected curated benchmark in `results_curated_validation_v3_fast/`, not on the historical uncurated `results_full_capacity_fast/` run.

## Curated benchmark used for model selection

Primary residue-level pool: 12 positive conditions across 7 unique proteins:
- Fibronectin III8-10: FN_COO, FN_NH3, FN_CH3_HEAD, FN_CH3_SIDE, FN_CH3_BETA
- Cytochrome c: CYTC_CH3, CYTC_COOH
- Lysozyme: LYZ_SWCNT
- Alpha-chymotrypsin: CHT_CNT_HYDRO
- Myoglobin: MB_CIT_AUNP (citrate-anionic proxy caveat)
- WNV E-domain III: WNV_GRAPHENE
- RNase A: RNASE_SILICA_4NM

Additional evidence is kept separate:
- curvature stress: CYTC_SILICA_4NM, CYTC_SILICA_11NM, RNASE_SILICA_11NM
- intentional chemistry negative control: CHT_CNT_PI_NEG
- Tier-2 orientation/activity evidence: RNase SAM, AChE SAM/citrate-AuNP, GB3, Proteinase K, CAII, fibrinogen
- region-level secondary evidence: hemoglobin and trypsin citrate-AuNP

Repeated conditions on a single protein are not treated as independent biological replicates. Protein-level macro summaries and protein-clustered bootstrap are the inferential summaries.

## Structure and source curation applied

- FNIII8-10: 1FNF chain A residues 1236-1509 only.
- Cytochrome c: 3NWV chain A only.
- Lysozyme: 1LYZ chain A.
- Alpha-chymotrypsin: chains A+B+C as one proteolytically split enzyme; duplicate D/E/F removed.
- Myoglobin: 1MBN chain A.
- WNV E-domain III: 2HG0 chain A residues 8-109 only.
- RNase A: 7RSA chain A.

Current curated structure QC: 1FNF 274 CA residues, 3NWV 104, 1LYZ 129, 4CHA 239 across A/B/C, 1MBN 153, 2HG0 102, 7RSA 124.

Source-condition corrections include:
- fibronectin source corrected to Liamas et al., IJMS 2018, DOI 10.3390/ijms19113321;
- fibronectin pH 7 is source matched; NaCl 1.0 M (COO), 0.8 M (NH3), 0.05 M (uncharged SAMs), 300 K;
- WNV ground truth restricted to domain III and mapped to 2HG0 residues 55-61 and 106-109;
- RNase 4-nm silica key residues frozen as 15,16,17,18,21,50,51,52,53,55;
- historical whole-deposit chain-duplicated full-capacity scores are excluded from final model-selection claims.

## Aggregate performance on the curated primary pool

Protein-macro median metrics:

| Method | AUROC | AP | Recall@10 | Recall@20 | Spatial@10A_top10 | Median nearest top10 (A) |
|---|---:|---:|---:|---:|---:|---:|
| Core | 0.7111 | 0.2820 | 0.3000 | 0.3636 | 0.6429 | 7.1760 |
| Core_APBS | 0.7111 | 0.2116 | 0.3000 | 0.3636 | 0.6364 | 7.1760 |
| Core_GNM | 0.7262 | 0.1432 | 0.1111 | 0.3333 | 0.5000 | 7.4533 |
| Core_GNM_APBS | 0.6981 | 0.1432 | 0.1538 | 0.3571 | 0.6364 | 6.6048 |
| IS_candidate = Core*GNM*APBS*radial | 0.6399 | 0.1274 | 0.1538 | 0.3077 | 0.5714 | 7.8346 |
| IS_RIN_degree | 0.6607 | 0.1457 | 0.1538 | 0.3077 | 0.6154 | 8.9424 |
| IS_RIN_closeness | 0.6435 | 0.1265 | 0.1538 | 0.3077 | 0.5000 | 9.1892 |
| IS_RIN_betweenness | 0.6458 | 0.1373 | 0.1818 | 0.3077 | 0.6154 | 7.8346 |

Interpretation: GNM can raise AUROC in some proteins, but does not improve the full retrieval objective. The universal GNM multiplier reduces AP and top-k/spatial recovery. Adding APBS and radial prominence does not recover a consistent protein-level gain. RIN variants likewise do not justify additional complexity in the primary score.

## Protein-level component consistency versus Core

Selected results from `component_deltas_by_protein.csv`:

- Core_GNM AUROC: median delta +0.0281; 5 proteins improved, 2 worsened.
- Core_GNM AP: median delta -0.0500; 3 improved, 4 worsened.
- Core_GNM Spatial@10A_top10: median delta -0.1667; 1 improved, 6 worsened.
- Core_GNM_APBS AUROC: median delta +0.0329; 5 improved, 2 worsened.
- Core_GNM_APBS AP: median delta -0.0500; 3 improved, 4 worsened.
- IS_candidate AUROC: median delta +0.0397; 4 improved, 3 worsened, but its protein-macro median AUROC is lower than Core because the distribution of changes is heterogeneous.
- IS_candidate AP: median delta -0.0659; 3 improved, 4 worsened.
- IS_candidate Recall@10: median delta -0.1667; 3 improved, 4 worsened.
- IS_candidate Spatial@10A_top10: median delta -0.1364; 2 improved, 5 worsened.

The confidence intervals for most component deltas cross zero, except some spatial degradation signals. There is no robust evidence that any auxiliary descriptor should be forced into every surface/protein score.

## Protein-clustered bootstrap

Protein-clustered bootstrap (10,000 resamples) for Core:
- protein-macro median AUROC = 0.7111; 95% CI 0.6143-0.7225
- AP = 0.2820; 95% CI 0.0818-0.4644
- Recall@10 = 0.3000; 95% CI 0.0909-0.3846
- Recall@20 = 0.3636; 95% CI 0.2727-0.5385
- Spatial@10A_top10 = 0.6429; 95% CI 0.6364-0.9231
- median nearest top10 = 7.176 A; 95% CI 3.691-9.189 A

Core_GNM has slightly higher median AUROC (0.7262) but substantially poorer AP/top-k/spatial summaries, so AUROC alone must not determine the model architecture.

## Mechanistic lessons from individual systems

- FN_COO: Core is already extremely strong (AUROC ~0.991); APBS can perfect the single-anchor ranking, whereas GNM strongly degrades it.
- FN_NH3: GNM/APBS can help substantially, showing that auxiliary physics may be condition-specific.
- FN_CH3 modes: GNM/radial can improve some hydrophobic orientation-specific labels, but one chemistry-only model cannot represent orientation state explicitly; the same chemistry does not encode head/side/beta approach history.
- WNV_GRAPHENE: Core performs well (AUROC ~0.756), but GNM severely degrades it (~0.383), directly arguing against a universal dynamics multiplier.
- RNASE_SILICA_4NM: Core AUROC ~0.713 and AP ~0.484; GNM raises AUROC slightly (~0.726) while lowering AP (~0.309), again showing a metric tradeoff.
- Lysozyme-SWCNT similarly showed that static chemistry/patch signal can already be strong and dynamic/geometric multiplication can hurt.
- Citrate-AuNP conditions remain chemically coarse when represented by a generic anionic proxy; APBS effects there must not be overgeneralized.

## Primary-score freeze candidate

The primary score should therefore be the mechanistically minimal **Core IS** based on:
1. material-chemistry compatibility;
2. side-chain solvent accessibility;
3. pH-dependent residue state/compatibility;
4. multiscale 5/8 A spatial patch persistence.

GNM, APBS, radial prominence, and RIN remain output descriptors/diagnostics and may be reported as optional contextual analyses, not universal multiplicative score terms.

This is a **model-development freeze candidate**, not yet a final prospective validation freeze.

## Remaining work before definitive freeze

1. Run explicit pH/ionic sensitivity where source conditions are unresolved or continuum ionic strength is a computational reference.
2. Complete Tier-2 orientation/spatial scoring under source-corrected conditions.
3. Treat hemoglobin/trypsin as region-level spatial evidence unless exact contact probabilities are digitized/frozen independently.
4. Resolve fibrinogen large-system APBS only if electrostatics are used diagnostically; Core primary score does not require APBS.
5. Preserve BSA systems as excluded/secondary until exact ground truth is frozen.
6. After these checks, freeze the Core rule set and do not revisit it using prospective holdout outcomes.
7. Evaluate a new untouched protein/material holdout set prospectively.

## Claim discipline

The current 12-condition/7-protein benchmark is a development benchmark inspected during model selection. It may support statements about internal benchmark performance, robustness, chemistry specificity, and feature ablation. It must not be described as blind, independent, or prospective validation.
