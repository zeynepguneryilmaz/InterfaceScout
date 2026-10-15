# InterfaceScout frozen V1: integrated assessment (2026-10-16)

## Frozen 285-panel execution
285 source proteins; 282 successful model runs (98.95%); 3 with no AF model; 281 full/near-full and 1 partial success. Eleven canonical chemistry maps per success generated 3102 landscapes. Median seven distinct rank-1 patch centers per protein. Average pairwise top-patch Jaccard 0.131. Cationic vs H-bond donor Jaccard 1.000 and anionic vs H-bond acceptor 0.877. Surface-residue-count vs zero-overlap Spearman 0.422; avoid accuracy claims based on zero overlap. These are applicability and structural-landscape findings only.

## Original experimental benchmark
Original frozen benchmark_case_results.csv: primary (n=7) top-3 near-8A hit 7/7; top-5 near-8A hit 7/7; top-5 direct-overlap hit 4/7. Primary top-5 near-8A recall by case: BSA/silica 0.128; HSA/PAA-Fe3O4 0.150; HCAII/silica pH6.3 1.000; HCAII/silica pH8.3 1.000; ubiquitin/citrate Au 1.000; beta2-microglobulin/MPA Au 0.800; asparaginase/AlumOH 0.164. Secondary structural-response cases transferrin/PS and catalase/PS: 0/2 near-8A top-3 hits and 0/2 direct top-5 hits. HCAII pH conditions are same-study cluster, not two independent studies.

Important source caveats: HSA original benchmark pH6.2 represents MES crosslinking, not necessarily prior adsorption in deionized water. HCAII fluorescent orientation probes use engineered constructs with WT PDB as structural surrogate. Report qualifiers without changing frozen labels/scores.

## Separate prospective experimental Level-3
L3E001 GB3/citrate Au: K13 at pH6.5, direct-rank=2; near8-rank=2; top5 direct hit; near8 top3 recall=1.0.
L3E025 human apo-transferrin/IONP: protected residues 47-60, 237-251, 554-564 (40 residues), pH7.4; first direct rank=12, near8 rank=3; no top5 direct hit; near8 top3 recall=0.075, top5 recall=0.10; best front1 recall=0.25. Source proteolysis at 60 C is a limitation.
Combined prospective counts: n=2; top3 near8 hit=2/2; top5 direct hit=1/2. Do not pool with original benchmark without strata.

## Computational Level-C
Previously published MD contacts remain a separate concordance stratum, not experimental ground truth. Three computational source families have a separate pre-score registry. No condition has yet cleared the exact residue, pH, structure, material-map and source-independence requirements; thus no validated Level-C numerical score exists. Do not report invented computational agreement.

## All 285 reconciled
See PANEL_285_V1_AND_EVIDENCE_AUDIT_2026-10-16.csv for all source IDs joined to frozen model/QC data and first-pass evidence statuses. This is NOT a comprehensive completed source review of 285 proteins. Serum albumin (LP025) overlaps legacy HSA family; serotransferrin (LP028) overlaps the prospective apoTf family conditionally. Human lysozyme C and hen lysozyme, and ubiquitin-eS31 fusion and free ubiquitin, are not interchangeable references.

## Evidence-limited final interpretation
Frozen V1 executes broadly and yields chemistry-associated candidate patches. Present external residue-level comparisons show heterogeneous agreement, including weak partial-region recall and secondary-case misses. No evidence presently supports accuracy generalization across 285 proteins, adsorption amount, affinity, corona occurrence, kinetics, unfolding, or a unique orientation. No frozen parameters, thresholds, benchmark labels or scores changed.
