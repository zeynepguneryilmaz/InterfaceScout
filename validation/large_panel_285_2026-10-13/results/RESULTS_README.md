# InterfaceScout 285-protein large-panel results

Batch run: GitHub Actions 38614780417

- Frozen source denominator: **285**
- Successful frozen V1 runs: **282**
- No AlphaFold DB model: **3**
- Runtime analysis errors: **0**
- Full/near-full structure runs: **281**
- Partial structure runs: **1**
- Chemistry channels: **11**
- Protein x chemistry landscapes: **3102**

This is a breadth/applicability and chemistry-specific surface-landscape analysis.
It is not an adsorption amount, affinity, corona-membership, or experimental
interface-accuracy benchmark. Experimental region-level accuracy remains the
separate Level-3 validation layer.

AlphaFold pLDDT is QC metadata only and never enters V1 ranking.

Chemistry-specificity (rank-1 patch members across 55 channel pairs/protein):
- overall median Jaccard: **0.000**
- overall mean Jaccard: **0.131**
- median distinct top-patch centers: **7.0/11**
- median zero-overlap pair fraction/protein: **0.836**

Files:
- protein_status_all.csv
- channel_summary_all.csv
- chemistry_specificity_by_protein.csv
- chemistry_pair_overlap_summary.csv
- panel_summary.json

Complete per-protein JSON outputs remain in the 12 artifacts from Actions run
38614780417. No protein/channel was selected based on favorable InterfaceScout
output or PC-DB adsorption labels.
