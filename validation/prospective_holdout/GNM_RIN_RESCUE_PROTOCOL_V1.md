# GNM/RIN accessibility-rescue ablation V1

Status: exploratory, pre-specified before rescue scores are computed. This does not modify frozen IS-v1-Core-primary.

## Systems
- GB3 / 2OED chain A / citrate-AuNP / anionic / pH 6.4 / experimental K4,K13,K50.
- beta-2-microglobulin / 1JNJ chain A / citrate-AuNP / anionic / pH 7.7 / experimental residues 2,3,26,28,29,30,33,55,56,58,59.
- HSA / 2VUF chain A / PAA-Fe3O4 / anionic / pH 7.0 modeling convention / experimental peptide regions 373-389,403-410,414-428.

GB3 pH 6.4 is a source-condition correction made before rescue scoring; the earlier 7.0 registry value was not source-matched.

## Fixed structural descriptors
- side-chain RSA: current InterfaceScout Shrake-Rupley implementation, 200 points, probe 1.4 A.
- canonical exposure threshold: scRSA >= 0.05.
- GNM: C-alpha Kirchhoff network, 10 A cutoff; mobility from pseudoinverse diagonal, normalized within protein.
- RIN: C-alpha graph, 8 A cutoff; primary rescue descriptor is degree. Closeness and betweenness are diagnostic only.

## Rescue rules (fixed before seeing rescue results)
- Core: canonical exposed residues only.
- GNM-rescue: canonical exposed residues OR residues in the top quartile of GNM mobility.
- RIN-rescue: canonical exposed residues OR residues in the bottom quartile of RIN degree (peripheral packing).
- GNM+RIN-rescue: canonical exposed residues OR residues satisfying BOTH top-quartile GNM mobility and bottom-quartile RIN degree.

For a rescued residue below the canonical exposure threshold, effective accessibility is floored at 0.05; otherwise its measured scRSA is retained. Chemistry compatibility and pH-state logic are unchanged. Patch persistence is recomputed on each variant's candidate set using the canonical 5/8 A radii.

This floor is an exploratory accessibility-rescue convention, not a new frozen primary model. No thresholds or formulas may be changed after results are observed within this V1 experiment.

## Evaluation
- residue-level systems: AUROC, AP, Recall@5/10/20, spatial recovery at 5/8/10 A, and mapping of experimental residues into the candidate set.
- HSA: peptide-region spatial recovery is primary; residue-level metrics are secondary/descriptive because the experimental evidence is peptide-level.
- Report whether rescue improves GB3 without materially degrading beta-2-microglobulin and whether it helps or fails on HSA.
