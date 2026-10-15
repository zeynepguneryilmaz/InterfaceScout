# Independent computational concordance — protocol draft (2026-10-15)

This is **NOT** experimental Level-3 accuracy or an update to V1. Original seven-condition benchmark, prospective two-condition benchmark, and 285-protein breadth outputs remain frozen.

## Purpose and safeguards
Compare *previously published independent* residue-level MD/simulation material contact annotations with the locked V1 protein-side patch predictions. Report as **Level-C computational concordance** only. Literature simulation outputs cannot become experimental gold-standard labels. Lock sources, conditions, residue annotations, protein sequence/PDB mappings, exposed surface chemistry mapping and all inclusion/exclusion decisions **before inspecting InterfaceScout predictions**. Freeze labels on a new branch, then use the existing frozen scoring ordering and spatial hit/recall metrics on a separate scoring branch; no post hoc threshold changes or channel fitting.

## Inclusion criteria
1. Published protein/material pair and explicit experimental or modeled surface material chemistry; distinguish coating from core.
2. Source-defined protein-contact residue list or contact frequency trace from simulations, not qualitative explanation or merely unfolded regions.
3. Matching species, construct, sequence indexing and experimental structure suitable for exact residue mapping. Verify protonation/pH.
4. Adsorbed state native-like or initial-contact frame; separately flag rearranged/unfolded simulations.
5. Simulation independent of InterfaceScout; any paper that calibrates against the same experimental benchmark case gets a **related-evidence** tag and is not claimed to be fully independent validation.
6. Cluster study × protein × material family and separate conditions (size, chemistry, pH) without inflating independent N.

## Primary descriptors (no claim of experimental accuracy)
- First direct-overlap rank and first near-5/near-8 Å rank.
- Top-3/top-5 near-8 Å hit and recall; top-5 direct-overlap hit and recall; best Pareto-front-1 near-8 Å recall; source contact-site propensity mean/max.
- Agreement by contact-frequency or source confidence *only* if source publishes reproducible numeric frequencies. Do not use pseudo-precision on qualitative sites.
- Protein/material family-clustered summaries, contact label size and surface-area controls, evidence dependence table; discordant examples remain included.

## Literature leads — UNSCORED pending pre-score audit
- C001: HEWL–silica, MD + MP-SPR, pH 7 context, DOI 10.1039/C5CP03910J. Need exact source MD contact residue table, PDB and defined silica model; do not use SPR alone as residue labels.
- C002: cytochrome c/RNase A/lysozyme–silica 4 and 11 nm, published atomistic MD, DOI 10.1371/journal.pone.0107696. Paper has explicit cytochrome c residue contacts; verify PDB constructs, detailed condition pH, 4 nm vs 11 nm and potential overlap with experimental Shrivastava silica family. Count related clusters conservatively.
- C003: lysozyme–SWCNT, Scientific Reports 2025, DOI 10.1038/s41598-025-96435-3. Original paper names Arg61 Trp62 Trp63 Arg68 Thr69 Pro70 Arg73 Asn103 Asn106 Ala107 Trp111 Arg112 Asn113 Lys116; verify source species PDB numbering tube chirality/functionalization, simulation conditions and frozen chemistry mapping before any scoring.

## Status
Discovery only. Zero additional ADMIT/ground truth/scoring in either experimental Level-3 or computational concordance. This file authorizes no changes to frozen thresholds/weights/claims.
