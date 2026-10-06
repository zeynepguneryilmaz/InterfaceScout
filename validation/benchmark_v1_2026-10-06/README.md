# InterfaceScout independent benchmark freeze — 2026-10-06

This directory freezes a literature/repository-derived validation panel **before any InterfaceScout scoring is consulted**. InterfaceScout V1, its current benchmark scores, PRS/ANM results, and future CG/flexibility-branch results are not eligibility variables.

## Frozen eligibility rules

1. Source must be a peer-reviewed experimental publication or a trustworthy public repository directly linked to experimental data.
2. Protein and material identities must be traceable to the original experiment.
3. LEVEL 1 accepts experimentally observed corona/adsorption occurrence or source-reported quantitative abundance. Model-generated labels are excluded.
4. LEVEL 2 requires adsorption-associated structural response measured by LiP/limited proteolysis, HDX, NMR/CSP, XL-MS/covalent capture, or a comparably direct structural method. **Increased accessibility/new exposure and decreased accessibility/protection are stored as different labels.**
5. LEVEL 3 requires protein-side residue or sequence-region localization that can be mapped to a protein structure. Broad domain evidence is allowed only with an explicit lower-resolution flag.
6. Experimental PDB structures are preferred. Predicted structures may be used only when separately flagged; predicted and experimental structures are never silently pooled.
7. Material-to-chemistry mapping is rule-based and auditable. If the exposed surface chemistry is not defensibly identifiable, the record is `ambiguous_unmapped`.
8. Duplicate publications/re-analyses of the same experimental family are not counted as independent datasets.
9. Missing values are left missing. No values are estimated from figures and no accession/structure IDs are guessed.
10. Eligibility is frozen independently of downstream InterfaceScout performance.

## Benchmark levels

- **LEVEL 1:** large-scale protein adsorption/corona. The currently vendored PC-DB human matrix contains 597 nanoparticle experiments × 1,513 protein columns. The published PC-DB itself reports 817 formulations and 2,497 proteins; only the directly accessible repository subset is frozen here.
- **LEVEL 2:** direction-resolved structural response. Current directly curated set is from Duan et al. 2019 limited-proteolysis/LC-MS/MS, with protected/contact-compatible regions kept separate from newly exposed regions. Additional verified datasets remain deferred until raw localization tables are acquired.
- **LEVEL 3:** high-quality residue/region localization from NMR, site-directed fluorescence, crosslink/covalent capture MS, ssNMR, and limited-proteolysis evidence.

## Files

- `raw_sources/` — vendored source files that were directly accessible.
- `harmonized_master.csv` — normalized Level-1 protein index. Experiment-by-protein labels and source values remain losslessly encoded in the frozen raw matrix to avoid a >900k-row duplicate expansion.
- `proteins.csv` — protein registry for the accessible PC-DB subset. UniProt accessions are intentionally blank where the repository does not provide a safe row-aligned crosswalk.
- `experiments.csv` — one row per PC-DB nanoparticle experiment.
- `material_to_interfacescout_map.csv` — transparent mapping to the 11 InterfaceScout chemistry channels with rule/confidence/status.
- `residue_level_ground_truth.csv` — LEVEL 3 contact/orientation/localization evidence.
- `structural_response_ground_truth.csv` — LEVEL 2 direction-specific protected vs newly exposed evidence.
- `exclusions.csv` — exclusions and deferred sources with reasons.
- `source_manifest.csv` — source-family provenance and duplicate handling.
- `benchmark_audit.xlsx` — audit workbook delivered as a binary artifact outside this text-only GitHub branch.
- `METHODOLOGY_PROVENANCE.md` — detailed provenance and anti-cherry-picking policy.

## Important PC-DB provenance note

The accessible human frequency matrix has 1,513 protein columns, while `proteinIDs.csv` contains 2,497 accessions. No row-aligned crosswalk between those two objects is exposed in the repository. Consequently, this freeze **does not** assign those accessions by position.

## Chemistry map

The mapping columns correspond to InterfaceScout's canonical channels: anionic, cationic, hydrophobic, pi/aromatic, H-bond donor, H-bond acceptor, metal-oxide/oxide-like, calcium/phosphate charged sites, transition-metal coordination, soft-metal sulfur affinity, and phosphate-rich. Multiple channels can be true for a chemically multifunctional surface. Explicit coatings take precedence over inaccessible core chemistry.

## What is not frozen yet

BAD 2.0, PRIDE PXD017052, the full Payne NP metadata workbook, NanoPro-3M, and the 2025 Fe3O4-wheat XL-MS/TLiP localization table are not silently reconstructed. Their current status is recorded in `source_manifest.csv` and `exclusions.csv`.

## Scoring rule

Do not edit eligibility after viewing InterfaceScout results. Any future addition must be versioned as a new benchmark release with a dated manifest and a source-based reason for inclusion.
