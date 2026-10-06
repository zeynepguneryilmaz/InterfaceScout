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
10. Eligibility is frozen independently of downstream InterfaceScout performance. Future additions require a new dated release.

## Current frozen LEVEL 1

Two independent experimental families are currently ingested:

### PC-DB / Chou
- 597 nanoparticle experiments.
- 1,513 human-protein columns.
- 903,261 experiment × protein cells.
- 66,792 cells satisfy the source repository's own detection rule `value > 0`.
- 836,469 are non-positive under that rule.
- 1,293 source values are negative; these are retained quantitatively but **not** relabeled as detected.
- The peer-reviewed PC-DB publication reports a larger database (817 formulations / 2,497 proteins). Only the directly accessible repository subset is frozen here.

### Payne / PXD053700
- 19 corona proteomics experiment columns.
- 401 source protein-ID groups representing 446 accession strings before group ambiguity handling.
- 376 protein groups contain one UniProt accession; 25 contain multiple accessions and therefore are not forced to a single accession.
- 7,619 experiment × protein-group cells.
- 3,536 cells have source Intensity > 0; 4,083 are non-positive.
- Exact NP metadata was parsed from `NP_Database_BovOnly_v5.xlsx` for samples 31–48.
- Proteomics sample 49 has no distinct metadata row in the source workbook. The observation is retained but zeta/size/concentration metadata is intentionally blank.

Combined frozen LEVEL 1: **616 experiment columns, 910,880 observation rows, 70,328 detected/present rows, and 840,552 non-positive rows.**

The current protein registry contains 1,914 source protein entities (1,513 PC-DB names + 401 Payne protein groups). This is **not** reported as 1,914 unique biological proteins because safe cross-source identity deduplication is not yet possible.

## LEVEL 2

The current frozen structural-response panel contains direction-resolved limited-proteolysis evidence. Protected/contact-compatible regions are kept separate from newly exposed/conformational-response regions. Additional verified structural-response sources remain deferred until their raw localization tables can be acquired.

## LEVEL 3

The residue/region localization panel contains experimentally localized protein-side regions from NMR, site-directed fluorescence, crosslink/covalent-capture MS, solid-state NMR, and limited-proteolysis evidence. Broader domain evidence carries an explicit lower-resolution flag.

## Material chemistry mapping

The mapping corresponds to InterfaceScout's 11 canonical channels: anionic, cationic, hydrophobic, pi/aromatic, H-bond donor, H-bond acceptor, metal-oxide/oxide-like, calcium/phosphate charged sites, transition-metal coordination, soft-metal sulfur affinity, and phosphate-rich.

Mapping is coating-first. An internal Au or iron-oxide core is not automatically assigned a soft-metal or oxide channel when an explicit PVP/PEG/PEI/citrate coating defines the exposed surface. Multiple channels may be true. Ambiguous surfaces remain unmapped.

Across the current LEVEL-1 experiment columns, 445 PC-DB experiments plus all 19 Payne experiment labels can be assigned at least one chemistry channel (**464 mapped**); 152 PC-DB experiments remain ambiguous/unmapped. Payne sample 49 has chemistry inferred from its explicit source surface label but lacks exact physical metadata.

## Files in the repository

- `raw_sources/` — source files copied without changing source values. The Payne binary NP workbook is preserved losslessly as base64 and its reproducibly parsed CSV is included.
- `harmonized_master.csv` — repository-scale compact index retained for inspection. The **full 910,880-row long-format file** is produced by `.github/workflows/benchmark_build_full.yml` and distributed in the benchmark artifact/package because it is too large for convenient GitHub Contents editing.
- `proteins.csv` — source protein registry with ambiguity explicitly preserved.
- `experiments.csv` — 616 Level-1 experiment/sample rows.
- `material_to_interfacescout_map.csv` — transparent chemistry mapping table.
- `residue_level_ground_truth.csv` — LEVEL 3 localization evidence.
- `structural_response_ground_truth.csv` — LEVEL 2 direction-specific structural-response evidence.
- `exclusions.csv` — exclusions, deferrals, identifier ambiguities, and metadata gaps.
- `source_manifest.csv` — source-family provenance and duplicate handling.
- `benchmark_audit.xlsx` — audit workbook delivered with the final package.
- `METHODOLOGY_PROVENANCE.md` — detailed provenance and anti-cherry-picking policy.
- `checksums.sha256` and `benchmark_metrics.json` — generated with the full benchmark artifact.

## Important identifier note

The PC-DB human frequency matrix has 1,513 protein columns, while the repository's `proteinIDs.csv` contains 2,497 accessions. No safe row-aligned crosswalk is exposed. The benchmark therefore does **not** assign these accessions by position.

## Verified but deferred sources

- **BAD 2.0**: peer-reviewed DOI 10.1021/acsami.4c06759; live database table and ACS Supporting Information are verified. The SI explicitly contains BAD 2.0 in XLSX form. The full machine-readable SI bytes were not materialized reproducibly in this run, so no BAD records are silently copied or manually reconstructed.
- **Blume / PXD017052**: independent Nature Communications/PRIDE proteomics cohort verified; binary processed/raw matrix not yet materialized into this release.
- **NanoPro-3M**: publication/project verified, but the claimed raw multi-million-record dataset was not found in the accessible GitHub project; not counted.
- **Fe3O4-wheat TLiP/XL-MS**: publication verified, raw residue/crosslink localization table not acquired; no inferred rows added.

## Scoring rule

Do not edit eligibility after viewing InterfaceScout results. InterfaceScout batch scoring must consume this frozen package as-is. Any later data-source addition requires a new version/date and an auditable source-based reason.
