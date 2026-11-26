# Submission verification records

The paired Poisson–Boltzmann numerical audit reads the original immutable ISV76 archive, rather than a new calculation. The source audit run is [45733233254](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45733233254). Its archived result is [paired_PB_summary.json](paired_PB_summary.json). The 76 structure-specific paired CSVs contain 35,140 finite matched residue potentials. The reported aggregate statistics agree with the manuscript and Supporting Information at the reported precision. The 8C50 numerical exception is explicitly retained.

[doi_metadata.csv](doi_metadata.csv) contains 166 distinct DOI metadata records retrieved from Crossref/DataCite or the cited publisher/PubMed fallback. This is a bibliographic metadata audit, not a claim that each publication's full text was independently re-read. Online-first dates for CATH v4.4 and the AIMP2-DX2 paper differ from issue years; the manuscript uses the issue years (2025 and 2020). The published PLOS ONE title spells “Accessibilites”; the original publisher title for the phage p4-DNA paper ends with “DNA Binding”.

## Identifier correction

The original ISV76 aggregate exports contain an identifier serialization error: `1E05` appears as `100000.0` in 11 map records, 417 residue records and 417 APBS records. Structure-specific 1E05 outputs confirm the identifier and values. The correction affects the identifier only; no scientific result was changed. Original source-run archives remain immutable. See [identifier_corrections.csv](identifier_corrections.csv). The complete local submission addendum also contains corrected aggregate overlays, exact derived residue-region CSVs, and the retrospective rescoring source package; it is not represented as deposited here until the corresponding ZIP has actually been uploaded.

Historical retrospective results are distinct from the updated publication rescoring: HCAII pH 8.3 anchors [10,37], ubiquitin CSP positions [2,3,15,17,18]. These corrected label definitions do not change calculated patches or ranks.
