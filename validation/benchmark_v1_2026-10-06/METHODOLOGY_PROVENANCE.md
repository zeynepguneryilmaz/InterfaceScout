# Methodology and provenance

## Scope and anti-cherry-picking design

The benchmark tests whether InterfaceScout's material-chemistry-conditioned structural descriptors relate to experimentally observed adsorption/corona occurrence, quantitative abundance, protected/contact-compatible regions, and adsorption-associated conformational response.

The benchmark was frozen before InterfaceScout batch scoring. Existing InterfaceScout benchmark scores, PRS, harmonic ANM, and future structure-based/breakable CG results were not used to decide eligibility, labels, structures, material mappings, or source inclusion.

## Source discovery and independence

Searches covered protein-corona databases, biomaterial adsorption databases, nanoparticle-protein proteomics repositories, PRIDE, publication Supporting Information, GitHub data repositories, limited-proteolysis/LiP-MS, NMR/CSP mapping, crosslink/covalent-capture MS, and solid-state NMR.

Sources are counted by independent **experimental data family**, not by number of papers. A later analysis of the same underlying PC-DB observations is provenance, not a new dataset.

## Data ingestion rule

Raw source values are preserved. Derived/harmonized tables are separate. Source-row identifiers permit backtracking to the original matrix or experiment. Missing fields remain missing. Figures are not digitized to manufacture values. Protein accessions, structures, pH values, residue boundaries, and NP metadata are not guessed.

A source can be verified but deferred when its binary/machine-readable table cannot be materialized reproducibly. Deferred is not equivalent to excluded.

## Frozen LEVEL-1 families

### PC-DB / Chou

The peer-reviewed PC-DB publication reports 817 nanoparticle formulations and 2,497 proteins from 83 studies. The public GitHub repository exposes a directly analyzable human matrix with 597 experiments × 1,513 protein columns. The frozen benchmark uses exactly that accessible subset.

The source repository's `Protein_Frequency.R` defines detected as `value > 0`. Therefore:
- 903,261 total experiment-protein cells;
- 66,792 detected/present;
- 836,469 non-positive;
- 1,293 negative source values remain quantitative source values but are not converted into positives.

### Payne / PXD053700

The open Payne repository contains a Perseus table with 401 Bos taurus protein-ID groups and 19 corona experiment columns. The frozen matrix contains 7,619 experiment × protein-group cells, of which 3,536 have Intensity > 0.

The NP workbook `NP_Database_BovOnly_v5.xlsx` was fetched from the source repository by exact blob SHA, stored losslessly as base64, decoded in a reproducible GitHub Actions workflow, and parsed to CSV. Samples 31–48 have exact zeta potential, Dtem, functionalized hydrodynamic diameter, NP incubation concentration, and source incubation-concentration metadata. The proteomics table includes sample 49 with the same textual surface label as sample 48, but the NP workbook has no separate sample-49 row. Its physical metadata is therefore left blank rather than copied.

Among the 401 protein groups, 376 have a single accession and 25 contain multiple accessions. Multi-accession groups retain the full source group string and do not receive a guessed single accession.

### Combined Level 1

The full long-format master contains:
- 616 experiment columns;
- 910,880 observation rows;
- 70,328 present/detected rows;
- 840,552 non-positive rows.

The source protein registry has 1,914 entities before cross-source deduplication. This is deliberately not described as 1,914 unique biological proteins because PC-DB does not expose a safe accession crosswalk for the 1,513 human matrix columns.

## Protein identifiers and structures

PC-DB protein names are retained verbatim. The accompanying 2,497-accession list is not assumed to align positionally with the 1,513 human matrix columns.

Payne single-accession groups preserve the source accession. Multi-accession groups remain ambiguous.

For LEVEL 2/3, experimental PDB structures are explicitly marked as experimental. For LEVEL 1, experimental-PDB vs AlphaFold enrichment remains a later metadata-enrichment stage. Until that enrichment is completed, no predicted structure is counted as experimental and no AlphaFold-required count is fabricated.

## Material mapping

Mapping is multi-label and exposed-surface-first. Explicit coating chemistry takes precedence over core composition. For example, an Au-containing or iron-oxide-containing hybrid coated with PVP/PEG/PEI is not automatically assigned soft-metal sulfur affinity or exposed metal-oxide coordination solely because those cores exist internally.

Uncertain exposed chemistry is left `ambiguous_unmapped`. The exact mapping rule and confidence are stored in `material_to_interfacescout_map.csv`.

## LEVEL-2 label semantics

Structural response is directional and never collapsed:
- increased proteolytic accessibility → `newly_exposed/conformational_response`;
- decreased accessibility → `protected/contact_candidate`;
- direct NMR/covalent-capture/other localization → direct localization evidence where methodologically justified.

Protection is contact-compatible but not automatically identical to direct contact because stabilization can reduce proteolytic accessibility.

## LEVEL-3 quality hierarchy

- **A:** residue/peptide localization by direct NMR/CSP, covalent-capture/crosslink MS, site-directed fluorescence, or solid-state NMR with interpretable folded-state evidence.
- **B:** experimentally localized broader region/domain evidence.

A broader region is not silently scored as if every residue were a directly observed atomistic contact.

## BAD 2.0

BAD 2.0 is peer-reviewed (DOI 10.1021/acsami.4c06759). Its live database exposes fields for protein/PDB, surface concentration, solution concentration, adsorbing surface, contact angle, surface tension, pH, ionic strength, temperature, measurement method, reference, DOI, validation, and explanatory notes. The ACS Supporting Information explicitly states that the actual BAD 2.0 is provided in table/XLSX format.

Because the complete SI XLSX bytes were not materialized reproducibly in the present ingestion environment, BAD 2.0 remains a verified-deferred source and contributes zero frozen rows. This prevents partial/manual web copying from introducing hidden sampling bias.

## Blume / PXD017052

The Nature Communications/PRIDE source is a verified independent high-dimensional plasma-proteomics cohort. Its downloadable supplementary/PRIDE binary files were located, but complete processed data bytes were not materialized successfully in this release. It remains verified-deferred and contributes zero frozen rows.

## Full-master generation

The repository contains compact inspection files, but the definitive Level-1 `harmonized_master.csv` is generated reproducibly by `.github/workflows/benchmark_build_full.yml`. The workflow emits every experiment × protein observation, including non-positive records, and packages the raw sources, derived tables, source manifest, exclusions, methodology, metrics, and SHA-256 checksums.

No InterfaceScout prediction code is invoked by this build workflow.

## Release immutability

Once InterfaceScout scoring begins, this dated benchmark package is immutable for primary analysis. New sources, corrected crosswalks, or structure enrichment must be released under a new version/date; they cannot be conditionally added after seeing model performance.
