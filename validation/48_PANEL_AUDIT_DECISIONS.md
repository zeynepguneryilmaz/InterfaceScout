# InterfaceScout 48-structure coordinate audit — screening decision (2026-11-22)
This is a **screened candidate development panel, NOT 48 final accepted structures**. No adsorption results or InterfaceScout map rankings were used in replacement decisions.

## Executed independent data checks
- GitHub Actions run: https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/44968292207
- 48 original PDB entries and 23 unique backup entries examined using RCSB metadata plus downloaded asymmetric unit and assembly-1 mmCIF coordinates.
- 48/48 original entries were downloadable and parseable, 47/48 original X-ray resolution ≤2.5 Å; 5CSC is a disclosed 2.8 Å historical exception.
- Pairwise sequence alignment found **zero qualifying matches between the 48 original structures** with local identity ≥30% **and coverage ≥70% of the shorter sequence**. This is a specific screening rule, NOT proof of non-homology; domain-level and structural classification still needed. Short local alignments can generate false high-identity hits; none should be interpreted without length and coverage.
- R-free missing from early entries is flagged **missing**, not counted as a structural failure. DCC/wwPDB validation statistics should be checked separately.
- Biological assembly labels and mmCIF coordinates retrieved; chain completeness requires per-chain mapping of SEQRES to observed coordinates, not only totals.

## Clear changes to the original list (4 swaps)
| Remove | Reason | Replace with (screened candidate) | Replacement rationale |
|---|---|---|---|
| 1MSO | A/B insulin entities only 21/30 residues; annotated dodecamer | 1C5E | 1.1 Å, R-free 0.133; trimeric bacteriophage lambda head protein D; check assembly |
| 1FNF | R-free 0.308 beyond conservative >0.30 review/exclusion threshold | 1PDB | 2.2 Å, R-free 0.250, monomeric DHFR; antifolate/NADPH bound, chemical state annotation required |
| 1HML | Observed CA fraction 0.866 (<0.90 screening threshold) | 1A1X | 2.0 Å, R-free 0.253, CA fraction 0.981; dimeric MTCP-1; check assembly |
| 1CON | Cd-substituted concanavalin A; changed coordination chemistry | 1ECA | 1.4 Å, monomeric erythrocruorin; heme-bearing globin, avoid treating cofactor environment as generic metal-free surface |

The replacements are provisional pending domain-family balance and exact per-chain assembly selection. 1ECA overlaps *globin fold* coverage already present in the panel, even if sequence homology did not trigger the screening threshold; a more structurally diverse replacement is preferable if available.

## Remaining critical review decisions
- **1ATP, 1PPB**: complexes include 20- and 36-residue auxiliary protein entities; either define a reproducible isolated-chain analysis (not the original whole assembly) or replace them. Do not silently count these as passed.
- **1DPS, 1RBL**: 12-/16-chain biological assemblies are important for geometric stress tests, but require assembly-level surface processing and tractability checks.
- **1BRS, 1CSE, 1A4Y**: heteromeric complexes require consistent treatment of protein partners and bound cofactors; decide whether to keep as explicit 'complex' stratum.
- **5CSC**: retain only as historical 2.8 Å exception or replace for a strict-quality 48 set, while keeping historical analysis in the appendix.
- Residue-level domain homology, side-chain completeness, biological assembly choice, ligand handling, and wwPDB validation report remain to be completed. Results in this audit should not be portrayed as final QA acceptance.

## Scientific gate before parameter selection
Every included PDB must have: selected biological unit/chain IDs, experimental method and quality, coordinate coverage by selected chain, sequence nonredundancy at domain level, all nonstandard residues/ligands policy, and a reproducible prepared input PDB/CIF with SHA-256. Only then freeze the 48-row manifest and begin parameter optimization.
