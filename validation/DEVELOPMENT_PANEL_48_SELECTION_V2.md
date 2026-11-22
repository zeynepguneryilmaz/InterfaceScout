# InterfaceScout: fixed 48-ID development panel proposal (2026-11-22)

**Status:** 48 unique selected candidate PDB IDs; **not** yet 48 fully accepted/preprocessed molecular structures. This selection is frozen against InterfaceScout model-output inspection. Final structural acceptance still requires wwPDB validation, complete residue and side-chain QC, biochemical annotation, and input conversion for all entries.

## Replacements adopted (6)
| Removed | Reason | New selection | Independent reason |
|---|---|---|---|
| 1MSO | insulin peptides 21/30 aa; dodecamer | 1C5E | 1.1 Å, R-free .133, lambda head protein D trimer |
| 1FNF | R-free .308 | 1PDB | 2.2 Å, R-free .250, monomeric DHFR |
| 1HML | 86.6% CA coverage | 1A1X | 2.0 Å, R-free .253, 98.1% CA coverage |
| 1CON | cadmium substituted lectin, undesirable metal artefact | 1P32 | 2.25 Å, R-free .234, distinctive acidic toroidal homotrimer |
| 1ATP | 20-aa secondary chain in regulatory complex | 2CBA | 1.54 Å, carbonic anhydrase II with Zn; soluble monomer; audit pending |
| 1PPB | 36-aa secondary chain in enzyme complex | 1STN | 1.7 Å, soluble staphylococcal nuclease monomer; audit pending |

1A32 was considered and rejected (2.1 Å but R-free 0.318). 1ECA was a previous temporary candidate, not selected because globin class already represented. 1F8E not selected because viral neuraminidase is a membrane-anchored enzyme; membrane cases are excluded from this project phase. 1GCI/1G6X/1A2P rejected as redundant to original proteins/complexes.

## Important pre-existing exceptions (unchanged)
- 5CSC is one of eight previously used development structures and is 2.8 Å. Report explicitly as a legacy stratum; repeat final comparison with it excluded. If strictly uniform <2.5 Å is required, replace it for primary analysis and retain legacy sensitivity separately.
- 1A4Y R-free .286 triggers review (not automatically excluded).
- 1DPS and 1RBL are large 12- and 16-chain assemblies. Their exterior and solvent-accessibility calculations cannot be treated as isolated chains; include dedicated assembly-length runtime checks.
- 1BRS, 1CSE, 1A4Y are intentionally labeled heteromeric complexes; report member proteins and separate biological-assembly-based interpretation. They do not represent a single isolated free protein.
- R-free not reported in old deposits is not counted as fail; check DCC/validation metrics instead.

## What is frozen now
The candidate **PDB IDs** and the **no-data-leakage selection logic**. Not frozen: chain IDs, alternate occupancy handling, missing sidechains, biological assembly converted to PDB, disulfide/metal/ligand interpretation; these must be documented per entry before any parameter tuning.

## Data sources and provenance
https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/44968292207
https://www.rcsb.org/
https://www.wwpdb.org/validation/validation-reports

## Next hard pass gates
1. Run coordinate audit with 1STN and 2CBA and compare all 48 final sequence entities with >=30% and >=70% coverage AND domain class overlap. Any confirmed duplicate is replaced from new prespecified candidates before running InterfaceScout.
2. For each structure choose biological assembly and explicit selected chain IDs, compute mapped observed heavy atom/CA residue completeness, terminal vs internal missing sections, and total sidechain heavy atom missingness (ideal sidechain model).
3. Remove non-protein HETATM/ligands consistently from generated **protein-only** input, but record ligand-bound crystallographic origin and test exposed-pocket sensitivity separately. No claims that stripping ligand reproduces apo biophysics. For Zn-dependent proteins flag model's omission of explicit metal-coordination chemistry.
4. Freeze prepared structure SHA-256, code commit, software versions; verify InterfaceScout handles all 48 selected units in full assembly form without chain ID collision or memory failure.
5. Only after these hard passes change status to FINAL_ACCEPTED and begin parameter optimization.
