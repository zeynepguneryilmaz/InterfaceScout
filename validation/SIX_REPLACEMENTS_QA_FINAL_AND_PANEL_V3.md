# InterfaceScout 48-entry development panel — selected roster v3

Selection status: 48 unique PDB IDs; the six recently included structures passed real mmCIF coordinate screening and all 11 chemistry-map software smoke tests. The whole 48-entry panel has not yet passed the full final quality audit.

| PDB | Resolution (Å) | R-free | Assembly | Observed/sequence | Missing standard sidechain atoms | Altloc residues | 11 maps |
|---|---:|---:|---|---|---:|---|---|
| 1C5E | 1.10 | 0.1328 | trimer | 95/95 x3 | 0 | 2,3,1 | PASS |
| 1PDB | 2.20 | 0.250 | monomer | 186/186 | 0 | 0 | PASS |
| 1A1X | 2.00 | 0.253 | dimer | 106/108 x2 | 0 | 0 | PASS |
| 1STN | 1.70 | not reported in RCSB refine field | monomer | 136/149 | 0 | 0 | PASS |
| 2CBA | 1.54 | not reported in RCSB refine field | monomer | 258/260 | 0 | 4 | PASS |
| **1RWZ** | **1.80** | **0.238** | PCNA ring trimer | **244/245 x3** | **0** | **0** | **PASS** |

1RWZ replaces 1P32, which was excluded because per-chain coordinate coverage was 182/209, 171/209 and 176/209, all below a prespecified 90% cutoff. 1RWZ offers comparable ring-like geometry without the missing-residue problem.

## Reproducible GitHub Actions reports
- Metadata, assembly and pairwise similarity screen: https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/44976381278
- Detailed coordinate, sidechain and altloc audit: https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/44976377537
- Prepared full biological assembly to InterfaceScout 11-map smoke tests for all six retained candidates: https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/44976528666

## Important scope and remaining gate
- R-free absent in queried refinement field for older 1STN and 2CBA is not a failed experimental metric; consult deposited wwPDB/DCC validation reports.
- Histidine, cysteine, pH and surface compatibility are approximate chemistry maps. A successful software execution does not prove independent physical accuracy.
- 1C5E and 2CBA contain alternate conformers: predeclare reproducible occupancy selection. 1PDB is ligand-bound; 2CBA Zn-bound. The protein-only PDB preparation drops hetero atoms, so neither should be presented as a fully physical ligand-aware model.
- End-to-end smoke included all protein chains of biological assembly 1 with unique single-character chain IDs; ensure a permanent prepared-input hash is retained.
- Domain-level homology, wwPDB geometry/clash validation, comprehensive residue completeness checks and exact prepared-input reproducibility across all 48 structures still remain. The legacy 5CSC exception (2.8 Å) should be separately sensitivity-tested.

Scientific decision: retain five quality-screened candidates, replace 1P32 with 1RWZ, and use the v3 roster for the next global 48-protein acceptance audit. Do not label the whole 48 as fully validated until those remaining tests are run.