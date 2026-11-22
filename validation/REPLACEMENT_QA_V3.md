# InterfaceScout 48-protein development panel — replacement QA revision v3
Date: 2026-11-22.

## Decision
The original six provisional replacement candidates were inspected using RCSB experimental metadata, coordinate mmCIF (asymmetric unit and biological assembly 1), polymer sequence lengths, per-chain C-alpha coverage, presence of standard side-chain heavy atoms among *observed* standard residues, alternate-location status (after correcting a null-altloc detection bug), and sequence screening.

**Five candidates retained after this documented quality screen:** 1C5E, 1PDB, 1A1X, 1STN, 2CBA. **1P32 rejected** because chain coverage was 182/209, 171/209, 176/209 (87.1%, 81.8%, 84.2%). **Selected replacement:** 1RWZ (archaeal PCNA ring-shaped trimer).

## Candidate QA
| PDB | X-ray resolution (Å) | R-free | Assembly | Observed / polymer length (per selected protein chain) | Missing standard side-chain heavy atoms in observed residues | Altloc residues | QC decision |
|---|---:|---:|---|---|---:|---|---|
| 1C5E | 1.10 | 0.1328 | trimer | 95/95 each x3 | 0 | 2, 3, 1 | Pass structural screening |
| 1PDB | 2.20 | 0.250 | monomer | 186/186 | 0 | 0 | Pass structural screening; bound state recorded |
| 1A1X | 2.00 | 0.253 | dimer | 106/108 each x2 | 0 | 0 | Pass structural screening; symmetry-generated chain ID |
| 1STN | 1.70 | not available in refinement field | monomer | 136/149 | 0 | 0 | Pass coordinate coverage with R-free/legacy QC caveat |
| 2CBA | 1.54 | not available in refinement field | monomer | 258/260 | 0 | 4 | Pass coordinate coverage with R-free/metal caveat |
| **1RWZ** | **1.80** | **0.238** | **PCNA trimer** | **244/245 each x3** | **0** | **0** | **Pass structural screening; selected for original 1CON slot** |

1RWZ is not a membrane protein. Its PCNA (sliding clamp) fold contributes a distinct trimeric, ring-shaped surface architecture, preserving the reason a toroidal protein was sought. RCSB SCOPe annotation identifies a DNA-clamp fold. RCSB: https://www.rcsb.org/structure/1RWZ ; https://www.rcsb.org/experimental/1RWZ ; https://www.rcsb.org/annotations/1RWZ.

## Evidence
- Latest corrected candidate-specific QA: https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/44976797742
- Most recent coordinate+sequence candidate screening: https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/44976805224
- QA code: validation/qa_six_replacements.py; sequence screen: validation/audit_rcsb_coordinates.py.
- 1RWZ biological assembly 1 mmCIF SHA-256 and remaining candidate records are in the workflow artifact and logs. The first local-sequence screen reports no qualifying 1RWZ-to-original-panel match under ≥30% aligned identity and ≥70% coverage of the shorter polymer; **this does not replace domain/fold-level CATH/SCOPe redundancy screening**.

## Important non-passes / incomplete work
This manifest is a **48-entry selected, screened PDB panel**, not a set of 48 fully validated prepared inputs. Specifically:
1. Per-chain observed vs declared sequence registration and internal missing loop locations still need a full structural alignment and inspection across 48.
2. wwPDB per-structure validation (clashes, torsions, Rfree archival caveats) and domain/fold classification must be reviewed. Missing Rfree is *not* silently interpreted as passing a numeric Rfree threshold.
3. 1RWZ and 1A1X biological assemblies contain generated multi-character chain IDs in mmCIF; prepared PDB input must uniquely remap and verify chain identities.
4. All alternate conformation residues require reproducible occupancy-based selection. Four 2CBA residues show alternate conformers.
5. Protein-only preparation drops ligands/metals. The corresponding physical interpretation is reduced; keep cofactor/ligand state as metadata, without pretending the protein was measured in an apo state.
6. Full InterfaceScout parameter benchmarking, wwPDB quality acceptance, and independent validation are **not yet performed**. No assertion about parameter performance is made.

The six PDB IDs in the v3 manifest were chosen **without inspecting InterfaceScout physicochemical map or benchmark outcomes**.

Panel currently contains 48 distinct IDs. Original 8 development structures retained (5CSC historical quality exception; analyze 5CSC-excluded sensitivity).
