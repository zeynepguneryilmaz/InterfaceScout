# Targeted quality audit of the six InterfaceScout replacement candidates

Date: 2026-11-22. Source: GitHub Actions https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/44969195735 (first successful targeted mmCIF audit); altloc correction rerun https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/44969257801.

**Scope:** X-ray method, reported resolution/R-free, assembly-1 mmCIF coordinate availability, entity-level sequence lengths and observed amino-acid/C-alpha coverage, missing canonical side-chain heavy atoms among observed residues, SHA-256 for assembly 1. This is NOT a comprehensive wwPDB validation or evidence of model performance.

| PDB | Resolution Å | R-free | Assembly | Observed / sequence (per chain) | Missing standard side-chain heavy atoms among observed residues | Decision |
|---|---:|---:|---|---|---|---|
| 1C5E | 1.10 | 0.1328 | trimer | 95/95 for A,B,C | 0 detected | Retain, conditional on structure and assembly prep |
| 1PDB | 2.20 | 0.250 | monomer | 186/186 | 0 detected | Retain, annotate NADPH/ligand-bound state |
| 1A1X | 2.00 | 0.253 | dimer | 106/108 both chains | 0 detected | Retain, differentiate symmetry-generated A-2 chain |
| 1P32 | 2.25 | 0.234 | trimer | 182/209, 171/209, 176/209 | 0 detected among resolved residues | **Reject from main panel**: none of the three chains meets the prespecified 90% coordinate completeness threshold; do not impute missing residues |
| 1STN | 1.70 | unavailable in RCSB refine endpoint | monomer | 136/149 (91.3%) | 0 detected | Conditional retain; missing termini and DCC validation to inspect |
| 2CBA | 1.54 | unavailable in RCSB refine endpoint | monomer | 258/260 (99.2%) | 0 detected | Conditional retain; zinc-bound catalytic-site interpretation |

**Interpretation of quality gates:**
- No missing canonical side-chain atoms *among the observed residues* does not mean no missing residues. For 1P32, entire stretches lack coordinates and cannot be accepted as complete even though observed side chains are intact.
- R-free null means **not available in the queried refinement field**, not R-free=0 or an unacceptable experimental quality.
- An initial altloc detector mistakenly flagged all residues due to a null-character representation issue. This field is excluded from this report until the corrected follow-up completes.
- For 1A1X the second generated chain is named A-2 in mmCIF; a PDB writer with one-character chain IDs must remap chains and maintain unique residue keys.
- Biological assembly selected by RCSB is a testable analysis choice; it is not guaranteed to be physiological in all conditions.
- Domain homology beyond the recorded pairwise sequence screen, wwPDB clash/geometry reports, selected-chain per-residue mapping, residue altloc treatment, ligand/metal treatment, prepared-input SHA and full InterfaceScout run have not all been completed.

**Scientific decision:** Of the six selected replacement candidates, *five remain provisionally usable* and **1P32 should be replaced** under the current coverage gate. No claim of six fully certified PDB structures is justified.

Replacements need to pass the same audited gates. Do not silently substitute 1ECA (globin redundancy) or 1A32 (R-free 0.318). 
