# ISV active panel selection — definitive removal of four PDB entries

**Decision date:** 2026-11-24. **Active ISV candidate list:** [ISV_active_candidates_76_v3.csv](ISV_active_candidates_76_v3.csv).

The 80-candidate structural screening is preserved as historical provenance, but **only 76 PDBs are active in ISV from this decision onward**. The requested five-record review resulted in **keeping 1STP** and **definitively excluding the other four**. No automatic substitution with a different PDB for these four has been authorized.

| PDB | Protein | Previous ISV classification | Final ISV selection |
| --- | --- | --- | --- |
| **3GHG** | Human fibrinogen | Serum / corona | **EXCLUDED** |
| **1QJS** | Hemopexin | Serum / corona | **EXCLUDED** |
| **1STP** | Streptavidin | Biointerface / immobilization | **RETAINED**, conditionally interpreted as biotin-bound-derived protein-only tetramer |
| **4DN0** | PelD cytosolic c-di-GMP-binding domain | Biofilm | **EXCLUDED** |
| **6LQH** | CsgG–CsgF outer-membrane secretion complex | Biofilm | **EXCLUDED** |

**Interpretation:** Exclusion means *removed from all future active ISV panel aggregates*, not that the original experimental structures are intrinsically invalid or that their prior computation must be deleted. Their original row data, surface QA, structure hashes and production results remain visible in the archived [80-entry candidate file](ISV_protein_candidate_pool_80_v2.csv), [80-structure full screening report](ISV_80_FULL_SCREENING_REPORT.md) and [surface/cofactor audit](ISV80_SURFACE_ASSEMBLY_COFACTOR_RESULTS.md).

The retained **1STP** has one missing-sidechain residue in the union of all 11 maps' Top-5 patch members and requires the **biotin-bound structure, biotin-removed representation** distinction. Retention does not constitute full scientific acceptance of ligand-sensitive conclusions.

## Current counts and statistics policy

- **76** distinct active PDB candidates (original 80 minus four excluded).
- **12 → 10** biofilm-related candidates retained, because 4DN0 and 6LQH were in that subgroup.
- **80 historical technical executions and their 880 map outputs remain archived**, not recalculated as 76.
- The 76 active candidates have previously passed the 11-map computation **as part of the original 80-run screen**, but their **76-only aggregate must be computed separately** before reporting any panel summary.
- **TV48** remains unchanged and remains the other active dataset.
- Three near-sequence-similar pairs among the original panel should be recognized in group-level statistics; the active list retains all those pairs.
- **Neither 76 active candidates nor 76 model-compatible structures implies 76 independent scientifically accepted validation cases.** Final source matching, exposed-atom chemistry and biological scope checks remain.

**Machine-readable rule:** all downstream ISV analysis scripts and manuscript tables should load `validation/ISV_active_candidates_76_v3.csv` as the active inclusion list. The 80-entry historical file is strictly historical and must not be used to calculate current ISV sample sizes without explicit labeling.
