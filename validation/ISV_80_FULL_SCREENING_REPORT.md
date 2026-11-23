# ISV — Full 80-candidate structural and InterfaceScout screening

**Scope:** Dataset 2 (ISV), protein-side initial scientific validation. **Screening, not final acceptance.** The final target is to retain **every scientifically usable protein**, not artificially limit the panel to 50 structures. TV48 remains unchanged.

## Source provenance and computed evidence

The screening covers **80 distinct proposed PDB entries**, including **12 biofilm-related candidates**. All 80 were audited at RCSB entry level, evaluated for deposited C-alpha coordinate coverage, processed as reproducible protein-only biological-assembly representations and analyzed by the unchanged production `analyze_all_maps` engine.

- [Original 80-PDB source list](ISV_protein_candidate_pool_80_v2.csv)
- [Merged per-protein decision table](evidence/ISV80_TRIAGE.csv)
- [Machine-readable summary](evidence/ISV80_TRIAGE_SUMMARY.json)
- [Complete preserved raw evidence archive](evidence/ISV80_COMPLETE_SCREENING_EVIDENCE.zip)
- [Full-screening executable](isv_full_80_qa.py), [structural triage logic](isv_triage_80.py)
- [Eight successful production-compute shards](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45270447667) and [successfully merged/archived evidence](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45272188031).

The original parent screening workflow's aggregation failed because the aggregation runner lacked an imported dependency; **all eight computation shards passed**. The archive/triage workflow successfully merged the actual saved outputs and passed all 80-record consistency checks. This difference matters: a failed parent workflow is not described as an end-to-end PASS.

## Actual numerical findings

| Screen | Result |
| --- | ---: |
| Structure records with archived 11-map result | **80 / 80** |
| 11 chemistry maps computed | **880 / 880** |
| Map instances producing at least one candidate patch | **871 / 880** |
| Map instances with no candidate patch | **9 / 880** |
| Deposited polymer best principal chain with model-1 C-alpha coverage ≥90% | **75 / 80** |
| Deposited coordinate entries flagged for coverage review | **5 / 80** |
| No automatic technical incompatibility found | **80 / 80** |

**Coverage review cases:** **3GHG** (fibrinogen; heterogeneous partial chain representation), **1QJS** (hemopexin), **1STP** (streptavidin; common truncated engineered mature construct), **4DN0** (PelD) and **6LQH** (membrane-embedded CsgG–CsgF complex). The deposited-sequence coverage number does **not** equal full UniProt protein coverage, and a legitimate, explicitly defined domain can remain suitable despite low full-length coverage.

**Missing chemical atoms:** All 80 deposited protein-only representations have at least one standard residue flagged as incomplete in its expected side-chain heavy-atom set. This is a **residue-level structural annotation**, not an indication that all 80 proteins are unusable. Incomplete side chains should be cross-referenced against surface membership and channel relevance before admission. Three structures had any missing backbone heavy atoms under the current all-chain audit: **1TEN, 1GFL and 1GCT**. These require localized inspection, not an automatic wholesale rejection.

**No experimental adsorption-site validation** was performed for this 80-panel candidate screening. Absence of adsorption labels does not bar inclusion in the general structural ISV panel.

## Preliminary structural-review tiers (NOT final inclusion or exclusion)

| Tier | Count | Interpretation |
| --- | ---: | --- |
| **Priority structural QA candidates** | **10** | Passed the conservative automatic metadata and structural warning screens; still require chain identity, functional state and sequence-family checks |
| **Conditional QA review** | **67** | 11-map engine passes; one or more caveats about resolution, C-alpha/side-chain coverage, assembly, domain scope or cofactors |
| **Special structural challenge** | **3** | Successful engine execution but high-complexity geometry or method requires separate scope/assumption review |
| **Technical execution failures** | **0** | None of the 80 proposed PDB entries failed the final 11-map execution |

The **10 priority IDs** are: **1F41, 1B09, 1HEL, 7RSA, 1S0Q, 9PAP, 1CEL, 2H5Q, 2ACY and 1CLL**. These are **not ten final accepted structures**; for example, 1S0Q requires confirmation of its deposited trypsinogen/processed trypsin representation.

The **three special-challenge IDs** are **1HC1** (large hemocyanin assembly), **6LQH** (membrane secretion complex) and **8C50** (low-resolution cryo-EM curli-like fibril). Each is computationally executable, but it may not be scientifically comparable with ordinary soluble globular proteins. Keep these in a distinct subset or reject them from core structural interpretation if assumptions fail.

The 67 remaining candidates are **not rejected**. A single missing terminal/remote sidechain, a relevant cofactor, NMR conformer ambiguity, or legitimate domain construct can be compatible with InterfaceScout when explicitly documented. Conversely, a partly absent *surface-exposed* chemical group can invalidate a local map irrespective of high global structural coverage.

## Biofilm-relevant 12-protein subset

All **12 / 12** biofilm-related proposed PDB entries completed the production 11-map run. Functional scope must remain distinct:

- **Surface and matrix association:** SdrG, BslA, TasA, RbmA, Bap1, LecB.
- **Matrix metabolism/regulation:** PelA, PelD, PslG, RbmB. These are relevant to biofilm research but do not independently evidence adsorption to abiotic substrates.
- **Challenging membrane/amyloid representations:** CsgG–CsgF and curli-like fibril. These should not be analyzed under unstated soluble-globular assumptions.

Several entries encode partial constructs; preserve their domain identity rather than describing all maps as full-length proteins. [Biofilm structural source audit](ISV_BIOFILM_SUBPANEL.md).

## Sequence and family-dependence pre-screen (completed)

[Sequence audit workflow](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45272945468) **PASS**: retrieved polymer sequences for **all 128 PDB records (80 ISV + 48 TV48)**, comprising 142 deposited protein entities. A k-mer-prefiltered local sequence-alignment screen evaluated 32 candidate entity pairs and flagged **three strongly similar ISV pairs**:

| Protein/PDB pair | Aligned identity | Coverage of shorter deposited sequence | Policy |
| --- | ---: | ---: | --- |
| Bovine albumin **4F5S** / human albumin **1AO6** | **75.8%** | **99.8%** | Retain both as distinct experimental species; cluster by albumin family in group-level statistics |
| Transferrin **3V83** / lactoferrin **1LFG** | **63.2%** | **96.7%** | Retain both, account for transferrin-family dependence |
| C-reactive protein **1B09** / serum amyloid P **1SAC** | **52.3%** | **97.5%** | Retain both, account for pentraxin family dependence |

No ISV–TV48 high-sequence-similarity pair met the prespecified flag criteria, but this **does not prove absence of remote/domain-level structural homology**. The 4-mer prefilter is a screening mechanism and could miss distant relationships. Do not treat the full 80-structure panel as 80 independent protein families.

The complete [source entity sequences and alignment audit](evidence/ISV80_SEQUENCE_EVIDENCE.zip), [sequence-pair data](evidence/ISV80_SEQUENCE_PAIRS.csv) and [summary](evidence/ISV80_SEQUENCE_SUMMARY.json) are permanently archived.

**Current use policy:** retain **all 80** as an **ISV computational working panel** because all are executable, with **10 priority, 67 conditionally interpreted and three separately analyzed structural challenges**. Only a verified subset should enter strict residue-level mechanistic comparisons after the remaining quality and family checks; the three special structures must not be silently pooled with soluble globular proteins.

## What is still needed before final acceptance?

1. **Surface-relevant side-chain completeness:** locate missing atoms per residue and cross-check whether they affect exposed chemical features or ranked patches. Prefer a complete experimental alternative where genuinely missing atoms alter interpretation; do not silently reconstruct.
2. **Species, protein, construct and chain:** validate exact identity and biological assembly, including partner chains in complexes, fibrinogen fragments and ligand-bound enzyme states.
3. **Independence and family structure:** compare the 80 entries to each other and the frozen TV48 at sequence/domain-family level. Do not claim that different PDB codes or bovine/human albumin are independent proteins/families.
4. **Quality and chemical-context audit:** use wwPDB validation reports for local side-chain geometry, high R-free/resolution cases, NMR model/ensemble policy, metals/glycans/ligands and membrane/fibril context.
5. **Freeze a final ISV manifest:** include all passed proteins with verified assembly and input hash; disclose any conditional strata, reject only scientifically invalid choices, and rerun the final frozen ISV panel.

**Current scientific conclusion:** The 80-candidate pool is **computationally feasible**, including the biofilm subgroup. The data **do not justify** treating all 80 as structurally/biologically validated or as independent experimental adsorption sites. The aim is to maximize **qualified, interpretable coverage**, not to reach a preconceived count.
