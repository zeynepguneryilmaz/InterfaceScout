# InterfaceScout — Active validation datasets

**Active scientific plan (v1.1.0): exactly two datasets.** This is the authoritative index for current validation and manuscript planning.

| Dataset | Purpose | Membership | Status |
| --- | --- | --- | --- |
| **Dataset 1 — Technical verification (TV48)** | Verify numerical implementation, structure processing, all 11 chemistry maps, dual tRSA/scRSA, glycine handling, patch construction, Pareto ranking, criterion dependence, and stability | **48 frozen, quality-screened protein structures** from [the accepted 48-structure manifest](strict_apo_48_selected_candidates_v2_final.csv) | **Completed for the tested v1.1.0 code:** 48/48 structures, 528 maps, 12,068 residue records, 2,112 objective ablations, 21,000 perturbation trials |
| **Dataset 2 — Initial scientific validation (ISV)** | Assess how interpretable protein-side surface descriptors and candidate regions are under independently documented structural/experimental evidence, using fair reference methods | **To be selected, checked and frozen**; no final protein list or final sample size is claimed yet | **Not yet constructed or validated** |

## Dataset 1 — Technical verification: TV48

**Frozen structure list:** [strict_apo_48_selected_candidates_v2_final.csv](strict_apo_48_selected_candidates_v2_final.csv)

**Acceptance checks:** biological-assembly source SHA-256, residue identity and geometry, tRSA surface gate, non-backbone scRSA, correct GLY treatment, 11 maps per structure, candidate patches and reproducible ranking, objective correlations, rank sensitivity and full exports as applicable.

**Numerical records:**
- [Production calculations and objective ablation](FINAL_V110_RESULTS.md)
- [Extended correlation, robustness and region-analysis results](EXTENDED_V110_SCIENTIFIC_RESULTS.md)
- [Permanent reproducible evidence bundle](evidence/extended_v110_full_evidence.zip)

**Interpretation:** passing TV48 establishes tested computational correctness and structural consistency. It does **not** establish experimental adsorption-site prediction accuracy.

## Dataset 2 — Initial scientific validation: ISV

**Selection and prospective audit rules (before interpreting outcomes):**

1. Specify the actual objective: interpretation of map-localized physicochemical regions and surface-region analysis, **not** an adsorption-energy, adsorption-pose or kinetic prediction contest.
2. Require original publication support for each experimental localization assertion, with the cited figure/table/supplement, evidence resolution and all relevant conditions. A signal that only measures total adsorption cannot be labeled a residue-resolved ground truth.
3. Confirm exact species, sequence, construct, chain, residue numbering, relevant glycans/cofactors, missing segments, conformational state and PDB correspondence. Record any justified transfer mapping.
4. Choose sufficiently varied proteins, region chemistries and structural classes; prevent near-duplicates/homologs and the same experimental system from acting as independent replicates.
5. Avoid evidence leakage: do not select channels, tune thresholds or choose examples *after* observing how a candidate ranks against experimental labels. Freeze evaluation rules before producing performance tables.
6. Disclose overlap with TV48 and any previous development or retrospective literature sets. For claims of independent performance, use truly held-out evidence and protein-family-aware independence; exclude overlapping systems from such a claim.
7. Prespecify baselines and statistical unit. Use comparable region/residue coverage for geometry- and accessibility-based methods, report Top-k recovery and spatial overlap with uncertainty clustered by protein/source, and separately describe any broad evidence regions.
8. Record both positive and negative cases, uncertainty and failure modes. Absence of evidence is not evidence of a true negative interface.

**Deliverables once selected:** a frozen `ISV_manifest.csv` (PDB, assembly SHA, original experimental reference and evidence coordinates, chain and residue mapping, evidence tier, selection decision), per-protein outputs, reference-method comparisons, data-derived figures and a validated reporting table.

**Expanded ISV candidate pool:** [80 proposed protein/PDB entries](ISV_protein_candidate_pool_80_v2.csv), including **12 biofilm-related candidates** spanning surface adhesins, extracellular matrix proteins, polysaccharide-processing enzymes, intracellular regulators and challenging amyloid/membrane assemblies. These are provisional candidates, **not 80 accepted structures**. The original [68-candidate RCSB record audit](ISV_68_RCSB_PRELIMINARY_AUDIT.md) remains historical screening provenance; [biofilm structure checks](isv_biofilm_rcsb_screen.py) are separate. The old eight-candidate evidence-focused shortlist is not an additional active dataset.

**Current structural metadata audit:** [68-protein RCSB audit and screening status](ISV_68_RCSB_PRELIMINARY_AUDIT.md). The source records were verified for all 68 proposed PDBs: 55 records had no initial metadata flags, 13 need review. Neither count indicates final per-chain coordinate QA or admission.

**Current state:** [Initial ISV source-supported candidate shortlist](ISV_CANDIDATE_SCREENING.md) and [candidate screening table](ISV_candidate_screening_v1.csv) have been created: **eight distinct protein candidates and one alternate structure for an existing candidate**. None is finally accepted; the exact experimental residue annotations and PDB/construct matching remain to be audited. No numerical ISV validation metrics should be claimed before the final manifest is frozen and analyzed.

### ISV80 structural-screening results (current)

The [80-candidate production and structural screening report](ISV_80_FULL_SCREENING_REPORT.md) documents **80/80 successful 11-map engine executions**, 75 deposited principal chains above the preliminary 90% C-alpha coverage screen, and three conservative quality-review strata (10 priority, 67 conditional, 3 special challenge). **All 80 remain in the candidate pool, but no 80-member final ISV manifest is asserted.** The [per-protein triage table](evidence/ISV80_TRIAGE.csv) and [archived raw evidence](evidence/ISV80_COMPLETE_SCREENING_EVIDENCE.zip) provide reproducible provenance. Additional manual structure/side-chain chemical relevance and protein-family independence review is required before final inclusion decisions.

### Current ISV working-panel decision

**Keep all 80 computationally compatible protein structures in the working panel**, without limiting the sample to 50. The [80-structure audit](ISV_80_FULL_SCREENING_REPORT.md) and [80-row versioned decision table](evidence/ISV80_TRIAGE.csv) record 880 successful chemistry maps and scientific caveats for every PDB.

These are **working structures, not 80 formally accepted independent biological validation systems**. Ten are conservative priority candidates, 67 require structural/contextual review, and three (1HC1, 6LQH, 8C50) are segregated structural challenges. Three closely related pairs ([sequence audit](evidence/ISV80_SEQUENCE_SUMMARY.json)) must be clustered by family when appropriate: 4F5S/1AO6, 3V83/1LFG, and 1B09/1SAC. There was no flagged near-sequence ISV-to-TV48 match under the documented prefilter; domain/fold independence remains subject to confirmation.

All 80 can contribute to **descriptive software applicability** where their coordinate caveats are disclosed; no structure should contribute to strict residue-level quantitative conclusions if its missing interface chemistry, assembly or construct mismatch undermines the result.

### ISV80 targeted chemical/assembly QA (completed)

[Surface chemistry, assembly and cofactor report](ISV80_SURFACE_ASSEMBLY_COFACTOR_RESULTS.md) uses verified production maps and source structure audits for all 80 candidates, including the five previously flagged structures (**3GHG, 1QJS, 1STP, 4DN0 and 6LQH**). All 80 have at least one exposed missing-sidechain residue intersecting some Top-5 candidate region. Results are **structural uncertainty flags**, not quantified map-score effects or proof of unacceptable experimental structures.

- [Per-protein summary table](evidence/ISV80_SURFACE_QA.csv)
- [Machine-readable complete result](evidence/ISV80_SURFACE_ASSEMBLY_SUMMARY.json)
- [Archived full residue and RCSB evidence](evidence/ISV80_SURFACE_ASSEMBLY_COMPLETE.zip)

**Decision:** retain all 80 in the computational working panel, but restrict strict mechanistic interpretation for the specified conditions. The biological assembly-1 model and omitted cofactor inventory have been documented; optimal in vivo assembly and energetics were not validated for all 80.

## Previous studies — preserved as historical supporting material, not active datasets

The older **eight-structure parameter-development set** is historical methodological provenance only; it is **not an active validation dataset**, and is not part of the current two-dataset work plan.

The **nine retrospective literature conditions** and earlier historical screening panels are retained for transparency and reproducibility but are **not Dataset 2**. Their labels have been consulted previously, so they are not an independent new holdout. Their results must not be silently substituted for initial scientific validation.

Historical output files and scripts remain accessible so past work can be audited; historical storage does not make these datasets active. All new manuscript validation claims must identify **TV48** or **ISV** explicitly.

## Scope control

- **Do not add a third primary dataset** without first identifying a concrete scientific deficiency that neither TV48 nor ISV can resolve.
- **Do not change the 48 accepted structures** without documented scientific and provenance justification.
- **Do not modify numerical production behavior to improve retrospective experimental agreement.**
- Future accepted ISV outputs should be generated from a frozen source commit and retained with checksums and scripts.
