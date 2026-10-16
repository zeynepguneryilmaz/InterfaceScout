# Residue-level structural QC execution state — 2026-10-16

- Frozen 285 source / 282 V1 successes unchanged.
- 239 originally STANDARD, but this original tag uses median pLDDT >=70 and <30% residues with pLDDT <50; it is not a high-confidence patch guarantee.
- A stricter post-hoc screening (median pLDDT >=90, fraction below50 <=5%, coverage 0.95–1.05) retains **97 provisional proteins**.
- The original 12 shard artifacts were successfully acquired and 282 protein JSON / 3102 rank-1 channel records audited for positions.
- In provisional 97, 320/1067 rank-1 centers are in first 30 precursor residues. HSA precursor signal/propeptide overlaps 6/11 top patches; real mature-chain suitability needs residue-specific checks.
- A separate GitHub Actions workflow `.github/workflows/posthoc-structural-qc-2026-10-16.yml` has been committed. It downloads archived input artifacts and invokes `scripts/validation/patch_residue_confidence_qc.py` to fetch original AlphaFold PDB and UniProt processed-region annotations and calculate residue-wise local pLDDT QC.
- **Workflow execution and its result were NOT verified at the time of this record.** Do not report patch-wise PASS counts or publish a final high-confidence N until its result artifact is inspected.
- Offline local environment lacked DNS to AlphaFold and UniProt. Existing metadata and patch JSON do not contain residue-wise pLDDT for all selected patch members.
- This remains a **post-hoc structure suitability analysis**, not a change to frozen V1 weights, thresholds or validation ground truth.
- Important limitation: passing local pLDDT and processing annotations still cannot certify relative domain orientations (PAE), oligomeric state, mature-chain conformation, or experimental adsorption-interface correctness.
