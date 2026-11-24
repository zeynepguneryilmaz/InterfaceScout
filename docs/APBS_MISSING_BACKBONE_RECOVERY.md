# APBS missing-backbone handling — conservative N-only recovery

**Policy implemented in production** `backend/electrostatics.py`. This does not
alter the 11 canonical physicochemical maps, patch geometries or Pareto scores.

## Trigger
1. First attempt PDB2PQR/PROPKA on the unchanged normalized input.
2. On an actual PDB2PQR error, identify **strictly isolated protein residues
   with exactly one heavy atom, N, and no CA/C/O or side-chain heavy atoms**.
3. Permit automatic APBS-only omission if there are at most **four such
   records and at most 1% of original profile residues**. This limit is
   conservative engineering policy, not a fitted physical threshold.
4. Retry PDB2PQR once with only the identified ATOM records removed from
   the **temporary electrostatics input**. Do not reconstruct coordinates,
   delete other atoms, or modify original PDB or 11 chemistry maps.
5. Require a perfect PQR residue-ID correspondence for all remaining
   original profile residues. Any additional missing or extra identity fails.
6. Run real APBS and sample the actual solvent-accessible shell. All original
   residue rows stay present; omitted residues have `null` charges and
   potentials and `pqr_mapped=false`. Both the GUI and exported
   Electrostatics_Settings include the explicit affected identifiers and
   the original PDB2PQR error.
7. Annotate `quality_status=partial_input_recovery`. The solver status is
   `computed` because genuine PB solved, but the calculation is **not
   equivalent to electrostatics on the full original structure**.
   Nearby potentials carry additional structural uncertainty.

## Explicit refusal cases
- Missing atoms in a partially defined backbone (e.g., N+CA, N+CA+C,
  missing carbonyl O) are **not** deleted automatically.
- A large number or fraction of N-only residues is rejected; no local
  threshold inflation based on a specific PDB.
- A second PDB2PQR failure or any unaccounted PQR/PDB mapping gap
  leaves APBS in an explicit error state. No invented potentials or
  automatic alternate PDB substitution.

## 1GCT test case
- Original gamma-chymotrypsin 1GCT model has N-only **SER A:11 and
  SER E:11**. These are target eligible entries, subject to source and
  identity checks.
- Production integration regression `tests/real_orphan_backbone_recovery.py`
  compares identical original chemistry maps with and without APBS,
  validates exactly two omitted residue IDs and null potential values,
  checks source PDB SHA256 is unchanged, and uses 1HEL as the complete
  control.
- Run via
  [real-backbone-recovery.yml](../.github/workflows/real-backbone-recovery.yml)
  with real external APBS and PDB2PQR. Do not claim the end-to-end case
  passed until that GitHub Action finishes successfully.

## Interpretation
The fallback is a **model-input repair strategy** for exceptional,
unparameterizable fragments, not a reconstruction of experimental
coordinates. Preserve the structure quality caveat in any research
paper. The 11 maps remain independently calculated from the original
input; recovery does not tune chemistry channels or validation labels.
