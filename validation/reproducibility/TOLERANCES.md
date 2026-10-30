# Numerical acceptance contract (v1)

Defined before reference generation or P1 execution. These are verification
tolerances, not changes to the publication model or its ranking rules.

- Scalar residue/patch scores, SASA (A^2), distances (A):
  `abs(actual - expected) <= 1e-10 + 1e-10 * abs(expected)`.
  This accommodates float64 accumulation-order differences, not changes at the
  old 3/5/6-decimal rounding scales. A positive expected local score must remain
  strictly positive, even below the absolute tolerance.
- Non-finite values fail. Discrete identities, surface flags, membership,
  seed membership, counts, Pareto fronts and rank order must match exactly.
- PDB B-factor fields: absolute error <= `0.0050000001` on the 0..100 scale
  (half a two-decimal quantum plus float representation allowance). Atom identity,
  coordinates, occupancy and all non-B-factor fields must be preserved exactly
  relative to the prepared analyzed PDB. No tolerance authorizes extra models.
- Analytic isolated-sphere SASA: the scalar tolerance above. Independent sampled
  SASA uses the same declared 200-point quadrature; it is a verification of
  sampling/occlusion, not convergence to an exact molecular surface.
- Neighborhood boundaries retain the model's `<=` radius and strict positive
  hemisphere dot product. Verification tolerances must never alter membership.
- Ranking is checked from the actual unrounded descriptors after independently
  validating their values. Near ties are reported; tolerances are never used to
  create ties or change Pareto dominance. Exact independent scalar rank matching
  is not claimed across different floating-point summation orders.

No output becomes a reference merely because the application generated it.
The offline verifier recomputes from hash-checked PDB inputs on every run,
checks a separate scalar chemistry/geometry calculation, reconstructs discrete
patch selection and dominance ordering, and exercises the real JS exporter.
Small structures also receive an independent all-pairs sampled-SASA check.
Analytic examples and deliberately corrupted outputs test these verifiers.

The scalar oracle shares the declared chemistry/reference-ASA tables and the
PDB parser with the application. Large-structure SASA shares Biopython 1.85.
These dependencies are recorded and are not independent biological validation.
Frozen old benchmark aggregates remain historical data, not test oracles.
