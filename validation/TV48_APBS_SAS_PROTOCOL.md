# TV48: real APBS potential on a solvent-accessible sampling shell

This study is a **post-hoc technical sensitivity analysis**, not an adsorption-site
accuracy experiment. Its input is the accepted, SHA-256-frozen TV48 manifest.

## Motivation

The existing optional APBS layer samples the potential at each PQR heavy
side-chain atom and averages by residue. Those coordinates can lie inside the
low-dielectric protein. This is not equivalent to the electrostatic potential
on a protein's solvent-facing boundary. Very large site values must therefore
not be interpreted as material binding energies or adsorption probabilities.

## What is actually computed

The separate audit `tv48_apbs_sas_grid_audit.py`:

1. Prepares the same frozen biological assembly and checks its assembly hash.
2. Uses real PDB2PQR with AMBER charges and PROPKA titration at the specified pH.
3. Runs real APBS LPBE at 298 K, protein dielectric 2 and solvent dielectric
   78.54, with the documented ionic conditions and auto-focused grid.
4. Uses a deterministic 42-point Fibonacci sphere per PQR heavy atom,
   expanded by a 1.4 Å probe radius and an additional 0.25 Å **outward
   displacement**. Points overlapping another expanded heavy-atom sphere
   are discarded. This approximates a **solvent-accessible probe-center
   shell**, **not a triangulated solvent-excluded surface** (SES).
5. Interpolates the APBS OpenDX grid at valid shell points; calculates
   residue-level surface-area-weighted means. Glycine can have a valid
   shell potential despite lacking a heavy side chain.
6. Compares the shell potential with the previous within-atom side-chain
   estimate using paired residue IDs, without changing any of the canonical
   11 maps or the Pareto ranking.
7. Samples 48/48 frozen TV48 entries at baseline pH 7.4, ionic strength
   150 mM, and 298 K. Separately, eight proteins, selected by
   **predeclared structure-size quantiles**, are checked at nominal grid
   target 0.6, 0.8, 1.1 Å; ionic 50/150/300 mM; and pH 6.5/7.4/8.5.
   Effective grid spacings and actual dimensions are recorded. The
   APBS grid dimension cap of 257 must not be confused with demonstrated
   grid convergence: a capped refinement is **not a finer-grid test**.

The external-solver artifacts, real per-residue results and representative
figures are emitted only by a successful run of
[`tv48-apbs-sas-grid.yml`](../.github/workflows/tv48-apbs-sas-grid.yml).
A GitHub workflow in progress is **not** equivalent to a PASS or final results.

## Interpretation gates

- No adsorption free energy, contact geometry, or measured binding site follows
  from the potential values.
- Only real, non-interpolated-to-target-number grid dimensions and the
  **actual** effective spacing can support numerical convergence.
- A finite 42-direction shell is an approximation. A subsequent repeat with
  84/162 directions would measure surface-sampling discretization sensitivity.
- The same 11 maps remain the canonical output. APBS is a separately flagged
  physical annotation and must not be silently substituted for a map score.
- Correlation of two derived descriptors is *not* proof that APBS provides
  experimentally useful extra information. Such claims require an independent
  evidence endpoint and matched controls.
- TV48 remains technical validation. The independent scientific evidence
  panel is the separately audited active ISV candidate set.
