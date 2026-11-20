# Offline numerical reproducibility

Install `backend/requirements.txt` in Python 3.10–3.12 and provide Node 22.
From the repository root, no network is used by either command:

```sh
python validation/reproduce.py --suite quick --output reproducibility-results/quick
python validation/reproduce.py --suite full --output reproducibility-results/full
```

Quick: four conditions / 44 maps, including raw RCSB 1CRN and 1SHG snapshots,
and 1CRN at pH 0 and 14. Full: 21 conditions / 231 maps, including the original
19-condition / 209-map P0 comparison. The older 19 inputs are versioned example
coordinate PDBs with annotated B-factors; these annotations are **not** used as
reference scores or inputs to scoring. Two added raw RCSB files retain their
original records. Source provenance, conditions and byte SHA-256 values are in
`inputs.json`. Inputs are never downloaded or refreshed during verification.
PDB byte content is protected from Git line-ending conversion.
Manifest SHA-256 values identify exact committed LF bytes. The additional
`p0_capture_sha256` fields identify the historical Windows CRLF snapshots;
the causal audit explicitly proves that their difference is line endings only.

Each run checks the input and specification hashes, recalculates SASA and all
11 channels, writes unrounded result JSON, checks a separate scalar chemistry
and geometry implementation, reconstructs patch selection/membership and Pareto
depth via a dominance DAG, and executes the real frontend B-factor exporter.
The report records Python, installed packages, Node, platform, source hashes,
Git state, input/output hashes, model version, parameters and numerical errors.
Any failed check gives a nonzero exit code. Existing output files in the selected
output directory are overwritten. `verification.json` is the run manifest;
generated score files are evidence, not automatically accepted golden fixtures.

The acceptance contract in `TOLERANCES.md` was committed before implementation
and execution. Exact decisions use verified unrounded application descriptors;
independent scalar summation can differ at machine precision. Reports include
minimum positive descriptor gaps and counts of adjacent gaps within tolerance;
these are sensitivity diagnostics, never a rule for changing ties or rankings.

Independence has explicit limits: the scalar oracle shares the published chemical
tables, reference ASA values and PDB parser. Two small raw structures additionally
use a separate all-pairs point-occlusion SASA implementation with the same declared
200-point quadrature. Large structures share pinned Biopython SASA. Hand-derived
examples and deliberate mutations test the verifier itself. This does not prove
biological adsorption accuracy, quadrature convergence, or cross-platform bitwise
identity. The historical aggregate checker `validation/verify_publication.py`
continues to check frozen report consistency and has a different purpose.

Optional historical causal audit (requires full Git history and original P0
capture directories from the P0 evidence package):

```sh
python validation/reproducibility/audit_rounding.py --before numeric-before --after numeric-after --recomputed reproducibility-results/full --output rounding-audit.json
```

This verifies that the historical kernel differs only by 20 rounding calls
(excluding its module docstring), replays it against identical checked inputs,
compares all 209 maps exactly to the archives, and compares current recomputation
to the P0 captures. It reports changed ranks, score descriptors, dominance margins,
members and seeds for the three changed maps. Archives are causal comparison
evidence; numerical correctness is established separately by `reproduce.py`.
