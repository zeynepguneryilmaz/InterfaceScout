# Verified real APBS/PDB2PQR/PROPKA integration — 1UBQ

**Result: PASS (24 November 2026, TRT).** This report covers the **optional, physically grounded electrostatics layer only**, not an APBS-informed reparameterization of the eleven canonical InterfaceScout chemistry maps.

## Exact execution record

- GitHub Actions [real-solvers run 45290088566](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45290088566), job `real-solvers`, conclusion **success**.
- Run commit: `11b0e180ac5675682c9bf7559de1199bda13201f` (subsequently merged into `main` through [PR #2](https://github.com/zeynepguneryilmaz/InterfaceScout/pull/2)).
- Runner: Ubuntu 24.04, Python 3.11.
- External solver packages: APBS 3.4.1 (Ubuntu package `3.4.1-5.1ubuntu2`), PDB2PQR 3.7.1, PROPKA 3.5.1.
- Structure: RCSB **1UBQ**, chain A, retrieved at test execution; pH **7.4**, ionic strength **150 mM**, temperature **298 K**.
- Force field and model: AMBER PQR parameters; PROPKA titration states via PDB2PQR; APBS **linearized** Poisson–Boltzmann (`lpbe`), protein dielectric 2, solvent dielectric 78.54, symmetric 1:1 electrolyte (0.150 M), solvent probe 1.4 Å.
- Atomic potential was read from a *real APBS OpenDX grid* and interpolated at side-chain heavy atoms. Values are in **kT/e** and are descriptive, not adsorption free energies.

## Verified assertions

| Test | Result |
| --- | --- |
| Install and invoke actual APBS, PDB2PQR and PROPKA | PASS |
| Five unit/regression tests (including zero-radius hydrogen handling) | **5/5 PASS** |
| Full 1UBQ external-solver pipeline, OpenDX parsing, residue mapping | **PASS** |
| PQR atoms parsed | **1,231** |
| PDB/PQR matched residues | **76/76** |
| InterfaceScout public maps present | **11/11** |
| Before/after canonical map objects exactly equal | **PASS** |
| Before/after residue profiles exactly equal | **PASS** |
| At least one finite, populated side-chain potential readout | **PASS** |
| JSON evidence uploaded as GitHub Actions artifact | **PASS** |

Artifact name: `InterfaceScout-1UBQ-APBS-PROPKA-evidence` (run 45290088566; artifact ID **14872858660**). It contains the per-residue electrostatics results, input provenance, and unchanged-canonical-model assertion. GitHub Actions artifacts expire; permanently archive a copy for publication.

## Important scope limitations

This is a **real numerical end-to-end integration test on one protein**. It is not an external experimental validation of electrostatic predictions, nor a 48-protein APBS stability screen. No APBS-derived weights, fitting, or automatic reranking of the existing 11-map/Pareto model has been introduced. PROPKA changes protonation assignment for the PQR used by APBS, but the optional layer does not replace the existing canonical pKa rules inside the original 11 maps.

Next independent checks: realistic pH/ionic-strength perturbation tests; sensitivity of the APBS grid; larger/heterogeneous biological assemblies; TV48 full electrostatic screen; and cross-platform solver availability. Do not label those checks completed based solely on this 1UBQ pass.
