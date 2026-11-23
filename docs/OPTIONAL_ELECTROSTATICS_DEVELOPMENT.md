# Optional physics-aware electrostatics — development module

> **Status:** Development/experimental. Not a replacement for any of the canonical 11 maps, not adsorption prediction, and not validated against an APBS executable in the current development environment.

## Purpose and separation

InterfaceScout's 11 existing chemical compatibility maps are calculated as before. The optional panel calls **PDB2PQR + PROPKA** to assign protonation/charges/radii on a protein-only input and can call **APBS** to solve the linearized Poisson–Boltzmann equation for an electrostatic potential grid. The result contains both the original and the versioned electrostatics-aware 11-map views, together with the raw DX/PQR grids and audit trail. Charge-mediated channels have constrained, explicit electrostatics-aware re-evaluation; neutral hydrophobic/aromatic chemistry is not changed by APBS.

### Dependencies

Install PDB2PQR (which includes the PROPKA integration) and APBS using their official packages into the same environment / system PATH. Do not list these external tools in core `backend/requirements.txt`, because the base application must remain fast and available without them. The status endpoint `GET /electrostatics/status` reports availability.

Example for a separate tested environment (adapt to your OS and exact tool installation):

```bash
python -m pip install pdb2pqr
# Install the APBS CLI from the official APBS distribution separately.
# Confirm both executables with:
pdb2pqr --help
apbs --version
```

On Windows run the base `run_local.bat` once, then `install_optional_electrostatics_windows.bat`; install the APBS executable separately and ensure `apbs` is on PATH before reopening the app. The optional UI is at the bottom of the sidebar; it reads pH, ionic strength and temperature from the solution inputs. Salt is NaCl or KCl, both modeled as monovalent ions. The ionic concentration for each +/-1 species is `ionic_mM / 1000` M. This is **not** a specific-ion binding model; NaCl and KCl yield the same continuum ionic strength description.

## Numerical scope and limits

- Input PDB: first model, selected protein chains, protein-only heavy atoms, using the existing prep policy.
- PROPKA: optional titration-state assignment at the requested pH via PDB2PQR. Disabling PROPKA means pH does **not** determine PDB2PQR titration.
- PDB2PQR: selected forcefield (PARSE, AMBER or CHARMM), atom partial charges and radii, and parameterized net molecular charge. Do not confuse net partial-charge sum with surface potential or net formal charge.
- APBS: linearized Poisson–Boltzmann (`lpbe`), 1:1 salt, solvent/protein dielectric, supplied temperature, multi-grid, native potential in `kBT/e`. Potential mV conversion factor is `kB T / e`. The grid is exported; no validated molecular surface interpolation or region-specific potential summary has yet been implemented.
- Protein-only structural policy **removes metal ions, bound ligands, glycans, post-translational modifications and other HETATM components**. The report warns if they were present. Their electrostatic effects are not modeled.
- APBS grid resolution is automatically selected subject to a safety/memory cap. Scientific production requires mesh convergence checks; the present prototype alone does **not** establish convergence or accuracy.\n- A separate 11-map view is generated only if both PDB2PQR/PROPKA and APBS succeed, and the mapped residue coverage checks pass. The earlier 11 maps remain accessible and exportable.
- Errors (missing executables, unsupported forcefield, invalid PQR, grid too large, nonzero exit) must stop the calculation rather than silently fallback to the canonical maps.

The returned archive contains:
```text
metadata/results.json
metadata/pdb2pqr.log
metadata/prepared_protein.pdb
data/prepared.pqr
metadata/apbs.in             # APBS enabled only
metadata/apbs.log            # APBS enabled only
data/potential.dx            # APBS enabled only\ndata/canonical_11maps.json   # APBS enabled only\ndata/physics_aware_11maps.json # APBS enabled only\ndata/charged_residue_surface.csv\ndata/standard_region_APBS_overlay.csv
```

## Verification before release

1. Install actual PDB2PQR/PROPKA and APBS binaries in CI and check command compatibility.
2. Run 1UBQ and 1AKI against known PDB2PQR/APBS reference cases and verify protonation, net charge and DX grid values.
3. Confirm pH dependence of PROPKA assignment, ionic strength and dielectric parameter transfer, temperature conversion.
4. Test missing HETATM structures and large biological units (including out-of-memory refusal).
5. Validate the implemented probe-expanded atom-sphere potential sampling and patch overlays against independent mesh-based surface controls; current sampling is not an exact SES/triangulated molecular surface.
6. Validate the separate 2.5.0-physics-context-prototype 11-map model on the frozen 48-PDB development panel and held-out 40-PDB evaluation panel. Compare original and modified ranks; do not conflate the 48 and 40 or claim experimental adsorption prediction.

References: https://pdb2pqr.readthedocs.io/en/latest/using/index.html and https://apbs.readthedocs.io/en/latest/using/
