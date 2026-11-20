# InterfaceScout v1.1.0-rc1 | candidate change log

**Prior publication release:** v1.0.0 (2026-10-05), unchanged. `main` is not overwritten.

## Technical updates (not a new material-adsorption model)

- Identity-aware residue matching across renamed/renumbered PDB chains in same-protein comparison and multi-structure consensus.
- Coordinate-based tie-breaking for map/patch geometry; author identifiers are not numerical ranking features.
- Strict PDB input checks, model and altloc reporting, missing atom warnings, explicit mmCIF rejection.
- Floating-point precision retained until result presentation or legacy PDB B-factor formatting.
- B-factor exports now operate on the exact prepared first structural model and selected chains, retaining the originally supplied PDB separately.
- 21-condition / 231-map offline numerical reproduction, independently written scalar checks and two independent SASA spot checks.
- Archived cross-platform real-browser exports, checksum manifests and strict publication package assembly.

**Unchanged:** 11 canonical map definitions and chemistry tables; 1.40 Å probe; 200 points per atom; scRSA 0.05; 6/9 Å multiscale radii; 8 Å coarse patch; descriptive/non-predictive scope.

## Compatibility note

Exact legacy rankings and previously frozen numerical output files are historical artifacts. Low-rank near-ties and scores affected by prior intermediate rounding may change. Do not compare arbitrary source snapshots without recording their Git hashes and runtime environments.
