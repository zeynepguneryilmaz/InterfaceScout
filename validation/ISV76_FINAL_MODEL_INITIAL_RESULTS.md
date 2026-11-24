# ISV76 fixed-parameter model: initial full-panel numerical readout

Date: 2026-11-24. This is a **descriptive structural candidate panel**, not
76 experimentally validated protein/material adsorption interfaces.

## Frozen protocol
- Active list: `ISV_active_candidates_76_v3.csv` (76 PDBs; historical 80 is excluded from current aggregates).
- The TV48-frozen engine produced 11 canonical maps per structure and
  optional real PDB2PQR/PROPKA/APBS with area-weighted SAS-shell sampling.
- Conditions: pH 7.4; ionic strength 150 mM; 298 K; APBS 0.8 Å target grid;
  42 deterministic SAS directions/heavy atom, 1.4 Å probe + 0.25 Å offset.
- No change to chemistry-map score, patch algorithm, Pareto objectives or
  electrostatic parameters based on these results.

## Actual verified run
- [ISV76 first full scientific run](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45385665397), all 76 structures accounted for.
- 76/76 canonical 11-map analyses: **836** map instances.
- 64,592 patches over 11 maps, 51,931 observed standard residues and
  36,273 tRSA-exposed residues.
- 75/76 actual production APBS runs and SAS outputs: **34,853** exposed
  residues with SAS potentials.
- Median surface coverage across 75 successful APBS structures: **0.97868**.
- Nine map instances had no patches (all on `soft_metal_sulfur`).
  This is not an error and should not be remedied by inventing a patch.

**Overall workflow status is deliberately INCOMPLETE**, because the
PDB2PQR solver cannot repair residue **SER A 11 in 1GCT**: essential backbone
coordinates are absent. This is an input-structure integrity problem, not a
reason to omit 1GCT or claim 76/76 APBS. Its 11-map results remain available.

## Outliers and scope restrictions

**8C50** (non-globular curli-like fibril) has potential values far outside
the usual ISV76 range: minimum approx **-3331.67 kT/e**, median
approx **-380.69 kT/e**, 683 mapped exposed residues, and total PQR
charge approx **-219 e**. These values demand specific geometry,
charge, grid and solver plausibility checks. Until then 8C50 must be
flagged and excluded from ordinary pooled APBS magnitude summaries.
The numerical execution itself succeeded: this caution is about physical
interpretation and validity, not a fabricated execution failure.

**1HC1** hemocyanin (six chains, 3,804 residues) completed real APBS;
large-assembly and cofactor/state caveats persist. **1STP** remains a
biotin-bound-derived, ligand-removed tetramer, not experimentally verified
apo streptavidin.

**4,730** source residues have at least one incomplete expected sidechain,
requiring localized examination for map-specific claims. Structural
suitability and experimental localization evidence have not been accepted
automatically.

## Aggregate-CSV string-preservation fix

The initial aggregate pandas read inferred the PDB identifier `1E05`
as floating-point scientific notation `100000.0`. All per-protein raw
CSVs and JSON preserved the correct string. The aggregation implementation
now explicitly reads PDB IDs as strings, with inclusion-list membership
assertions. The corrected conversation-side archive regenerates all three
large combined CSVs from their exact original per-protein CSVs without
altering electrostatic or map values.

## Scientific interpretation

TV48 remains the development/technical verification set. ISV76 is
an untouched-parameter generalization/application dataset; any
future experiment-to-structure residue localization comparisons must
use original publication labels and evidence tiers with appropriate
construct, chain, surface-context and family matching.

The 76 PDBs are **76 technical candidates**, not 76 independent
experimentally demonstrated adsorption interfaces.
