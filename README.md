# InterfaceScout

InterfaceScout is an open-source, deterministic framework for locating **plausible protein-side material-contact regions** from a protein structure and solution conditions.

A single analysis calculates the protein state once and then generates **all supported chemistry maps**. The user can switch among these maps instantly without rerunning the calculation.

InterfaceScout is intended for screening and interpretation of protein–material interactions. It does **not** predict adsorption free energy, adsorption capacity, a unique adsorption orientation, or an atomically exact contact map.

## Single-pass workflow

For one input protein structure, InterfaceScout:

1. calculates side-chain solvent accessibility using Shrake–Rupley SASA;
2. retains solvent-exposed residues using `scRSA >= 0.05`;
3. evaluates pH-dependent ionization-state availability;
4. computes the shared structural context once;
5. generates every supported chemistry map from the same prepared structure and solvent-exposure state;
6. aggregates compatible residues at **6 and 9 Å** spatial scales;
7. constructs non-transitive **8 Å coarse interface patches** around local chemistry maxima on the same outward-facing protein surface;
8. ranks patches by Pareto dominance across chemistry support, accessibility, patch coherence, and orientation coherence.

GNM and residue-interaction-network descriptors are structural context only and do not alter patch membership or ranking. Literature interaction-energy values and optional APBS electrostatic descriptors are auxiliary only.

## Chemistry maps

One analysis generates the following maps:

- anionic surface compatibility
- cationic surface compatibility
- hydrophobic surface compatibility
- π / carbon-like compatibility
- H-bond donor surface compatibility
- H-bond acceptor surface compatibility
- oxide surface compatibility
- hydroxyapatite / Ca²⁺ compatibility
- transition-metal coordination
- gold affinity
- phosphate surface compatibility

After the calculation finishes, these maps appear as selectable tabs in the interface. Switching maps changes the displayed 3D patches, patch cards, table, and selected-map export without recomputing the protein.

## Material presets

Material presets are provided only to help the user identify the most relevant chemistry map. They do not trigger separate calculations.

Examples include:

- silica / deprotonated silanol-rich → anionic map
- citrate-coated gold → anionic map
- MPA-coated gold → anionic map
- PAA-coated iron oxide → anionic map
- polystyrene → hydrophobic map
- aluminum hydroxide → cationic map
- iron oxide → oxide map
- hydroxyapatite → hydroxyapatite map
- calcium fluoride → cationic map

The selected material preset only determines which map is shown first. **All maps are still calculated in the same run.** No fitted material-specific numerical weight is used.

## Installation

### Windows

Double-click `run_local.bat`.

### macOS / Linux

First run:

```bash
bash run_local.sh
```

Later runs:

```bash
bash start.sh
```

On macOS, `start.command` can also be used.

The application opens locally at:

`http://localhost:8000`

## Inputs

- PDB ID or local PDB file
- optional protein chain
- pH
- ionic strength
- temperature
- optional material preset for choosing the initially displayed map

## Outputs

For every chemistry map, InterfaceScout reports ranked coarse interface patches. Each patch includes:

- patch center and member residues
- Pareto front and display rank
- chemistry support
- mean solvent accessibility
- multiscale patch coherence
- orientation coherence
- optional structural-network context

The complete run can be exported as JSON. The currently selected chemistry map can also be exported as CSV.

### How to interpret the output

**Pareto front 1** contains the primary plausible interface candidates for the selected chemistry map. Display rank orders the reported patches for inspection.

Patch membership is intentionally coarse. Comparison with experiments should therefore match the resolution of the experimental evidence: use residue overlap for residue-resolved measurements and spatial/region agreement for region-level measurements.

Scores and ranks are meaningful **within a chemistry map**. They should not be interpreted as quantitative binding affinities between different material chemistries.

## Canonical settings

- SASA probe radius: **1.40 Å**
- Shrake–Rupley sampling: **200 points/atom**
- surface threshold: **scRSA >= 0.05**
- multiscale aggregation: **6 and 9 Å**
- coarse patch radius: **8 Å**
- empirical fitted weights: **none**

The 6/9 Å aggregation pair was selected on an adsorption-label-free development panel before external experimental evaluation. Experimental adsorption labels were not used to tune the model.

## Validation and reproducibility

Publication-related sensitivity and external-validation resources are stored under:

`backend/v2/validation/`

The folder name reflects development history only. There is one public program: **InterfaceScout**.

The development/sensitivity panel and the external experimental-validation panel are separate. Prediction-first records are retained so that experimental interface labels are not used to construct or rank patches.

## Repository structure

```text
InterfaceScout/
├── backend/
│   ├── app.py               single application entry point
│   ├── main.py              internal chemistry/SASA kernel
│   ├── requirements.txt
│   └── v2/                  internal model modules and validation resources
├── frontend/
│   ├── index.html
│   ├── logo.png
│   └── favicon.png
├── run_local.bat
├── run_local.sh
├── start.sh
├── start.command
├── LICENSE
└── README.md
```

Users do not select between versions or engines. `backend/app.py` starts the canonical InterfaceScout application.

## Scope

InterfaceScout addresses the question:

> Given a folded protein structure and solution condition, which exposed protein surface regions are plausible candidates for contact under each supported material-interface chemistry?

It does not model material porosity, interfacial hydration, mass transfer, nanoparticle aggregation, adsorption-induced unfolding, equilibrium coverage, or quantitative binding affinity.

## Citation

If you use InterfaceScout in published work, please cite the associated InterfaceScout publication when the publication details become available.
