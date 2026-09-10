# InterfaceScout

InterfaceScout is an open-source, deterministic framework for locating **plausible protein-side material-contact regions** from a protein structure, a generalized material-interface chemistry, and solution conditions.

It is intended for screening and interpretation of protein–material interactions. InterfaceScout does **not** predict adsorption free energy, adsorption capacity, a unique adsorption orientation, or an atomically exact contact map.

## What InterfaceScout does

For an input protein structure, the program:

1. calculates side-chain solvent accessibility using Shrake–Rupley SASA;
2. retains solvent-exposed residues using `scRSA >= 0.05`;
3. assigns residues compatible with the selected material-interface chemistry;
4. applies pH-dependent ionization-state availability where the interaction mechanism requires it;
5. aggregates compatible residues at **6 and 9 Å** spatial scales;
6. constructs non-transitive **8 Å coarse interface patches** around local chemistry maxima on the same outward-facing protein surface;
7. ranks patches by Pareto dominance across chemistry support, accessibility, patch coherence, and orientation coherence.

GNM and residue-interaction-network descriptors are structural context only and do not alter patch membership or ranking. Literature interaction-energy values and optional APBS electrostatic descriptors are also auxiliary only.

## Material-interface modes

The user selects the generalized mode that best represents the experimental surface:

- silica / deprotonated silanol-rich
- citrate-coated gold
- MPA-coated gold
- PAA-coated iron oxide
- polystyrene
- aluminum hydroxide
- iron oxide
- hydroxyapatite
- calcium fluoride
- generic hydrophobic
- generic anionic
- generic cationic

Each mode selects a predefined chemistry channel. No fitted material-specific numerical weight is used.

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
- material-interface mode
- pH
- ionic strength
- temperature

## Outputs

InterfaceScout reports ranked coarse interface patches. Each patch includes:

- patch center and member residues
- Pareto front and display rank
- chemistry support
- mean solvent accessibility
- multiscale patch coherence
- orientation coherence
- optional structural-network context

Results can be exported as JSON, CSV, or PDB for external visualization.

### How to interpret the output

**Pareto front 1** contains the primary plausible interface candidates. Display rank orders the reported patches for inspection.

Patch membership is intentionally coarse. Comparison with experiments should therefore match the resolution of the experimental evidence: use residue overlap for residue-resolved measurements and spatial/region agreement for region-level measurements.

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

> Given a folded protein structure, solution condition, and generalized material-interface chemistry, which exposed protein surface regions are plausible candidates for material contact?

It does not model material porosity, interfacial hydration, mass transfer, nanoparticle aggregation, adsorption-induced unfolding, equilibrium coverage, or quantitative binding affinity.

## Citation

If you use InterfaceScout in published work, please cite the associated InterfaceScout publication when the publication details become available.
