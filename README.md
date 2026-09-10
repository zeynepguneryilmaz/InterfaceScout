# InterfaceScout

InterfaceScout is an open-source, deterministic framework for locating **plausible protein-side material-contact regions** from a protein structure, a generalized material surface chemistry, and solution conditions.

It is designed as a lightweight screening and interpretation tool for protein–material studies. InterfaceScout does **not** predict adsorption free energy, adsorption capacity, a unique orientation, or an atomically exact contact map.

## What the program does

For an input protein structure, InterfaceScout:

1. calculates side-chain solvent accessibility using Shrake–Rupley SASA;
2. identifies solvent-exposed residues (`scRSA >= 0.05`);
3. assigns chemistry-compatible residue classes for the selected material interface;
4. applies pH-dependent ionization-state availability where required;
5. aggregates compatible residues at **6 and 9 Å** spatial scales;
6. constructs non-transitive **8 Å coarse interface patches** around local chemistry maxima;
7. ranks patches by Pareto dominance across chemistry support, accessibility, patch coherence, and orientation coherence.

GNM and residue-interaction-network descriptors are reported only as structural context. They do not change patch membership or ranking. Literature interaction-energy values and optional APBS electrostatics are also auxiliary descriptors only.

## Supported material-interface modes

The user selects the generalized interface chemistry that best represents the experimental material:

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

These modes select a predefined chemistry channel; no fitted material-specific numeric weight is applied.

## Installation and use

### Windows

Double-click:

`run_local.bat`

### macOS / Linux

First run:

```bash
bash run_local.sh
```

Later runs:

```bash
bash start.sh
```

or, on macOS:

```bash
./start.command
```

The application opens locally at:

`http://localhost:8000`

## Inputs

The interface accepts:

- a PDB ID or local PDB file;
- an optional chain selection;
- material-interface mode;
- pH;
- ionic strength;
- temperature.

## Outputs

InterfaceScout reports ranked coarse interface patches. For each patch, the output includes:

- patch center and member residues;
- Pareto front and display rank;
- chemistry support;
- mean solvent accessibility;
- multiscale patch coherence;
- orientation coherence;
- optional structural-network context.

Results can be exported as JSON, CSV, or a PDB file in which patch ranking is encoded in the B-factor field for visualization.

### How to interpret the ranking

**Pareto front 1** identifies primary plausible interface candidates. Lower display rank indicates a stronger ordering within the reported patch list.

Patch membership is intentionally coarse. Agreement with experimental data should therefore be evaluated at the same resolution as the experiment: exact residue overlap for residue-resolved measurements, and spatial/region overlap for regional or orientation-level measurements.

## Canonical model settings

- SASA probe radius: **1.40 Å**
- Shrake–Rupley sampling: **200 points/atom**
- surface threshold: **scRSA >= 0.05**
- multiscale aggregation: **6 and 9 Å**
- coarse patch radius: **8 Å**
- weighted empirical fitting: **none**

The 6/9 Å aggregation pair was selected on an adsorption-label-free development panel before external experimental evaluation. Experimental adsorption labels were not used to tune the model.

## Validation and reproducibility

Scripts and machine-readable files used for parameter sensitivity and external experimental evaluation are kept under:

`backend/v2/validation/`

The development/sensitivity panel and the external experimental-validation panel are separate. The repository preserves prediction-first evaluation records so that experimental interface labels are not used to construct or rank patches.

## Repository structure

```text
InterfaceScout/
├── backend/
│   ├── main.py              internal chemistry/SASA kernel used by InterfaceScout
│   ├── requirements.txt
│   └── v2/                  canonical InterfaceScout implementation
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

The `v2` directory name is retained only as an internal development path. The public software is a single application named **InterfaceScout**; users do not select between software versions.

## Scope

InterfaceScout answers the question:

> Given a folded protein structure, solution condition, and generalized interface chemistry, which exposed protein surface regions are plausible candidates for material contact?

It does not model material porosity, interfacial hydration, mass transfer, nanoparticle aggregation, adsorption-induced unfolding, equilibrium coverage, or quantitative affinity.

## Citation

If you use InterfaceScout in published work, please cite the associated InterfaceScout publication after publication details become available.
