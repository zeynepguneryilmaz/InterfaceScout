# InterfaceScout

InterfaceScout is an open-source, deterministic framework for locating **plausible protein-side contact regions for generalized surface chemistries** from a protein structure and solution conditions.

A single analysis prepares the protein once and generates every supported surface-chemistry map. Users can inspect individual maps or combine multiple chemistry classes as a derived overlay without rerunning the calculation.

InterfaceScout does **not** predict adsorption free energy, adsorption capacity, a unique adsorption orientation, or an atomically exact contact map.

## Workflow

For one protein structure, InterfaceScout:

1. calculates side-chain solvent accessibility with Shrake–Rupley SASA;
2. retains exposed residues using `scRSA >= 0.05`;
3. evaluates pH-dependent ionization-state availability where required by the interaction mechanism;
4. generates all canonical surface-chemistry maps from the same prepared structure;
5. aggregates compatible residues at **6 and 9 Å**;
6. constructs non-transitive **8 Å coarse interface patches** around local chemistry maxima on the same outward-facing protein surface;
7. ranks each single-channel map by Pareto dominance across chemistry support, accessibility, patch coherence, and orientation coherence.

There is no material catalogue, no material-specific fitted weight, and no GNM/RIN/APBS term in prediction or output.

## Surface-chemistry maps

One calculation generates:

- **Anionic surface**
- **Cationic surface**
- **Hydrophobic surface**
- **π / aromatic surface**
- **H-bond donor surface**
- **H-bond acceptor surface**
- **Metal-oxide surface**
- **Calcium/phosphate charged sites**
- **Transition-metal coordination**
- **Soft-metal sulfur affinity**
- **Phosphate-rich surface**

The chemistry class or classes should be selected according to experimentally supported surface composition, functional groups, charge state, coating, spectroscopy, or other interface characterization. A real interface may contain more than one chemistry; InterfaceScout therefore does not force the system into a named-material category.

## Combining multiple surface chemistries

After analysis, two or more chemistry classes can be selected to create a **Composite chemistry overlay**. This is deliberately not a fitted material score. For each residue, the overlay takes the **maximum normalized support** among the selected canonical maps, avoids double counting, and retains the dominant and supporting chemistry channels.

The composite avoids inventing unknown fractions such as `70% hydrophobic + 30% anionic`. If quantitative surface-composition fractions are independently known, they remain external experimental information.

**Single-channel maps are the canonical prediction outputs. The Composite chemistry overlay is a derived exploratory visualization, not a separately validated combined model or binding-affinity score.**

## How to interpret results

**Pareto front 1** contains the primary plausible interface candidates for a chemistry class. Display rank orders patches for inspection.

Each patch reports its center and member residues, chemistry support, mean accessibility, multiscale patch coherence, and orientation coherence. Patch membership is intentionally coarse: InterfaceScout localizes plausible protein-side surface regions rather than asserting that every member residue is an atomically exact adsorption contact.

Experimental comparison should match evidence resolution. Residue-resolved evidence can be compared by direct overlap and spatial proximity; peptide, helix, region, domain, or molecular-face evidence should be evaluated at the corresponding regional resolution. Scores and ranks are meaningful within a chemistry map and are not quantitative affinities between chemistry classes.

## Inputs and outputs

Inputs:
- PDB ID or local PDB file
- optional protein chain or comma-separated chains
- pH
- ionic strength
- temperature

Outputs:
- all canonical single-channel chemistry maps
- ranked coarse patches for every chemistry map
- residue-level chemistry support
- user-selected Composite chemistry overlay
- complete JSON export
- CSV export for the displayed chemistry or composite view

## Canonical settings

- SASA probe radius: **1.40 Å**
- Shrake–Rupley sampling: **200 points/atom**
- surface threshold: **scRSA >= 0.05**
- multiscale aggregation: **6 and 9 Å**
- coarse patch radius: **8 Å**
- empirical fitted weights: **none**

The 200-point SASA setting and 6/9 Å aggregation pair were selected on a separate adsorption-label-free development panel before external experimental evaluation. Experimental adsorption labels were not used to tune the prediction model.

## Installation

InterfaceScout supports **Python 3.10, 3.11, and 3.12**. The publication build is continuously smoke-tested on Linux, macOS, and Windows through GitHub Actions.

### Windows

Double-click:

```text
run_local.bat
```

The launcher validates Python, creates the local virtual environment when needed, installs the runtime dependencies, starts InterfaceScout, and opens the browser.

### macOS

First setup:

```bash
bash run_local.sh
```

Later runs can use:

```bash
bash start.command
```

### Linux

First setup:

```bash
bash run_local.sh
```

Later runs:

```bash
bash start.sh
```

The application is served locally at `http://localhost:8000`.

## Validation and reproducibility

The live application uses `backend/app.py`. Publication reproducibility resources are kept separately under `backend/v2/validation/`; the `v2` directory name is retained only as provenance of the frozen publication implementation and is not a user-selectable version.

Two independent workflows are retained:
- **parameter selection:** the adsorption-label-free eight-protein development panel;
- **prediction-first experimental evaluation:** saved predictions are generated without experimental localization labels, followed by a separate ground-truth comparison stage.

The validation-only material-to-chemistry mapping is not exposed by the public application and assigns no numerical material weights.

## Repository structure

```text
InterfaceScout/
├── backend/
│   ├── app.py
│   ├── main.py
│   ├── requirements.txt
│   ├── requirements-publication-lock.txt
│   └── v2/
│       ├── api.py
│       ├── interface_engine.py
│       ├── coarse_patch.py
│       ├── geometry.py
│       ├── model_settings.py
│       ├── prepare.py
│       └── validation/
├── frontend/
│   └── index.html
├── run_local.bat
├── run_local.sh
├── start.sh
├── start.command
└── README.md
```

## Scope

InterfaceScout addresses:

> Given a folded protein structure and solution condition, which exposed protein surface regions are plausible candidates for contact with one or more defined surface-chemistry classes?

It does not model material porosity, interfacial hydration, mass transfer, nanoparticle aggregation, adsorption-induced unfolding, equilibrium coverage, or quantitative binding affinity.

## Citation

If you use InterfaceScout in published work, please cite the associated InterfaceScout publication when publication details become available.
