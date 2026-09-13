# InterfaceScout

InterfaceScout is an open-source framework for prioritizing plausible **protein-side contact regions at material interfaces** from a folded protein structure, solution pH, and an interface-chemistry hypothesis.

A single analysis prepares the protein once and generates 11 canonical surface-chemistry maps. Each map combines side-chain solvent exposure, chemistry compatibility, pH-dependent state availability, multiscale spatial enrichment, coarse patch construction, and Pareto prioritization. The output is a short, interpretable set of candidate protein-surface regions that can guide residue-resolved experiments or focus subsequent material-specific simulations.

## Canonical surface-chemistry maps

One analysis generates:

- Anionic surface
- Cationic surface
- Hydrophobic surface
- π / aromatic surface
- H-bond donor surface
- H-bond acceptor surface
- Metal-oxide surface
- Calcium/phosphate charged sites
- Transition-metal coordination
- Soft-metal sulfur affinity
- Phosphate-rich surface

The chemistry map is selected from independent knowledge of the interface, such as surface composition, coating chemistry, charge state, spectroscopy, or other material characterization. Several maps can be inspected side by side when more than one interaction motif is plausible.

## What the outputs mean

For each chemistry map, InterfaceScout reports residue-level values and candidate patches.

- **Chemistry support**: fraction of patch members compatible with the selected chemistry class.
- **Mean accessibility**: average side-chain relative solvent accessibility of patch members.
- **Patch coherence**: persistence of the local chemistry-compatible signal across the 6 and 9 Å aggregation scales.
- **Orientation coherence**: degree to which patch members occupy a common coarse outward-facing protein surface.

Pareto fronts define multiobjective priority classes. A deterministic within-front order provides a reproducible sequence for inspecting candidates.

The optional **Composite chemistry overlay** summarizes residue-level support across two or more selected canonical maps. It is intended for qualitative inspection of mixed chemistry hypotheses; quantitative propensity values are interpreted within each individual map because every map is normalized independently.

## Inputs

- PDB ID or local PDB file
- optional chain or comma-separated chains
- pH
- ionic strength
- temperature

In the publication model, pH enters the scoring calculation. Ionic strength and temperature are recorded with the run for reproducibility.

## Outputs

The browser application provides:

- all 11 canonical chemistry maps;
- residue-level chemistry values;
- all candidate patches and Pareto fronts;
- an interactive 3D structure viewer;
- a complete Excel workbook;
- B-factor PDB property tracks for canonical maps and the optional Composite overlay.

The Excel workbook contains `Run_Method`, `Map_Definitions`, `Residue_Values`, `Patches`, and `Patch_Members`; `Composite_Selected` is added when a Composite overlay is active.

For a canonical map, the downloadable PDB writes normalized residue propensity (0–100) into the standard B-factor field. `REMARK 900` records identify these values as InterfaceScout scores rather than experimental crystallographic B-factors.

## Publication settings

The frozen publication settings are:

- Shrake–Rupley probe radius: **1.40 Å**
- Shrake–Rupley sampling: **200 points/atom**
- surface threshold: **scRSA ≥ 0.05**
- multiscale aggregation radii: **6 and 9 Å**
- coarse candidate-patch radius: **8 Å**
- fitted material-specific weights: **none**

The 200-point sampling density and 6/9 Å aggregation pair were selected on an independent eight-protein development panel using numerical convergence, spatial consensus, and locality. Experimental adsorption-localization labels were not used for parameter selection.

## Installation

InterfaceScout supports Python **3.10, 3.11, and 3.12**.

### Windows

Run:

```text
run_local.bat
```

### macOS / Linux

First setup:

```bash
bash run_local.sh
```

Later launches can use `start.command` on macOS or `start.sh` on Linux.

The local application is served at `http://localhost:8000`.

## Publication reproducibility

The repository contains a single frozen publication workflow under `backend/v2/validation/` and a generated `publication_data/` package.

Two stages are kept separate:

1. **Development / parameter robustness** — eight structurally diverse proteins are used only for numerical and spatial setting selection.
2. **Literature-derived benchmark** — predictions are saved first, then compared with independently encoded experimental protein-side localization.

The final benchmark contains seven primary conditions across six proteins and two secondary challenge cases. Experimental annotations are stored at the resolution supported by each source: residue anchor, localized residue set, peptide/segment, or broader region.

To reproduce the complete publication package locally:

```bash
cd backend
python -m v2.validation.reproduce_publication
```

This regenerates the development outputs, benchmark predictions, exact prepared PDB snapshots, comparison tables, and verification summary in `publication_data/`.

`publication_data/verification.json` checks the principal manuscript-level reproducibility targets, including the selected 200-point SASA setting, the 6/9 Å aggregation pair, seven primary Top-3 near-8 Å hits, four primary Top-5 direct-overlap hits, and median primary Top-5 near-8 Å recall of 0.80.

See `backend/v2/validation/README.md` and `backend/v2/validation/PROTOCOL.md` for the frozen workflow and interpretation rules.

## Repository structure

```text
InterfaceScout/
├── backend/
│   ├── app.py
│   ├── main.py
│   ├── requirements.txt
│   ├── requirements-publication-lock.txt
│   └── v2/
│       ├── interface_engine.py
│       ├── coarse_patch.py
│       ├── geometry.py
│       ├── model_settings.py
│       ├── prepare.py
│       ├── surface_modes.py
│       └── validation/
├── frontend/
├── publication_data/
├── run_local.bat
├── run_local.sh
├── start.command
├── start.sh
├── LICENSE
└── README.md
```

## Scope

InterfaceScout addresses the practical question:

> Given a folded protein structure, solution pH, and a defined interface-chemistry hypothesis, which exposed protein surface regions should be prioritized for experimental or higher-resolution follow-up?

It is a chemistry-aware screening and experiment-design layer. Material-specific free-energy calculations, explicit interfacial solvent, ligand-density effects, nanoscale curvature, and adsorption-induced structural adaptation can be examined in subsequent specialized simulations or experiments when needed.

## License

InterfaceScout is released under the MIT License. See `LICENSE`.

## Citation

Please cite the associated InterfaceScout publication once bibliographic details are available.
