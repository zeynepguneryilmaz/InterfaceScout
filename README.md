# InterfaceScout

InterfaceScout is an open-source, deterministic framework for locating **plausible protein-side contact regions for generalized surface chemistries** from a protein structure and solution conditions.

A single analysis prepares the protein once and generates every supported surface-chemistry map. Users can inspect individual chemistry maps, combine multiple chemistry classes as a derived overlay, and view GNM/RIN structural-context maps without rerunning the calculation.

InterfaceScout does **not** predict adsorption free energy, adsorption capacity, a unique adsorption orientation, or an atomically exact contact map.

## Workflow

For one protein structure, InterfaceScout:

1. calculates side-chain solvent accessibility with Shrake–Rupley SASA;
2. retains exposed residues using `scRSA >= 0.05`;
3. evaluates pH-dependent ionization-state availability where required by the interaction mechanism;
4. generates all canonical surface-chemistry maps from the same prepared structure;
5. aggregates compatible residues at **6 and 9 Å**;
6. constructs non-transitive **8 Å coarse interface patches** around local chemistry maxima on the same outward-facing protein surface;
7. ranks each single-channel map by Pareto dominance across chemistry support, accessibility, patch coherence, and orientation coherence;
8. calculates GNM and residue-interaction-network descriptors as separate structural context.

There is no material catalogue and no material-specific fitted weight.

## Surface-chemistry maps

One calculation generates these generalized chemistry classes:

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

A user should select and interpret the chemistry class or classes that are justified by the known surface composition, functional groups, charge state, coating, spectroscopy, or other experimental characterization of the interface under study.

A real interface may contain more than one chemistry. For example, a heterogeneous or functionalized surface can simultaneously expose charged, hydrophobic, hydrogen-bonding, aromatic, or coordinating groups. InterfaceScout therefore does not force a surface into one material label.

## Combining multiple surface chemistries

After analysis, select two or more chemistry classes to create a **Composite chemistry overlay**.

The composite is deliberately not a fitted material score. For each residue, the default overlay takes the **maximum normalized support** among the selected canonical chemistry maps. A residue that is supported by more than one chemistry is not counted repeatedly. The interface also preserves:

- the number of selected chemistry channels supporting the residue;
- the dominant chemistry channel at that residue;
- the underlying individual chemistry maps.

This avoids inventing unknown fractions such as `70% hydrophobic + 30% anionic`. If quantitative surface-composition fractions are independently known, they should be treated as external experimental information rather than silently inferred by InterfaceScout.

**Important:** single-channel chemistry maps are the canonical prediction outputs. The Composite chemistry overlay is a derived exploratory visualization and is **not a separately validated combined model or binding-affinity score**.

## Structural-context maps

The same run also provides separate maps for:

- **GNM fluctuation** — normalized native-state Cα mobility;
- **RIN degree** — local connectivity in the 4.5 Å heavy-atom residue interaction network;
- **RIN betweenness** — participation in shortest network paths;
- **RIN closeness** — network proximity to the rest of the protein.

These maps help interpret whether a predicted surface region is relatively mobile, rigid, locally connected, or network-central. **GNM and RIN do not alter patch membership or Pareto ranking.**

## How to interpret a chemistry-map result

**Pareto front 1** contains the primary plausible interface candidates for that chemistry class. Display rank orders patches for inspection.

Each patch reports its center and member residues together with chemistry support, mean accessibility, multiscale patch coherence, and orientation coherence. Patch membership is intentionally coarse: InterfaceScout is intended to localize a plausible protein-side surface region, not claim every member residue as an atomically exact adsorption contact.

Comparison with experiments should match the experimental resolution. Residue-resolved experiments can be compared by direct residue overlap and spatial proximity; peptide/helix/domain-level evidence should be compared at the corresponding regional resolution. Scores and ranks are meaningful within a chemistry map and should not be interpreted as quantitative affinities between different chemistry classes.

## Inputs

- PDB ID or local PDB file
- optional protein chain
- pH
- ionic strength
- temperature

No material name is required.

## Outputs

- all canonical single-channel chemistry maps
- ranked coarse patches for each chemistry map
- residue-level chemistry support used for map inspection/composite overlays
- Composite chemistry overlay for any user-selected set of two or more chemistry classes
- GNM fluctuation map
- RIN degree, betweenness, and closeness maps
- complete JSON export
- CSV export for the currently displayed chemistry, composite, or structural-context view

## Canonical settings

- SASA probe radius: **1.40 Å**
- Shrake–Rupley sampling: **200 points/atom**
- surface threshold: **scRSA >= 0.05**
- multiscale aggregation: **6 and 9 Å**
- coarse patch radius: **8 Å**
- empirical fitted weights: **none**

The 6/9 Å pair was selected on a separate adsorption-label-free development panel before external experimental evaluation. Experimental adsorption labels were not used to tune the model.

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

The application opens locally at `http://localhost:8000`.

## Validation and reproducibility

Publication-related sensitivity and external-validation resources are stored under `backend/v2/validation/`. The directory name reflects development history only; there is one public application, **InterfaceScout**.

The development/sensitivity panel and external experimental-validation panel are separate. Prediction-first records are retained so experimental interface labels are not used to construct or rank patches.

## Scope

InterfaceScout addresses the question:

> Given a folded protein structure and solution condition, which exposed protein surface regions are plausible candidates for contact with one or more defined surface-chemistry classes?

It does not model material porosity, interfacial hydration, mass transfer, nanoparticle aggregation, adsorption-induced unfolding, equilibrium coverage, or quantitative binding affinity.

## Citation

If you use InterfaceScout in published work, please cite the associated InterfaceScout publication when the publication details become available.
