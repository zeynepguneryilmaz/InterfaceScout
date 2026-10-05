# InterfaceScout

[![Cross-platform CI](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/workflows/publication-ci.yml/badge.svg)](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/workflows/publication-ci.yml)

InterfaceScout is an open-source tool for identifying and prioritizing plausible **protein-side contact regions** at material interfaces from a protein structure and solution conditions.

The program prepares the selected protein structure once and calculates 11 canonical surface-chemistry maps. Candidate surface patches are prioritized using chemistry support, side-chain accessibility, patch coherence, and coarse orientation coherence without material-specific fitted weights.

## Publication model

The frozen publication model uses:

- Shrake-Rupley probe radius: **1.40 Å**
- Shrake-Rupley sampling: **200 points/atom**
- exposed-side-chain threshold: **scRSA ≥ 0.05**
- multiscale aggregation radii: **6 and 9 Å**
- coarse patch radius: **8 Å**

The origin of these values is deliberately separated:

- **200 points/atom** was the lowest tested sampling level that met the predefined numerical-convergence criteria against a 1000-points/atom reference on an eight-protein, adsorption-label-free development panel.
- **6/9 Å** ranked first among ten candidate radius pairs in the label-free consensus analysis. The earlier 5/8 Å pair ranked second and is not used by the publication model.
- **8 Å** is a predefined coarse geometric contact-region scale. It was not optimized against the literature adsorption benchmark.

The frozen selection evidence is retained in `validation/parameter_selection.json`.

## Surface-chemistry maps

InterfaceScout calculates anionic, cationic, hydrophobic, π/aromatic, H-bond donor, H-bond acceptor, metal-oxide, calcium/phosphate charged-site, transition-metal coordination, soft-metal sulfur-affinity, and phosphate-rich surface maps in one run.

## Inputs

- PDB ID or local PDB file
- optional chain selection
- pH
- ionic strength
- temperature

In the current model, **pH enters the scoring calculation**. Ionic strength and temperature are retained as run metadata and should not be interpreted as explicit electrostatic-screening or thermal-response calculations.

## Intended scope and limitations

InterfaceScout is intended for coarse **contact-region prioritization** when a folded protein structure is an informative representation of the adsorbing state and the material interface can be related to one or more canonical chemistry classes.

The current validation most directly supports folded globular proteins and folded oligomeric proteins under conditions where large adsorption-induced structural rearrangements are not dominant. Use more cautiously for intrinsically disordered proteins, strongly unfolding/restructuring adsorption regimes, membrane-embedded systems, covalent attachment, or systems whose relevant chemistry depends critically on cofactors, glycans, lipids, non-standard residues, or other hetero components removed by the standard preparation policy.

InterfaceScout does **not** predict absolute adsorption free energy, adsorbed amount, atomistically unique binding orientation, adsorption kinetics, or adsorption-induced unfolding. Ranked patches are plausible protein-side regions for experimental follow-up or higher-resolution simulation, not exact atomistic binding sites.

The literature benchmark contains seven primary conditions and two secondary challenge cases. The transferrin-polystyrene and catalase-polystyrene cases are retained specifically to expose the limitation of a static folded-state model when adsorption-associated structural change is important. See `validation/README.md`.

## Outputs

1. **Complete Excel workbook** with run/method information, map definitions, residue-level values, ranked patches, and patch members for all 11 maps.
2. **B-factor PDB** for the displayed map, with the 0-100 InterfaceScout residue propensity written to the standard B-factor field.

## Installation and use

InterfaceScout supports **Python 3.10-3.12**. The pinned dependency set is tested on Windows, macOS, and Linux by GitHub Actions.

### Windows

First setup:

```text
run_local.bat
```

Subsequent launches:

```text
start.bat
```

### macOS

First setup:

```bash
bash run_local.sh
```

Subsequent launches can use `start.command`.

### Linux

First setup:

```bash
bash run_local.sh
```

Subsequent launches:

```bash
bash start.sh
```

The local application opens at `http://localhost:8000`.

The browser libraries required for 3D visualization and Excel export are vendored in `frontend/vendor/`, so a **local PDB file can be analyzed without CDN access after installation**. Fetching by PDB ID still requires network access to the RCSB Protein Data Bank.

## Examples

Parameter-robustness development panel:

```text
1CRN  1R69  1SHG  2PPN  4AKE  1OMP  1TIM  5CSC
```

Literature comparison:

```text
4F5S  2VUF  2CBA (pH 6.3 and 8.3)  1UBQ  1JNJ  3ECA  1SUV  1TGU
```

Each condition contains the same user-facing result types produced by InterfaceScout: a complete Excel workbook and a B-factor PDB.

## Reproducibility checks

```bash
python -m unittest discover -s tests -v
python validation/verify_publication.py
```

Cross-platform CI repeats these checks under Python 3.10, 3.11, and 3.12 on Windows, macOS, and Linux.

## License

InterfaceScout is released under the MIT License. See `LICENSE`.

## Citation

Citation metadata are provided in `CITATION.cff`. Please cite the associated InterfaceScout publication once bibliographic details are available.
