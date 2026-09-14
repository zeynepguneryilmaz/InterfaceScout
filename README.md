# InterfaceScout

InterfaceScout is an open-source tool for identifying and prioritizing plausible protein-side contact regions at material interfaces from a protein structure and solution conditions.

The program prepares the selected protein structure once and calculates 11 canonical surface-chemistry maps. Candidate surface patches are prioritized using chemistry support, side-chain accessibility, patch coherence, and coarse orientation coherence without material-specific fitted weights.

## Surface-chemistry maps

InterfaceScout calculates:

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

## Inputs

- PDB ID or local PDB file
- optional chain selection
- pH
- ionic strength
- temperature

In the current model, pH enters the scoring calculation. Ionic strength and temperature are retained as run metadata.

## Outputs

The browser application provides:

1. **Complete Excel workbook** containing run/method information, map definitions, residue-level values, ranked patches, and patch members for all 11 chemistry maps.
2. **B-factor PDB** for the currently displayed chemistry map, with the 0–100 InterfaceScout residue propensity written to the standard B-factor field for visualization in molecular-graphics software.

The Excel workbook contains `Run_Method`, `Map_Definitions`, `Residue_Values`, `Patches`, and `Patch_Members`. `Composite_Selected` is added only when a multi-map Composite overlay is selected.

## Installation and use

InterfaceScout supports Python 3.10–3.12.

### Windows

Double-click or run:

```text
run_local.bat
```

### macOS

First setup:

```bash
bash run_local.sh
```

After setup, `start.command` can be used to launch InterfaceScout again.

### Linux

First setup:

```bash
bash run_local.sh
```

After setup, `start.sh` can be used to launch InterfaceScout again.

The local application opens at `http://localhost:8000`.

## Examples

The `examples/` directory contains only the structures reported in the manuscript.

Structures used for parameter robustness:

```text
1CRN  1R69  1SHG  2PPN  4AKE  1OMP  1TIM  5CSC
```

Structures used for literature comparison:

```text
4F5S  2VUF  2CBA (pH 6.3 and 8.3)  1UBQ  1JNJ  3ECA  1SUV  1TGU
```

Each condition contains only the same user-facing files produced by InterfaceScout: the complete Excel result workbook and the B-factor PDB for the displayed chemistry map.

## Fixed settings used in the manuscript

- Shrake–Rupley probe radius: **1.40 Å**
- Shrake–Rupley sampling: **200 points/atom**
- surface threshold: **scRSA ≥ 0.05**
- multiscale aggregation radii: **6 and 9 Å**
- coarse patch radius: **8 Å**

## License

InterfaceScout is released under the MIT License. See `LICENSE`.

## Citation

Please cite the associated InterfaceScout publication once bibliographic details are available.
