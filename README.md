# InterfaceScout

[![Cross-platform CI](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/workflows/publication-ci.yml/badge.svg)](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/workflows/publication-ci.yml)

InterfaceScout is an open-source **physicochemical protein-surface mapping and candidate-region analysis** tool for preliminary experimental or computational study design. It is not an adsorption-site predictor.

The program prepares the selected protein structure once and calculates 11 canonical surface-chemistry maps. Candidate surface patches are ranked using chemistry support, side-chain accessibility, patch coherence, and coarse orientation coherence without material-specific fitted weights. The ranked candidates identify structural regions worth examining, not observed adsorption sites.

## Research-first workflows (new, opt-in interface)

The **Research question** selector offers five practical starting workflows: exploring an unfamiliar
protein, planning noncovalent protein–material experiments, selecting residues for surface
immobilization/labeling, preparing computational interface simulations, and comparing provided
protein conformations. These choices are **suggestions only**: users can select any of the
11 maps and never need to identify a particular material to begin.

The new decision-support layer keeps two physically distinct data types separate:

1. **Direct protein-side descriptors:** side-chain solvent accessibility, exposed residue
   identities and simple, explicitly defined amino-acid groups such as nonpolar, aromatic,
   cysteine, lysine, acidic and basic residues. These groups are descriptive; they do not
   calculate microscopic chemical reactivity, reaction rates or a full electrostatic potential.
2. **Eleven canonical compatibility maps:** unchanged chemistry rules and patch ranking.
   Users can inspect maps independently or combine two or more. Combination reports
   **exact residue intersection** separately from **spatial co-localization** of different
   compatible residues within the same existing 8 Å patch. Combined maps are not fitted
   adsorption scores.

**Conformational comparison** accepts a second PDB file or PDB ID and compares the
same canonical analyses under the same pH. Outputs include accessibility changes,
per-channel Top-1/3/5 region-membership Jaccard overlaps, and added/removed residues.
The user must provide matched residue numbering/chain identities (the program checks
mapping identity and rejects insufficient overlap); it does not silently match unrelated proteins.

**Out of scope:** GNM, ANM, residue-interaction-network analysis, adsorption pose prediction,
adsorption-induced unfolding, binding free energies and simulation of dynamics. If users
already have an MD trajectory, they can export representative structures as PDB files and
compare them using this workflow. Native XTC/TRR/DCD parsing is not implemented.

The research report can be exported to JSON; the original complete Excel and map-labelled
B-factor PDB exports remain available. The original chemistry engine is frozen and **has
not been recalibrated to maximize experimental adsorption contacts**.

## Multi-PDB, consensus, and complete exports

The original single-PDB Analyze, map visualization, composite overlay, complete Excel and current-map B-factor PDB export remain in place.

**Single structure:** load a PDB by ID or local file and select Analyze. The new **Complete export** button downloads a single ZIP containing all 11 independently colored map screenshots, 11 separate map PDF documents, 11 separate canonical-score B-factor PDBs, additional Top-1/3/5 region-membership B-factor PDBs, every original input structure, detailed per-residue/per-patch CSV tables, an Excel workbook, complete structured JSON and a PDF report that contains **actual 3Dmol viewer screenshots** for all 11 maps. If a selected multi-channel intersection is available, the ZIP also contains intersection tables and a combined-map screenshot.

**Multiple structures:** upload several local PDBs at once or add by PDB ID using the Multi-PDB panel (1–20 structures). Each PDB is analyzed separately using the same 11 chemistry channels, and each retains a complete independent export directory. The application shows pairwise comparison tables, including changes in side-chain accessibility and Top-1/3/5 residue-set overlap.

For two or more structures, consensus is calculated **only for residue positions that have the same chain, residue number, insertion code and residue identity**; insufficiently matched structures are rejected. The consensus shows per-residue accessibility, chemical-channel compatibility and Top-1/3/5 member frequency, as well as strict (all structures) and majority (more than half) shared residue sets. Map intersections can also be described across the structures. Consensus visualizations can be viewed in the original 3D viewer, and the ZIP includes a consensus report PDF, per-map consensus screenshots, detailed CSV/Excel files, and consensus-frequency B-factor PDBs (display scale 0–100).

**Interpretation:** consensus frequencies describe the *supplied structure collection*; they are not equilibrium occupancies, binding probabilities, contact predictions or conformational time courses. The downloaded B-factor column contains InterfaceScout visualization numbers, **not** experimental crystallographic B-factors. The ZIP generator refuses to substitute fake images when viewer screenshots cannot be captured.

**Intentionally not implemented:** raw MD trajectories, GNM, ANM and RINpy. Researchers can export previously simulated selected frames as PDB files and compare them, but no trajectory analysis or elastic-network conformer generation is claimed.

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
- pH (affects the current scoring)
- ionic strength (saved as metadata only)
- temperature (saved as metadata only)
- second PDB file or PDB ID for optional same-protein comparison

In the current model, **pH enters the scoring calculation**. Ionic strength and temperature are retained as run metadata and should not be interpreted as explicit electrostatic-screening or thermal-response calculations.

## Intended scope and limitations

InterfaceScout is intended for **reproducible characterization of accessible physicochemical protein-surface regions**, and for comparing candidate regions across user-provided structures. Material-compatibility maps support hypothesis generation rather than adsorption claims.

The existing chemistry-based publication benchmark is limited and does not establish general adsorption predictive performance. Structural characterization is most straightforward for folded globular proteins and folded oligomeric proteins with representative atomic structures. Use more cautiously for intrinsically disordered proteins, strongly unfolding/restructuring adsorption regimes, membrane-embedded systems, covalent attachment, or systems whose relevant chemistry depends critically on cofactors, glycans, lipids, non-standard residues, or other hetero components removed by the standard preparation policy.

InterfaceScout does **not** predict absolute adsorption free energy, adsorbed amount, atomistically unique binding orientation, adsorption kinetics, or adsorption-induced unfolding. Ranked patches are interpretable protein-side candidate regions for experimental follow-up or higher-resolution simulation, not exact atomistic binding sites or predicted final contacts.

The literature benchmark contains seven primary conditions and two secondary challenge cases. The transferrin-polystyrene and catalase-polystyrene cases are retained specifically to expose the limitation of a static folded-state model when adsorption-associated structural change is important. See `validation/README.md`.

## Outputs

1. **Complete Excel workbook** with run/method information, map definitions, residue-level values, ranked patches, and patch members for all 11 maps.
2. **B-factor PDB** for the displayed map, with the 0-100 InterfaceScout residue propensity written to the standard B-factor field.
3. **Research decision report (JSON)** including direct protein properties, chosen maps, optional co-localization analysis and optional PDB comparison.

## Installation and use

This is **one cross-platform application**, not three separate codebases. Use Python **3.10–3.12**.

| Platform | First setup | Subsequent launch | Instructions |
|---|---|---|---|
| Windows | `run_local.bat` | `start.bat` | [Windows guide](docs/INSTALL_WINDOWS.md) |
| macOS | `bash run_local.sh` | `bash start.command` | [macOS guide](docs/INSTALL_MACOS.md) |
| Linux | `bash run_local.sh` | `bash start.sh` | [Linux guide](docs/INSTALL_LINUX.md) |

The local interface is at **http://localhost:8000**. After dependencies are installed, subsequent launches do not run pip. For pip SSL certificate-chain errors, follow the relevant platform guide rather than disabling certificate verification. Local PDB files work without a network connection after installation; PDB ID retrieval still requires RCSB access.

When a new `v*` tag is published, the [release workflow](.github/workflows/platform-release.yml) builds Windows, macOS and Linux **source archives** containing the corresponding launchers. These are not standalone `.exe`/`.app` binaries.

## Examples

See the [examples index](examples/README.md) for the existing **16-protein, 17-condition** archive of real Excel and annotated-PDB exports, the [1LYZ reproducible example](examples/1LYZ/README.md), and the frozen benchmark source in [validation](validation/README.md). The example exporter needs a running InterfaceScout app and does not invent absent results.

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
