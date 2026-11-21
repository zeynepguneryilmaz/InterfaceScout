# InterfaceScout

[![Cross-platform CI](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/workflows/publication-ci.yml/badge.svg)](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/workflows/publication-ci.yml)

InterfaceScout is an open-source **physicochemical protein-surface mapping and candidate-region analysis** tool for preliminary experimental or computational study design. It is not an adsorption-site predictor.

The program prepares the selected protein structure once and calculates 11 canonical surface-chemistry maps. Candidate surface patches are ranked using chemistry support, side-chain accessibility, patch coherence, and coarse orientation coherence without material-specific fitted weights. The ranked candidates identify structural regions worth examining, not observed adsorption sites.

## Quick start — three clear workflows

1. **Single protein:** Load a local PDB file or retrieve a PDB ID, adjust pH if relevant and select **Analyze**. Inspect any of the 11 chemistry maps independently or select several to investigate exact residue overlap and co-localization.
2. **Compare two:** Analyze the first structure, select a second PDB (file or ID), then choose **Compare 11 maps**. Both structures are independently analyzed across all 11 maps. Matched-residue comparisons and consensus are provided when the structures permit valid matching.
3. **Multi-PDB:** Upload or queue a set of PDB structures, then select **Analyze selected PDBs**. Explore each structure, matched-residue comparisons and descriptive consensus where applicable.

**One complete download:** **Download complete results · ZIP** always includes all **11 maps for each analyzed structure**, regardless of the individual map or combination selected on screen. For each structure it includes 11 per-map PDFs, 11 genuine browser-captured 3D PNG screenshots embedded in the PDF reports, a complete Excel workbook covering all channels, per-residue/per-patch CSVs, the original PDB, 11 canonical B-factor PDBs, Top-1/3/5 membership B-factor PDBs and detailed JSON. For Compare and Multi-PDB, each protein gets its own complete export directory, with pairwise/consensus outputs when defined. Choosing multiple channels only adds overlap/co-localization content; it never filters out any of the 11 canonical channels.

Complete Results now assembles a ZIP **on the local InterfaceScout Python server**
and starts a native browser attachment download. A visible **Download ZIP now**
link points to that same file and can be clicked again if the browser did not
save it automatically. This is a normal HTTP download, not a short-lived browser
Blob link. Prepared downloads expire after one hour, and the local server keeps
at most 20 recent exports. All 11 maps are included for every structure.

Advanced solution metadata and individual-file downloads are available from expandable sections. pH affects the current computation; ionic strength and temperature are recorded as metadata only.

## Research-first workflows (new, opt-in interface)

The **Quick views** panel applies a compact set of chemistry views immediately with one click:
independent viewing of all eleven maps, hydrophobic + aromatic, positive + H-bond donor,
and negative + H-bond acceptor neighborhoods.
Each choice explains the chemistry it juxtaposes and is **not an algorithmic recommendation
or an adsorption prediction**. Every one of the eleven canonical maps remains available.
The H-bond donor and acceptor channels are intentionally **not** joined by an assumed
exact same-residue intersection: donor/acceptor complementarity and spatial co-location
must be interpreted separately.

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

## Branding and desktop launchers

The approved InterfaceScout icon is used consistently as `frontend/logo.png` (in-app header), `frontend/favicon.png` (browser tab), and `assets/InterfaceScout.ico` (Windows shortcut). The macOS app builder uses the same logo to make an `.icns` icon, while the Linux desktop entry uses the PNG. The surrounding square is part of the app icon; everything outside its rounded corners is transparent. The original atom-style mark is replaced, not packaged as a second competing logo.


## Research PDF: Surface Atlas

The Complete Results ZIP retains **all 11 individual map PDFs** with a visual first-page summary, a categorical rank-color guide that matches the actual browser viewer, the first five ranked regions, and complete ranked-patch / residue appendices. Every analyzed protein also has both `full_report.pdf` (the legacy filename) and `Complete_Surface_Atlas.pdf` (the same complete 11-map visual report), keeping compatibility with existing scripts.

**Scientific interpretation:** Browser map colors indicate ranked patch membership, not a continuous propensity gradient. The PDFs must not be read as continuous numeric heatmaps. Numeric values (scRSA, local scores, propensity, persistence, accessibility) remain in the actual tables and Excel/CSV exports. The PDF includes true browser-rendered 3D screenshots, not server-generated molecular poses. Consensus Top-5 membership fractions are derived from the supplied structures, not measured binding probabilities. No user data are filtered from ZIPs when changing quick views.

Cross-platform browser CI exports a real-data example (`1CRN_hydrophobic_map.pdf`, `1CRN_Complete_Surface_Atlas.pdf`), renders the first pages, and attaches the proofs to the workflow run for visual review.

## Surface Atlas PDF navigation, identity, and publication graphics

- **Protein identity** uses the original PDB `HEADER` accession and `COMPND MOLECULE` annotations when present. Identifiers provided in a PDB-ID request are explicitly labeled as *requested*, not independently verified against a PDB header. No protein name is inferred from an uploaded filename. Every complete atlas records the selected chain, model version, UTC export time and input identity provenance. Source coordinates and analysis parameters remain unchanged.
- **Navigable reports:** the Complete Surface Atlas PDF contains a clickable contents page, native PDF sidebar bookmarks for all eleven chemistry channels and full-data sections, plus an appendix for all original model/solution settings. The consensus and comparison PDFs also provide native clickable navigation. Atlas map tables emphasize **chemistry support**, accessibility, patch size and ranked center residues; no score weights are changed.
- **Evidence-grounded P1–P3 annotations:** the real browser-rendered map images place small labels outside their patch-center atoms and join them with a short leader line to those atoms. The line is visual annotation, not a simulated contact.
- **Opposite-side image when useful:** only if two of the top three patch centers lie on opposite hemispheres of the submitted protein (center-radial cosine below −0.35), an additional actual viewer capture is made after 180° rotation. It is saved as `*_reverse.png` and included on a separate PDF page. An optional secondary capture never replaces or prevents the required eleven primary maps. Reverse views are of **the same protein**, not a second conformation or a superposition.
- **Publication PNG:** the expandable **Individual files and advanced downloads** section contains **Render publication PNG**. With an individual chemistry map open, this asks the current 3Dmol viewer to re-render the same structure at 1920×1080, verifies actual output pixel dimensions and offers a genuine PNG download. This is **not interpolation or an AI drawing**. The regular complete ZIP continues to contain *all eleven* canonical map PDFs, screenshots, PDBs and tabular files regardless of a selected Quick View or publication snapshot. The scientific values in XLSX/CSV/JSON are authoritative.
- **Interpretation:** individual map screenshot colors indicate *categorical patch rank*, not continuous score magnitude. Consensus gradients indicate the fraction of supplied structures with Top-5 membership, not adsorption probability. Comparison Top-5 Jaccard quantifies matched-residue set overlap, not model accuracy. There is no automatic conformational alignment inferred solely from two separate 3D snapshots.
