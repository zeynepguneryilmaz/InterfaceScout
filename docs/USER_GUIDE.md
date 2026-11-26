# InterfaceScout | User Guide

**Physicochemical protein-surface mapping · Candidate-region exploration · Structure comparison**  
**Software:** InterfaceScout 1.1.0 · **Computational model:** `2.5.0-dual-accessibility`  
**Guide edition:** 1.0 · **Updated:** 23 November 2026 · **Language:** English

> **At a glance.** InterfaceScout maps accessible, protein-side physicochemical features, identifies descriptive candidate surface patches, and compares compatible structures. It is a research decision-support tool. It **does not** predict a unique adsorption site or pose, binding free energy, adsorption kinetics, or adsorption-induced conformational changes.

**Source code:** [InterfaceScout on GitHub](https://github.com/zeynepguneryilmaz/InterfaceScout) · **Citation:** [CITATION.cff](../CITATION.cff) · **Examples:** [Examples index](../examples/README.md)

## Contents

1. [Quick start](#1-quick-start)
2. [Installation and requirements](#2-installation-and-requirements)
3. [Choose the right protein structure](#3-choose-the-right-protein-structure)
4. [Workspace and first analysis](#4-workspace-and-first-analysis)
5. [Understanding the 11 surface maps](#5-understanding-the-11-surface-maps)
6. [Combining maps and Regional Chemistry Profile](#6-combining-maps-and-regional-chemistry-profile)
7. [Comparing two conformations](#7-comparing-two-conformations)
8. [Multi-PDB analysis](#8-multi-pdb-analysis)
9. [Exporting and reporting results](#9-exporting-and-reporting-results)
10. [Worked examples](#10-worked-examples)
11. [Experimental comparison: evidence-aware use](#11-experimental-comparison-evidence-aware-use)
12. [Model details and limitations](#12-model-details-and-limitations)
13. [Troubleshooting](#13-troubleshooting)
14. [Reproducibility and support](#14-reproducibility-and-support)

---

## 1. Quick start

**Your first map in six steps**

1. **Launch locally.** Open InterfaceScout using the launcher for your operating system. If no browser opens, visit `http://localhost:8000`.
2. **Select `Single protein`.** Enter a four-character **PDB ID** and click **Fetch**, or load a local `.pdb` file. A small, well-resolved example is **1AKI** (lysozyme).
3. **Check the structure.** If appropriate, enter one or more chain identifiers, such as `A` or `A,B`. Confirm the selected representation is the biologically relevant construct or assembly.
4. **Set solution pH and click `Analyze`.** The default pH is **7.4**. All 11 surface-chemistry maps are computed in the same run.
5. **Explore.** Use **Cartoon**, **Surface**, or **Cartoon + surface**; inspect individual maps, patch cards, and **Quick Views**.
6. **Save the complete result.** Select **Download complete results · ZIP**. This is the recommended research archive; it includes numerical results and the 11 map visualizations.

**Before interpreting a highlighted region:** inspect its residue IDs, the accessible surface coverage, the map definition, and any missing structural components. Highlighted patches are **candidates for investigation**, not confirmed adsorption contacts.

## 2. Installation and requirements

InterfaceScout runs on your computer as a **Python-backed local application** with a browser-based 3D interface. It is not a precompiled stand-alone executable.

| Requirement | Supported configuration |
| --- | --- |
| Operating system | Windows, macOS, Linux |
| Python | 3.10, 3.11, or 3.12 |
| Browser | Modern browser with WebGL support |
| Connectivity | Required for initial dependency installation (unless a compatible offline wheelhouse is provided) and for fetching PDB IDs; local PDB files can be used after installation without retrieving structures online |
| Local address | `http://localhost:8000` |

**Installation procedure.** Download the repository using **Code → Download ZIP**, extract the ZIP fully to a permanent directory, and start InterfaceScout from the extracted folder.

| System | First launch | Subsequent launch |
| --- | --- | --- |
| Windows | Double-click `run_local.bat` | `start.bat` |
| macOS | `bash run_local.sh` | `bash start.command` |
| Linux | `bash run_local.sh` | `bash start.sh` |

If a launcher does not open the browser automatically, navigate to `http://localhost:8000`. Windows users can create an optional desktop shortcut using `Create Desktop Shortcut.bat`; this is a launcher and still requires the extracted software folder. For platform-specific details see [Windows](INSTALL_WINDOWS.md), [macOS](INSTALL_MACOS.md), and [Linux](INSTALL_LINUX.md).

**Installation security.** If `pip` reports an HTTPS certificate error, use your institution's verified certificate bundle or a trusted network. Do **not** disable TLS/certificate verification. The standard source ZIP should not be assumed to include offline dependency wheels.

## 3. Choose the right protein structure

Structure selection is the most important scientific quality-control step. A high-quality model of the **wrong construct** is less useful than a suitable structure with known, manageable limitations.

| Check | Practical question |
| --- | --- |
| Identity | Is this the intended protein, species, isoform, mutant, and experimental construct? |
| Chains / assembly | Is the relevant surface a monomer, dimer, heteromer, or higher-order biological unit? |
| Completeness | Are the regions of interest, side chains, termini, or domains absent? Is the file predominantly Cα-only? |
| Structural state | Is the structure apo/holo, open/closed, active-like/inactive-like, or ligand-bound? |
| Chemistry | Are glycans, cofactors, metals, heme, modified residues, and bound ligands relevant to your question? |
| Correspondence | If comparing with an experiment, can the literature sequence and residue numbering be mapped reliably onto this PDB? |

**Practical rules:**

- The uploaded structure is the **object being mapped**; it need not reproduce every experimental condition for exploratory work.
- The more specific the claim (especially a residue-level experimental comparison), the stronger the required protein/construct/sequence match.
- Choose and record the **biological assembly** explicitly when oligomerization changes exposed interfaces. A default single chain may not be the functional unit.
- Several unsupported components are removed or not chemically modeled during protein-only preparation; do not interpret a metal- or phosphorylation-driven site as though its original chemistry were retained.
- NMR models and snapshots from simulations can be analyzed individually, but the tool does not infer their statistical populations or time evolution.

## 4. Workspace and first analysis

The workspace has three workflows:

| Workflow | Best use | What it provides |
| --- | --- | --- |
| **Single protein** | Examine a structure before experiments or docking/MD | 11 maps, candidate patches, residues, chemistry combinations |
| **Compare two** | Compare two compatible structures of the **same protein** | Matched-residue accessibility changes and Top-1/3/5 map comparisons |
| **Multi-PDB** | Explore 2–20 supplied structures | Per-structure 11-map analyses and descriptive comparisons/consensus where residue correspondence permits |

### Input controls

- **PDB ID → Fetch:** retrieve a structure from RCSB (requires network access).
- **Load a local PDB file:** upload a `.pdb` with protein `ATOM` records.
- **Chain (optional):** single chain `A`, or several chains such as `A,B`; inspect the structure before excluding a chain.
- **Solution pH:** used by the current chemistry scoring.
- **Advanced conditions:** ionic strength and temperature are recorded for the canonical chemistry maps and used as physical inputs when optional continuum electrostatics is enabled.
- **Electrostatics:** enable the optional PDB2PQR/PROPKA and APBS calculation when the external executables are available on PATH. NPBE is the default; LPBE and paired mode (`both`) are alternatives. The primary residue-associated surface potential is an area-weighted approximate solvent-accessible-shell value; the side-chain atom-coordinate potential is a separate descriptor. Canonical map values and ranks remain unchanged. See the [method overview](../README.md#optional-continuum-electrostatics).
- **Analyze:** performs one full 11-channel calculation.

### Inspecting the 3D viewer

Drag to **rotate**, use the wheel to **zoom**, and hover to inspect residues. The **Home / all** control returns to the whole protein. Select **Cartoon**, **Surface**, or **Cartoon + surface** depending on whether you want to see fold context or the accessible external geometry. Use the map tabs and candidate patch cards to investigate ranked regions; examine numerical residue/patch tables before drawing scientific conclusions.

**Visual convention:** patch colors and ranks are *categorical priorities within a map*. They are **not** a calibrated affinity scale or experimentally measured heatmap.

## 5. Understanding the 11 surface maps

All 11 maps are computed on every run. A channel describes a **protein-side surface feature or material-compatibility rule**; it is not an adsorption probability for an actual nanoparticle, coating, or crystal face.

| Channel | What to inspect |
| --- | --- |
| `anionic` | Regions compatible with negatively charged interface chemistry |
| `cationic` | Regions compatible with positively charged interface chemistry |
| `hydrophobic` | Exposed nonpolar features |
| `pi_carbon` | Aromatic/π-associated features |
| `hbond_donor` | Compatibility with hydrogen-bond donating surface groups |
| `hbond_acceptor` | Compatibility with hydrogen-bond accepting surface groups |
| `oxide` | Chemistry rules associated with metal-oxide-like interfaces |
| `calcium_phosphate` | Features relevant to calcium/phosphate-rich charged sites |
| `metal_coord` | Environments involving potentially coordinating exposed side chains |
| `soft_metal_sulfur` | Sulfur-associated compatibility in soft-metal contexts |
| `phosphate` | Chemistry rules associated with phosphate-rich surfaces |

**Quick Views** provide one-click *descriptive* combinations: **All 11 maps**, **Hydrophobic + aromatic**, **Positive + donor**, and **Negative + acceptor**. Quick Views do not learn material preferences or supply predicted binding poses.

**Nonindependence matters.** The `cationic` and `hbond_donor` channels can share residue chemistry and therefore produce identical Top-5 residue sets under a given condition. Therefore 11 labeled maps must **not** be described as 11 independent sources of biological evidence. Check overlap and complementarity before using a multi-channel finding as corroboration.

## 6. Combining maps and Regional Chemistry Profile

Select **two or more maps** to explore chemical features in the same region. The combined view distinguishes two different ideas:

| Term | Meaning |
| --- | --- |
| **Exact residue overlap** | The **same residue identifier** belongs to the selected map lists. |
| **Spatial co-localization** | Features occur inside the **same center-defined 8 Å candidate patch**. Different supporting residues can be farther than 8 Å from each other. |

The **Regional Chemistry Profile** supplies three interpretable descriptors without introducing a new combined affinity score:

1. **Per-channel mean local support (0–1):** average original local support across patch residues, with absent support contributing zero.
2. **Supporting residue count/fraction:** number and proportion of patch residues with positive original channel support.
3. **Minimum distinct-residue Cα separation (Å):** the shortest geometrical Cα distance between **different** residues supporting different selected maps; multi-map exports also report pairwise distances.

`NA` for separation means no eligible distinct supporting pair or missing coordinates—not 0 Å. Cα distance is a coarse geometric descriptor, **not** an atomic contact distance, hydrogen bond, adsorption orientation, or energy.

**Recommended reading order:** inspect each map separately → evaluate exact overlap → evaluate co-localization → examine support fractions and distances → export the numerical residue list. Correlated channels do not count as independent support.

## 7. Comparing two conformations

Use **Compare two** only for structures of the **same protein** with reliably matched residue identifiers.

1. In **Single protein**, load and **Analyze** the first structure.
2. Select **Compare two** and enter the **Second PDB ID** or select a **Second local PDB**; specify the correct second chain if needed.
3. Click **Compare 11 maps**. Both structures are analyzed using the **same pH**.
4. Read the matched-residue count, changes in solvent-exposure classification, and per-channel **Top-1, Top-3, and Top-5 Jaccard overlaps**.
5. Export the complete result to retain **both** structures and their comparisons.

**Top-k Jaccard** is `|A ∩ B| / |A ∪ B|` for the relevant selected residue sets. A value of 1 indicates identical selections; 0 indicates no selected members in common. This is **not** an equilibrium probability or a prediction of changed adsorption.

| Example pair | Structural question | Interpretation caveat |
| --- | --- | --- |
| `1F6R` / `1F6S` | α-Lactalbumin, apo vs Ca-bound | Calcium is excluded from the protein-only chemistry model. |
| `2HAV` / `3QYT` | Transferrin, apo vs iron-bound | Iron is excluded; mapped differences describe the represented protein structures. |
| `4CFE` / `5ISO` | AMPK heterotrimer, different phosphorylation states | Nonstandard phosphorylated residues may be excluded; this does not isolate phosphorylation chemistry. |
| `2GS2` / `2GS7` | EGFR kinase, active-like vs inactive-like | Bound ligands, constructs, and mutation differences can confound the comparison. |

**Compatibility rule:** comparison uses chain ID, residue number, insertion code, and residue identity for matching. It does **not** perform hidden sequence alignment or structural superposition, and it does not ingest native MD trajectories.

## 8. Multi-PDB analysis

**Multi-PDB** accepts **2–20 structures per queued analysis**. Add local PDB files, add the current structure, or add another four-character PDB ID; verify the list and select **Analyze selected PDBs**.

Each structure retains its own **11 complete maps**. Descriptive residue consensus is interpretable only if matching positions represent the same biological residue. If a collection mixes different proteins or nonhomologous sequences, compare their maps **separately** rather than inventing a common residue-level consensus. A percentage across supplied frames describes the supplied set, **not** a Boltzmann population, MD time occupancy, or adsorption probability.

## 9. Exporting and reporting results

The recommended archive is **Download complete results · ZIP**. It is more complete and more reproducible than saving a single screen capture.

| Export | Contents and intended use |
| --- | --- |
| **Complete Surface Atlas / map PDFs** | All 11 map views and report/method context |
| **PNG screenshots** | Images captured from the real browser-based 3D viewer |
| **Complete Excel (`.xlsx`)** | Residue/patch tables, map values, settings, and selected Regional Chemistry Profile outputs |
| **CSV** | Machine-readable numerical tables for independent analysis |
| **JSON** | Full machine-readable results and provenance |
| **B-factor PDB files** | Map display values written into the crystallographic B-factor field for visualization software |
| **Source PDB and provenance** | Structure identity, input and reproducibility context |
| **Publication PNG** | Optional real 1920 × 1080 rendering of the selected map view via **Individual files and advanced downloads** |

**Important B-factor warning:** values written to a B-factor PDB by InterfaceScout are **software visualization scores**, not experimental atomic displacement parameters. Do not interpret them crystallographically.

### Minimum methods/reporting record

Record: InterfaceScout public and model version; input PDB and archive date; source/biological assembly/chains/model; known missing fragments and excluded components; pH; selected channel(s); region rank and residues; selection coverage; exact overlap versus spatial co-localization; and, for experimental comparisons, study citation, evidence resolution, and matched residue count. Preserve the original ZIP without manual changes to its numerical results.

## 10. Worked examples

### Example A · Start with lysozyme (`1AKI`)

**Question:** Which external regions of a small globular protein present hydrophobic or aromatic-compatible features?

1. Load `1AKI`, confirm the selected chain, and run **Analyze** at pH 7.4.
2. Inspect **Hydrophobic** and **π / aromatic** independently.
3. Open **Hydrophobic + aromatic** Quick View. Compare exact shared residues against nearby but distinct supporting residues.
4. Examine the Regional Chemistry Profile for the highlighted patches.
5. Export the complete ZIP and record candidate residue IDs. These are potential **follow-up investigation targets**, not lysozyme adsorption-site predictions.

### Example B · Vitronectin domain / hydroxyapatite (`6O5E`)

**Question:** Can a chemistry-informed candidate region be meaningfully compared with independently described Asp-associated experimental observations?

1. Load the `6O5E` **vitronectin domain** (not the whole protein).
2. Inspect `calcium_phosphate` and `oxide` separately, then inspect their regional overlap.
3. Cross-map experimental residues to this structure and calculate Top-k recovery **at equal accessible surface coverage**.
4. Compare with a simple exposed acidic-residue baseline before making any claim of superiority.

In the existing retrospective pilot audit, three of four accessible experimental Asp labels appeared among the selected Top-5 `calcium_phosphate` patch members; a simple charge-based control recovered four of four. This illustrates **plausible chemistry, not superior adsorption-site prediction**. The experimental label and structure-domain limitations must be disclosed.

## 11. Experimental comparison: evidence-aware use

A protein can be useful in InterfaceScout even when there is **no direct experimental adsorption map**. However, what can be claimed depends on the **resolution and independence of evidence**.

| Evidence tier | Source measurement | Supported use |
| --- | --- | --- |
| **A** | Experimentally localized residues/peptides with relevant structural matching | Quantitative region/residue overlap, with measurement ambiguity retained |
| **B** | Experimental domain or broad region | Regional consistency, not precise residue-level accuracy |
| **C** | Adsorption quantity, kinetics, corona membership, or an interaction signal without a localized region | Contextual material–protein case study; **not** binding-site validation |
| **D** | No matched protein–material adsorption evidence | Exploratory candidate-region generation and testable hypotheses |

For a fair comparison: select channels by chemistry **before viewing experimental outcomes**; verify species, construct, assay state, numbering and visibility; evaluate **Top-1/3/5 plus selection coverage**; compare with equal-coverage SASA, charge/amino-acid-class, and randomized-patch baselines; and cluster repeated outcomes by underlying protein/publication to avoid pseudoreplication. Report both supported and unsupported examples.

The publication study uses **TV48** (48 structures for technical verification) and **ISV76** (76 structures for application and surface-chemistry analysis). Literature-associated examples illustrate interpretation of the maps in real research contexts. Download the [research datasets](../datasets/README.md) and consult their [provenance index](../validation/PUBLICATION_EVIDENCE.md).

## 12. Model details and limitations

| Component | Current behavior |
| --- | --- |
| Public software version | `1.1.0` |
| Internal model version | `2.5.0-dual-accessibility` |
| pH | Used in chemical scoring; default 7.4 |
| Ionic strength, temperature | Canonical-map metadata; optional APBS inputs |
| SASA method | Shrake–Rupley; 1.40 Å probe, 200 sample points per atom |
| Surface exposure | Total-residue relative accessibility (`tRSA ≥ 0.05`), normalized using Tien et al. (2013) maxima |
| Chemical side-chain support | Non-backbone side-chain SASA / reference ASA (`scRSA`); glycine has no side-chain heavy atoms (`scRSA = 0`) |
| Patch mean accessibility | Mean geometric `tRSA`, distinct from chemical `scRSA` |
| Multiscale radii | 6 and 9 Å |
| Candidate patch neighborhood | 8 Å center-defined scale |

**Outside the model's demonstrated scope:** adsorption affinity/ΔG; unique protein pose or material contact geometry; nanoparticle curvature, real surface grafting density, and detailed solvent/salt effects; adsorption-induced unfolding; kinetics; bound cofactor/glycan/PTM chemistry; native MD trajectory parsing; and conformational population prediction. Missing residues, unresolved flexible segments and protein-only preparation can materially change a result.

## 13. Troubleshooting

| Symptom | Check / action |
| --- | --- |
| **Analyze** is disabled | Ensure a `.pdb` with protein `ATOM` records was successfully loaded. |
| PDB ID **Fetch** fails | Verify the four-character ID, internet/RCSB access, or load a local PDB. |
| 3D viewer is blank | Verify WebGL support, refresh, and test with a small protein; inspect the structure file and selected chains. |
| Unexpected residue omissions | Examine missing regions, noncanonical residues, side chains, chain filtering and structure preparation. |
| **Compare 11 maps** reports no matches | Check same-protein identity, chain/residue numbering, insertion codes and residue identity. |
| Large ZIP does not complete | Confirm sufficient browser/memory/disk resources and download permissions; test with one small structure. |
| SSL certificate installation error | Use a trusted network, verified CA, or approved offline wheelhouse. Never disable TLS verification. |
| Two maps look identical | This can be a real channel correlation; compare numerical sets and report nonindependence. |

For Windows startup failures inspect `backend/startup_log.txt` and validate the installed packages using `backend/.venv/Scripts/python.exe backend/dependency_check.py`.

## 14. Reproducibility and support

From the repository root, after installing the backend requirements, run:

```bash
python -m pip install -r backend/requirements.txt
python -m unittest discover -s tests -v
python validation/verify_publication.py
```

Run these commands from a full repository checkout. Technical checks and publication validation are described in the [validation folder](../validation/README.md). When citing InterfaceScout, use the citation metadata in [CITATION.cff](../CITATION.cff); cite the underlying experimental source separately when discussing literature comparisons. For reproducible problems open a [GitHub Issue](https://github.com/zeynepguneryilmaz/InterfaceScout/issues) with the OS, Python version, minimal shareable structure, and steps; **never** publish confidential/unreleased structural data in a public issue.

**Documentation note.** This guide describes InterfaceScout as implemented in the main-branch source and public README for the v1.1.0 publication model. The worked-example observations are illustrative retrospective research findings, **not** a training set or a built-in automated experimental validation feature. If code or model behavior changes, review this guide alongside the software release and model version.

---

© 2026 Ö. Zeynep Güner Yılmaz · InterfaceScout · Original software/documentation distributed under the repository [MIT License](../LICENSE). Third-party dependencies retain their respective licenses; see [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md).