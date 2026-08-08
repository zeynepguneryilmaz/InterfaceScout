# InterfaceScout

**Residue-level protein surface chemistry and surface-chemistry compatibility mapping.**

InterfaceScout is an open-source, deterministic computational framework for analyzing the solvent-exposed surface of protein structures and identifying residue-level compatibility with generic classes of material surface chemistry.

From a protein structure and user-defined solution conditions, InterfaceScout combines solvent-accessible surface area, local three-dimensional residue organization, pH-dependent protonation, and Poisson–Boltzmann electrostatic potential with literature-informed residue–chemistry reference strengths.

The framework generates residue-level interaction-propensity maps and chemistry-specific patch-density maps for eleven interfacial chemistry classes:

- cationic
- anionic
- hydrogen-bond donor
- hydrogen-bond acceptor
- π / carbon-like
- hydrophobic
- oxide
- hydroxyapatite / Ca²⁺
- transition-metal coordination
- gold
- phosphate

InterfaceScout provides a protein-centered description of chemically compatible surface motifs and spatially enriched residue patches that can support the interpretation and comparison of protein–material interactions across chemically distinct interfaces.

**No machine learning or fitting to experimental adsorption data is used.**

The application runs locally on your computer and opens in your web browser at:

`http://localhost:8000`

Protein structures and analysis results are processed locally.

---

## What's in this folder

    InterfaceScout/
    ├── backend/          The analysis engine (Python)
    ├── frontend/         The web interface
    ├── run_local.bat     First-time setup  (Windows)
    ├── start.bat         Daily launcher    (Windows)
    ├── run_local.sh      First-time setup  (macOS / Linux)
    ├── start.command     Daily launcher    (macOS)
    ├── start.sh          Daily launcher    (Linux)
    ├── interfacescout.ico / .png   App icon
    └── README.md         This file

**Keep all of these together in one folder.** The launchers sit next to `backend/` and find it automatically — don't move them into subfolders.

---

## Requirements

- **Python 3.11 or 3.12** with SSL support
- An internet connection during initial setup to download software dependencies and, where required, the APBS electrostatics binary
- Supported operating systems:
  - **Windows 10/11**
  - **macOS** — Intel or Apple Silicon
  - **Linux**

---

## Windows

1. **Double-click `run_local.bat`** the first time. It installs the required components, places an **InterfaceScout** icon on your Desktop, and starts the application.

2. On later runs, double-click the **InterfaceScout** Desktop icon. The existing installation is detected automatically and the application opens in your browser at:

   `http://localhost:8000`

> The Desktop icon runs `run_local.bat`, which can be used again after the first installation. To force a clean reinstall, delete the `backend\.venv` folder and run the setup again.
>
> If Windows SmartScreen displays a warning for the `.bat` file, choose **More info → Run anyway**.
>
> If the application fails to start, the command window remains open and a log is written to `backend\startup_log.txt`.

---

## macOS

1. Run **`run_local.sh`** once by right-clicking it and choosing **Open**, or run:

   `bash run_local.sh`

   in Terminal.

   The setup installs the required components and places an **InterfaceScout.command** launcher on your Desktop.

2. On later runs, double-click **`InterfaceScout.command`**. The backend starts locally and the application opens in your browser.

> Gatekeeper may request confirmation on the first run. Right-click → **Open** and confirm if needed.

---

## Linux

1. Run:

   `bash run_local.sh`

   once.

   The setup installs the required components and places an **InterfaceScout.desktop** launcher on your Desktop.

2. On later runs, double-click the InterfaceScout Desktop launcher.

> If the launcher is marked as untrusted, right-click and choose **Allow Launching** where supported.

---

## Using InterfaceScout

1. Enter a PDB ID, such as `4F5U`, and click **Fetch**, or upload a PDB file.

2. Define the analysis conditions:
   - **pH**
   - **ionic strength** (mM)
   - **temperature** (K)
   - **patch radius** (Å; default 12)

3. Click **Analyze Surface**.

4. Explore:
   - residue-level surface-feature mapping
   - generic surface-chemistry compatibility mapping
   - interaction-propensity maps
   - chemistry-specific patch-density maps

5. Inspect the three-dimensional protein representation and residue-level analysis tables.

6. Export results as:
   - **CSV**
   - **PDB**
   - **PDF**

---

## Analysis workflow

InterfaceScout analyzes the protein structure independently of an explicit material model.

The workflow includes:

1. Protein structure preparation
2. Solvent-accessible surface analysis
3. Condition-dependent protonation and electrostatics
4. Assignment of residue–surface chemistry interaction classes
5. Residue-level compatibility scoring
6. Chemistry-specific patch-density analysis
7. Three-dimensional visualization and data export

Solvent-accessible surface area is calculated using the **Shrake–Rupley algorithm**.

Condition-specific protonation states and atomic charges are assigned using **PDB2PQR / PROPKA**.

Electrostatic potentials are calculated using **APBS** and the Poisson–Boltzmann equation.

Residue-level compatibility is evaluated using literature-informed residue–chemistry reference strengths together with solvent exposure, local three-dimensional chemical context, electrostatic complementarity, and mechanism-specific protonation weighting.

---

## Surface feature mapping

InterfaceScout identifies chemically relevant features presented by solvent-accessible protein residues, including:

- positive charge
- negative charge
- hydrogen-bond donor
- hydrogen-bond acceptor
- hydrophobic character
- aromatic character
- metal-binding groups
- thiol groups
- carboxyl groups
- amine-containing groups

These features can be visualized directly on the three-dimensional protein structure.

---

## Surface-chemistry compatibility mapping

InterfaceScout evaluates exposed residues against eleven generic classes of surface chemistry:

| Surface chemistry | Main interaction motifs |
|---|---|
| Cationic | electrostatic attraction and hydrogen bonding |
| Anionic | electrostatic attraction and hydrogen bonding |
| H-bond donor | hydrogen-bond complementarity |
| H-bond acceptor | hydrogen-bond complementarity |
| π / carbon-like | π–π and cation–π interactions |
| Hydrophobic | hydrophobic and CH–π contacts |
| Oxide | carboxylate–oxide interactions and hydrogen bonding |
| Hydroxyapatite / Ca²⁺ | Ca²⁺ coordination and related interactions |
| Transition-metal coordination | His/Cys/Asp/Glu/Met coordination |
| Gold | sulfur–gold affinity |
| Phosphate | electrostatic and hydrogen-bond interactions |

For each chemistry class, InterfaceScout reports residue-level **interaction propensity** and spatially resolved **patch density**.

---

## Interaction propensity

Residue-level interaction propensity combines:

- literature-informed residue–chemistry reference strength
- solvent exposure
- local three-dimensional residue context
- electrostatic complementarity where applicable
- mechanism-specific pH-dependent protonation

The resulting normalized propensity values allow compatible residues to be ranked within each surface-chemistry map and analysis condition.

---

## Patch-density analysis

Patch-density analysis identifies spatial regions enriched in exposed residues compatible with the same surface chemistry.

For each surface position, InterfaceScout integrates compatible residue contributions within a user-defined three-dimensional neighborhood.

The default patch radius is **12 Å**.

This provides a complementary view of:

- strong individual residue contributors
- spatially coherent chemistry-enriched surface regions

---

## Theory and references

The **Theory** section within the application documents:

- the scoring framework
- mathematical expressions used in the analysis
- solvent-accessibility weighting
- local three-dimensional context
- electrostatic weighting
- protonation-dependent weighting
- patch-density calculation
- literature sources supporting residue–chemistry interaction assignments

---

## Security note

InterfaceScout is open source and runs as a local application.

The launcher scripts are not digitally signed, so Windows SmartScreen or some antivirus programs may occasionally display a warning for downloaded scripts.

To proceed:

- **SmartScreen:** click **More info → Run anyway**
- **Downloaded ZIP blocked:** right-click the ZIP → **Properties** → **Unblock** → **OK**
- **Antivirus quarantined a file:** inspect and restore the file if appropriate

All launcher scripts are plain-text files and can be inspected directly.

---

## Troubleshooting

- **"could not find backend\main.py"**  
  Keep `run_local.bat` / `start.bat` next to the `backend/` folder.

- **"No environment found"**  
  Run the setup script once first: `run_local.bat` on Windows or `run_local.sh` on macOS/Linux.

- **Browser didn't open**  
  Open `http://localhost:8000` manually.

- **APBS not found**  
  Check that APBS is installed and accessible to InterfaceScout.

- **Port 8000 busy**  
  Close the process using the port or restart the launcher.

---

## Citation

If you use InterfaceScout in your research, please cite the associated InterfaceScout publication.

Software repository:

`https://github.com/zeynepguneryilmaz/InterfaceScout`

---

*InterfaceScout · protein surface chemistry · residue-level compatibility mapping · Shrake–Rupley SASA · PDB2PQR · PROPKA · APBS electrostatics · surface-chemistry interaction propensity · patch-density mapping*
