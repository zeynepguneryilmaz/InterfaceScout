# InterfaceScout

**Residue-level protein surface chemistry & functional-group compatibility mapping.**

InterfaceScout analyzes the solvent-exposed surface of a protein once and scores
every residue for compatibility with a panel of generic surface chemistries
(cationic, anionic, hydrogen-bond donor/acceptor, π/graphitic, hydrophobic,
oxide, hydroxyapatite/Ca²⁺, metal-coordination, gold, phosphate). Scoring is
fully deterministic and physics-based — literature-derived interaction energies
weighted by solvent exposure, three-dimensional cluster context, APBS
electrostatics, and a pH-dependent protonation factor. **No machine learning.**

The app runs entirely on your computer and opens in your web browser at
`http://localhost:8000`. Nothing is uploaded to any server.

---

## What's in this folder

```
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
```

**Keep all of these together in one folder.** The launchers sit next to
`backend/` and find it automatically — don't move them into subfolders.

---

## Requirements

- **Python 3.11 or 3.12** (with SSL). On Windows the setup script finds it
  automatically; on macOS/Linux install from python.org or your package manager
  if needed.
- An internet connection **the first time only** (to download dependencies and,
  on Windows, the APBS electrostatics binary).
- Works on **Windows 10/11**, **macOS** (Intel or Apple Silicon), and **Linux**.

---

## Windows

1. **Double-click `run_local.bat`** (first time only). It installs everything,
   puts an **InterfaceScout** icon on your Desktop, and starts the app.
2. **From then on, double-click the InterfaceScout Desktop icon.** It detects the
   existing installation and starts immediately (no re-install); your browser
   opens at `http://localhost:8000`.

> The Desktop icon runs `run_local.bat`, which is safe to run any time: the first
> run sets everything up, later runs start the app straight away. To force a
> clean reinstall, delete the `backend\.venv` folder and run it again.
> If Windows SmartScreen warns about the `.bat` file, choose **More info → Run
> anyway**. If the app fails to start, the window stays open with the error and a
> log is written to `backend\startup_log.txt`.

---

## macOS

1. **Run `run_local.sh` once**: right-click it → **Open**, or in Terminal run
   `bash run_local.sh`. It installs everything and puts an
   **InterfaceScout.command** icon on your Desktop.
2. **From then on, double-click `InterfaceScout.command`** on your Desktop. The
   backend starts in the background and your browser opens; no Terminal stays
   open.

> Gatekeeper may ask for confirmation the first time — right-click → **Open**,
> then confirm. To give the icon the app picture: right-click → **Get Info**,
> drag `interfacescout.png` onto the icon at the top-left.

---

## Linux

1. **Run `run_local.sh` once**: `bash run_local.sh`. It installs everything and
   puts an **InterfaceScout.desktop** launcher on your Desktop.
2. **From then on, double-click the InterfaceScout Desktop launcher.** The
   backend starts with no terminal window and your browser opens.

> If your desktop marks the launcher "untrusted", right-click → **Allow
> Launching** (the installer marks it trusted on GNOME automatically).

---

## Using the app

1. Enter a PDB ID (e.g. `4F5U`) and click **Fetch**, or paste/upload a structure.
2. Set the environment: **pH**, **ionic strength** (mM), **temperature** (K), and
   **patch radius** (Å, default 12).
3. Click **Analyze Surface**. The protein is analyzed once; every surface
   chemistry is then available instantly.
4. Explore the 3D map, switch between **propensity** and **patch-density**
   colouring, and read the per-residue table.
5. Export results as **CSV**, **PDB**, or **PDF**.

The **Theory** page in the app documents the scoring scheme, every equation, and
the literature behind each interaction energy and the pH-dependent protonation
factor.

---

## Security note (antivirus / SmartScreen warnings)

InterfaceScout is **open source and runs entirely on your own computer** — it
does not upload data anywhere. Because the launcher scripts are not digitally
signed (code-signing certificates are costly and unusual for academic tools),
Windows SmartScreen or some antivirus products (Defender, Avast, Kaspersky,
McAfee) may show a one-time **"unrecognized app"** warning. This is a *false
positive* caused by the missing signature, not by any harmful content — every
script is plain text you can open and read.

To proceed:

- **SmartScreen:** click **More info → Run anyway**.
- **Downloaded ZIP blocked:** right-click the ZIP → **Properties** → tick
  **Unblock** → **OK**, then extract.
- **Antivirus quarantined a file:** restore it and, if needed, add the
  InterfaceScout folder to your antivirus exclusions. Downloading the release
  from the official GitHub repository (rather than email/USB) also reduces these
  warnings, because files from a known source accumulate reputation over time.

The scripts only: find Python, create a local virtual environment, install the
listed Python packages, download the APBS electrostatics binary (Windows), and
start the local server. Nothing else.

---

## Troubleshooting

- **"could not find backend\main.py":** the launcher was moved out of the
  InterfaceScout folder. Keep `run_local.bat`/`start.bat` next to `backend/`.
- **"No environment found":** run the setup script once first
  (`run_local.bat` on Windows, `run_local.sh` on macOS/Linux).
- **Browser didn't open:** open `http://localhost:8000` manually.
- **APBS not found (Windows):** the setup script installs it (via conda if
  available, otherwise a standalone download). Without APBS the app falls back to
  a Debye-Huckel electrostatic estimate.
- **Port 8000 busy:** the launchers free it automatically.

---

*InterfaceScout · residue-level surface chemistry mapping · Shrake-Rupley SASA ·
APBS electrostatics · literature-derived interaction energies.*
