"""
Protein Feature Affinity Analyzer — Backend v4
================================================
Stack  : FastAPI + pdb2pqr + APBS + biopython + numpy

Model  : Material-independent FEATURE AFFINITY MAPPING.
         The protein is analyzed once; each surface patch is mapped against
         all surface features and functional groups. No material is selected.

Structural analysis (per patch):
         - Shrake-Rupley SASA      (biopython, probe=1.4Å, n=100)
         - Henderson-Hasselbalch   (pH-dependent residue charges + PropKa)
         - APBS LinearPB           (pdb2pqr→PQR, APBS→.dx, trilinear φ)
                                     Fallback: Debye-Huckel if APBS unavailable
         - DSSP secondary structure (helix / sheet / loop)
         - PCA curvature           (convex / flat / concave)
         - Binary-search pI        (patch isoelectric point)

Feature affinity (per patch):
         - functional_group_affinity : ΔG for each surface group
                                       (−SH/Au, Ca2+, Ti-O, π, −NH2, −COOH, …)
                                       summed per-residue, SASA-weighted
         - surface_property_affinity : charge/hydrophobicity/H-bond/pI profile
         - patch_intrinsic_penalties : orientational entropy + denaturation cost

Output : ranked functional-group recommendations + property profile per patch.
"""

import os, sys, math, json, uuid, shutil, subprocess, tempfile, logging
from pathlib import Path
from typing import Optional, List, Dict, Any

import numpy as np
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

from Bio.PDB import PDBParser, DSSP
from Bio.PDB.SASA import ShrakeRupley

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
log = logging.getLogger("ppia")

# ── Binary discovery ────────────────────────────────────────────────────────
def _find(name: str) -> Optional[str]:
    # 1. Explicit environment override (e.g. APBS_PATH, PDB2PQR_PATH)
    env_key = f"{name.upper().replace('-','_')}_PATH"
    env_val = os.environ.get(env_key)
    if env_val and Path(env_val).exists():
        return env_val

    # 2. On PATH. A conda-installed APBS often lives under
    #    ...\Library\bin\apbs.exe (even inside a Schrodinger/PyMOL conda env) and
    #    is fully functional. Only skip an APBS that is bundled *inside* the
    #    PyMOL application itself (not a conda Library\bin install).
    p = shutil.which(name)
    if p:
        low = p.lower().replace("/", "\\")
        if name == "apbs":
            is_conda_libbin = "\\library\\bin\\" in low
            looks_bundled = ("pymol" in low or "schrodinger" in low) and not is_conda_libbin
            if looks_bundled:
                p = None
        if p:
            return p

    # 3. apbs_binary python package (Linux/macOS)
    if name == "apbs":
        try:
            from apbs_binary import APBS_BIN_PATH
            if Path(str(APBS_BIN_PATH)).exists():
                return str(APBS_BIN_PATH)
        except ImportError:
            pass

    # 4. Deep search only for APBS (the one binary we really want and that
    #    users commonly install manually). For pdb2pqr/mkdssp, PATH is enough
    #    (mkdssp is optional — there is a geometry fallback). This keeps
    #    startup fast and avoids scanning large folders for things not present.
    if name != "apbs":
        return None

    exe = name + (".exe" if os.name == "nt" else "")
    here = Path(__file__).parent
    roots = [
        here, here / "bin", here / "vendor",
        here.parent,                       # project root
        here.parent / "apbs-win",          # where run_local.bat extracts APBS
    ]
    home = Path.home()
    if os.name == "nt":
        roots += [home / "Downloads"]
    else:
        roots += [Path("/usr/local/bin"), Path("/opt")]

    def _scan(base: Path, max_depth: int = 4):
        """Yield files named `exe` under base, up to max_depth, cheaply."""
        base = base.resolve()
        try:
            stack = [(base, 0)]
            while stack:
                d, depth = stack.pop()
                try:
                    with os.scandir(d) as it:
                        for e in it:
                            if e.is_file() and e.name.lower() == exe.lower():
                                yield Path(e.path)
                            elif e.is_dir() and depth < max_depth:
                                # skip obviously irrelevant huge dirs
                                nm = e.name.lower()
                                if nm in ("node_modules", "__pycache__", ".git",
                                          "windows", "winsxs", "$recycle.bin"):
                                    continue
                                stack.append((Path(e.path), depth + 1))
                except (PermissionError, OSError):
                    continue
        except (PermissionError, OSError):
            return

    candidates = []
    for root in roots:
        if root and root.exists():
            candidates.extend(_scan(root))
        if candidates and name != "apbs":
            # for non-apbs, first hit is fine
            return str(candidates[0])

    if candidates:
        # For APBS prefer a path that contains 'bin' (the real solver) over
        # helper exes under tools/visualization/opendx etc.
        for c in candidates:
            if os.sep + "bin" + os.sep in str(c).lower():
                return str(c)
        return str(candidates[0])
    return None

PDB2PQR = _find("pdb2pqr")
APBS    = _find("apbs")
MKDSSP  = _find("mkdssp") or _find("dssp")
log.info(f"pdb2pqr={PDB2PQR}  apbs={APBS}  dssp={MKDSSP}")

app = FastAPI(title="InterfaceScout Backend", version="4.0")
app.add_middleware(CORSMiddleware, allow_origins=["*"],
                   allow_methods=["*"], allow_headers=["*"])

# ════════════════════════════════════════════════════════════════════════════
# CONSTANTS & DATABASES
# ════════════════════════════════════════════════════════════════════════════
KB  = 1.380649e-23
NA  = 6.02214076e23
E   = 1.602176634e-19
EPS0 = 8.8541878e-12

# AA: hydro, pKa, charge_sign(+1=base/-1=acid), hbD, hbA, sigC, sigN, sigO, sigS, category
AA = {
    "ALA": ( 0.25, None,  0, 0, 1, -22.5, 14.0,  6.0,    0, "hydrophobic"),
    "ARG": (-1.80, 12.5, +1, 5, 2, -22.5, 14.0,  6.0,    0, "pos_charged"),
    "ASN": (-0.64, None,  0, 2, 3, -22.5, 14.0,  6.0,    0, "polar"),
    "ASP": (-0.72,  3.9, -1, 1, 4, -22.5, 14.0,  6.0,    0, "neg_charged"),
    "CYS": ( 0.04,  8.3, -1, 1, 2, -22.5, 14.0,  6.0, -17.0,"polar"),
    "GLN": (-0.69, None,  0, 2, 3, -22.5, 14.0,  6.0,    0, "polar"),
    "GLU": (-0.62,  4.1, -1, 1, 4, -22.5, 14.0,  6.0,    0, "neg_charged"),
    "GLY": ( 0.16, None,  0, 1, 1, -22.5, 14.0,  6.0,    0, "special"),
    "HIS": (-0.40,  6.5, +1, 2, 3,  -9.0, 14.0,  6.0,    0, "polar"),
    "ILE": ( 0.73, None,  0, 0, 1, -22.5, 14.0,  6.0,    0, "hydrophobic"),
    "LEU": ( 0.53, None,  0, 0, 1, -22.5, 14.0,  6.0,    0, "hydrophobic"),
    "LYS": (-1.10, 10.5, +1, 3, 2, -22.5, 14.0,  6.0,    0, "pos_charged"),
    "MET": ( 0.26, None,  0, 0, 1, -22.5, 14.0,  6.0, -17.0,"hydrophobic"),
    "PHE": ( 0.61, None,  0, 0, 1,  -9.0, 14.0,  6.0,    0, "aromatic"),
    "PRO": (-0.07, None,  0, 0, 1, -22.5, 14.0,  6.0,    0, "special"),
    "SER": (-0.26, None,  0, 2, 2, -22.5, 14.0,  6.0,    0, "polar"),
    "THR": (-0.18, None,  0, 2, 2, -22.5, 14.0,  6.0,    0, "polar"),
    "TRP": ( 0.37, None,  0, 1, 2,  -9.0, 14.0,  6.0,    0, "aromatic"),
    "TYR": ( 0.02, 10.1, -1, 2, 2,  -9.0, 14.0,  6.0,    0, "aromatic"),
    "VAL": ( 0.54, None,  0, 0, 1, -22.5, 14.0,  6.0,    0, "hydrophobic"),
}

HAMAKER = {           # × 10⁻²¹ J, protein-water-material
    "titanium": 40, "hap": 15, "gold": 100, "pdms": 5,
    "graphene": 50, "plga": 3, "chitosan": 8, "zirconia": 20,
}

# ── L14: PEG brush parameters per material preset ─────────────────────────
# (L_brush_Å, sigma_chains_per_nm2)
# L = Flory radius estimate: RF ≈ a·N^0.6, a=3.5Å, N≈45 for PEG-2000
# sigma: grafting density (0 = no brush)
PEG_BRUSH: dict = {
    "plga":     (55.0, 0.20),   # PEG-PLGA ~2000 Da brush
    "gold":     (25.0, 0.08),   # thiol-PEG SAM (sparse)
    "chitosan": ( 0.0, 0.00),
    "titanium": ( 0.0, 0.00),
    "hap":      ( 0.0, 0.00),
    "pdms":     ( 0.0, 0.00),
    "graphene": ( 0.0, 0.00),
    "zirconia": ( 0.0, 0.00),
}

# ── L13: Dielectric constants ─────────────────────────────────────────────
EPS_PROTEIN = 2.0    # protein interior (Å-scale)
EPS_WATER   = 78.54  # bulk water at 310 K

# ── L15: B-factor scale ───────────────────────────────────────────────────
BFAC_SCALE = 0.40    # fraction of denaturation cost weighted by B-factor

# ════════════════════════════════════════════════════════════════════════════
# PYDANTIC MODELS
# ════════════════════════════════════════════════════════════════════════════
class MaterialParams(BaseModel):
    charge:         float = -10.0   # zeta potential mV
    surface_energy: float = 45.0    # mJ/m²
    hydrophilicity: float = 0.6     # 0–1
    roughness:      float = 5.0     # nm Ra
    dielectric:     float = 10.0    # ε_r
    functional:     str   = "-OH"
    preset_key:     str   = ""

class EnvParams(BaseModel):
    pH:    float = 7.4
    ionic: float = 150.0   # mM
    temp:  float = 310.0   # K
    patch_radius: float = 12.0   # Å, radius for SAP-style patch-density (user-set)

class AnalyzeRequest(BaseModel):
    pdb_id:   Optional[str] = None
    pdb_text: Optional[str] = None
    env:      EnvParams      = EnvParams()
    # NOTE: material input removed — the analyzer now maps each patch against
    # ALL surface features and functional groups, material-independently.

# ════════════════════════════════════════════════════════════════════════════
# LAYER 2 — Henderson-Hasselbalch protonation
# ════════════════════════════════════════════════════════════════════════════
def residue_charge(res_name: str, pH: float) -> float:
    aa = AA.get(res_name)
    if not aa or aa[1] is None: return 0.0
    pKa, sign = aa[1], aa[2]
    if sign == +1:
        return  1.0 / (1.0 + 10**(pH - pKa))
    else:
        return -(1.0 / (1.0 + 10**(pKa - pH)))


def _protonation_fraction(res_name: str, pH: float) -> float:
    """
    Fraction of the residue in its CHARGED protonation state at this pH,
    in [0, 1]. For cationic residues (Lys/Arg/His) this is the protonated
    (positive) fraction; for anionic residues (Asp/Glu/Cys/Tyr) it is the
    deprotonated (negative) fraction. Residues with no ionisable side chain
    return 1.0 (always "available"). Uses Henderson-Hasselbalch with the side
    chain pKa from the AA table.
    """
    aa = AA.get(res_name)
    if not aa or aa[1] is None:
        return 1.0
    pKa, sign = aa[1], aa[2]
    if sign == +1:                      # basic: charged = protonated
        return 1.0 / (1.0 + 10**(pH - pKa))
    else:                               # acidic: charged = deprotonated
        return 1.0 / (1.0 + 10**(pKa - pH))


def _protonation_factor(res_name: str, mechanism: str, pH: float) -> float:
    """
    pH-dependent scaling of a residue's contribution to a surface chemistry,
    based on its interaction MECHANISM and protonation state. Returns a factor
    (typically in [~0.3, 1.0]) that multiplies the literature base energy.

    Literature basis:
      - cation-pi (Lys/Arg/His -> aromatic/graphitic surface): requires a
        POSITIVE charge on the residue. His is the pH switch: protonated His+
        binds aromatic faces strongly, neutral His does not. Lys/Arg stay
        protonated across pH 4.7-7.4. (Dougherty, Science 1996, 271:163;
        Ma & Dougherty, Chem Rev 1997, 97:1303; Gallivan & Dougherty, PNAS
        1999, 96:9459; Liao et al., J Chem Theory Comput 2024, acs.jctc.4c00606)
      - pi-pi stacking (Phe/Tyr/Trp): aromatic rings are neutral and do NOT
        ionise in this range -> pH-independent. (Ma & Dougherty, Chem Rev 1997)
      - carboxylate coordination (Asp/Glu -> oxide / Ca / metal surfaces):
        requires the DEPROTONATED carboxylate; protonation near/below pKa
        weakens binding. (JACS Au 2021, 10.1021/jacsau.1c00565; bioceramics
        review, Phil Trans / PMC3363020)
      - hydrophobic / thiol-gold / plain H-bond donors-acceptors on neutral
        side chains: pH-independent in this range.
    A small floor (0.3) is kept so a residue never contributes exactly zero,
    reflecting residual/transient interactions and avoiding hard cut-offs.
    """
    m = (mechanism or "").lower()
    FLOOR = 0.3

    # --- cation-pi: scale by POSITIVE (protonated) fraction ---
    if "cation" in m:
        frac = _protonation_fraction(res_name, pH)   # positive fraction for basics
        return FLOOR + (1.0 - FLOOR) * frac

    # --- pi-pi stacking on aromatics: pH-independent ---
    if "π-π" in m or "pi-pi" in m or "stacking" in m:
        # His listed as "π-π / cation-π": treat its pH dependence via cation-pi
        # branch above; pure aromatic stacking (Phe/Tyr/Trp) is pH-independent.
        return 1.0

    # --- carboxylate-based coordination to oxide/Ca/metal: need deprotonated ---
    if any(k in m for k in ("coordination", "carboxylate", "ca", "oxide", "metal")):
        aa = AA.get(res_name)
        if aa and aa[1] is not None and aa[2] == -1:   # Asp/Glu (acidic)
            frac = _protonation_fraction(res_name, pH)  # deprotonated fraction
            return FLOOR + (1.0 - FLOOR) * frac
        return 1.0   # His/Cys metal coordination handled as ~pH-stable here

    # default: no pH scaling
    return 1.0

# ════════════════════════════════════════════════════════════════════════════
# LAYER 3 — pdb2pqr + APBS + DX parser
# ════════════════════════════════════════════════════════════════════════════
def run_pdb2pqr(pdb: Path, pqr: Path, pH: float) -> None:
    if not PDB2PQR:
        raise RuntimeError("pdb2pqr not found")
    cmd = [PDB2PQR, "--ff=PARSE",
           f"--with-ph={pH:.2f}",
           "--titration-state-method=propka",
           "--drop-water", "--keep-chain",
           str(pdb), str(pqr)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    if r.returncode != 0:
        raise RuntimeError(f"pdb2pqr failed: {r.stderr[-400:]}")
    if not pqr.exists():
        raise RuntimeError("pdb2pqr produced no PQR file")


def build_apbs_input(pqr: Path, dx_stem: Path, ionic_M: float, T: float) -> Path:
    """
    Builds APBS input using mg-manual with adaptive grid.
    Grid dimensions: nearest odd number >= (protein_extent + 2*pad) / spacing.
    Grid spacing: 0.5 Å; pad: 10 Å; minimum dime: 33.
    """
    inp = pqr.parent / "apbs.in"
    c   = ionic_M

    # Read PQR to determine bounding box
    xs, ys, zs = [], [], []
    for line in pqr.read_text().splitlines():
        if not line.startswith(("ATOM","HETATM")): continue
        try:
            xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
        except: pass
    if not xs:
        raise ValueError("PQR has no atoms")

    pad = 15.0; spacing = 0.5
    def grid_dim(lo, hi):
        n = int(math.ceil((hi - lo + 2*pad) / spacing))
        n = max(n, 33)
        if n % 2 == 0: n += 1   # must be odd
        return n, (lo+hi)/2.0

    nx, cx = grid_dim(min(xs), max(xs))
    ny, cy = grid_dim(min(ys), max(ys))
    nz, cz = grid_dim(min(zs), max(zs))
    gx = nx * spacing; gy = ny * spacing; gz = nz * spacing

    log.info(f"APBS grid: {nx}x{ny}x{nz}, glen {gx:.1f}x{gy:.1f}x{gz:.1f}")

    inp.write_text(f"""read
    mol pqr {pqr}
end
elec
    mg-manual
    mol 1
    dime  {nx} {ny} {nz}
    glen  {gx:.2f} {gy:.2f} {gz:.2f}
    gcent {cx:.3f} {cy:.3f} {cz:.3f}
    lpbe
    bcfl sdh
    ion charge +1 conc {c:.4f} radius 2.0
    ion charge -1 conc {c:.4f} radius 2.0
    pdie 2.0
    sdie 78.54
    srfm smol
    sdens 10.0
    chgm spl2
    srad 1.4
    swin 0.3
    temp {T:.1f}
    calcenergy total
    calcforce no
    write pot dx {dx_stem}
end
quit
""")
    return inp


def _apbs_env() -> dict:
    """
    Build the environment for running APBS.

    Windows: a conda-installed APBS needs its conda environment DLLs
    (Library\bin and the env root). A standalone APBS-3.4.1 build ships its own
    MinGW DLLs in its bin/; for that case we keep a minimal PATH so unrelated
    libgfortran/openblas from other tools cannot shadow it (which would cause
    0xc000007b).

    macOS/Linux: the apbs-binary package needs its LIB_DIR on the loader path.
    """
    env = dict(os.environ)

    if os.name == "nt":
        if APBS:
            bin_dir = str(Path(APBS).parent)
            sysroot = os.environ.get("SystemRoot", r"C:\Windows")
            low = bin_dir.lower().replace("/", "\\")
            paths = [bin_dir]
            if "\\library\\bin" in low:
                # conda layout: ...\<env>\Library\bin\apbs.exe
                # add the sibling conda DLL locations the binary depends on
                libbin = Path(bin_dir)                      # Library\bin
                lib    = libbin.parent                       # Library
                envroot= lib.parent                          # <env> root (has python3x.dll)
                paths += [
                    str(lib / "mingw-w64" / "bin"),
                    str(lib / "usr" / "bin"),
                    str(envroot),                            # python39/310 .dll lives here
                    str(envroot / "Scripts"),
                ]
            paths += [
                sysroot + r"\System32", sysroot, sysroot + r"\System32\Wbem",
            ]
            env["PATH"] = os.pathsep.join(p for p in paths if p)
        return env

    # macOS / Linux
    try:
        import apbs_binary
        lib_dir = str(getattr(apbs_binary, "LIB_DIR", ""))
        if lib_dir:
            if sys.platform == "darwin":
                env["DYLD_LIBRARY_PATH"] = lib_dir + os.pathsep + env.get("DYLD_LIBRARY_PATH", "")
            else:
                env["LD_LIBRARY_PATH"] = lib_dir + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    except Exception:
        pass
    if APBS:
        bin_dir = str(Path(APBS).parent)
        env["PATH"] = bin_dir + os.pathsep + env.get("PATH", "")
    return env


def run_apbs(inp: Path, workdir: Path) -> Path:
    if not APBS:
        raise RuntimeError("APBS not found")
    r = subprocess.run([APBS, str(inp)],
                       capture_output=True, text=True,
                       cwd=str(workdir), timeout=600,
                       env=_apbs_env())
    dxfiles = list(workdir.glob("*.dx"))
    if not dxfiles:
        # exit code 3221225595 (0xC000007B) = DLL/bitness mismatch on Windows
        hint = ""
        if r.returncode in (3221225595, -1073741701):
            hint = (" [0xC000007B: APBS could not load its runtime DLLs. "
                    "This usually means another libgfortran/openblas on PATH "
                    "is shadowing APBS's own. The minimal-PATH launcher should "
                    "prevent this.]")
        log.warning(f"APBS stdout: {r.stdout[-200:]}")
        log.warning(f"APBS stderr: {r.stderr[-200:]}")
        raise RuntimeError(f"APBS produced no .dx (exit={r.returncode}){hint}")
    return dxfiles[0]


def parse_dx(dx: Path) -> dict:
    lines = dx.read_text().splitlines()
    origin = delta = counts = None
    data: list[float] = []
    i = 0
    while i < len(lines):
        ln = lines[i].strip()
        if ln.startswith("object 1 class gridpositions counts"):
            p = ln.split(); counts = [int(p[-3]), int(p[-2]), int(p[-1])]
        elif ln.startswith("origin"):
            origin = list(map(float, ln.split()[1:4]))
        elif ln.startswith("delta"):
            if delta is None: delta = []
            delta.append(list(map(float, ln.split()[1:4])))
        elif ln.startswith("object 3"):
            i += 1
            while i < len(lines):
                row = lines[i].strip()
                if row.startswith(("attribute", "object")): break
                data.extend(map(float, row.split()))
                i += 1
            continue
        i += 1
    if origin is None or delta is None or counts is None:
        raise ValueError("DX parse failed — missing header fields")
    return {"origin": origin, "delta": delta, "counts": counts,
            "data": np.array(data, dtype=np.float32)}


def phi_at(grid: dict, x: float, y: float, z: float) -> float:
    """Trilinear interpolation of φ (kT/e) at position (x,y,z) Å."""
    o = grid["origin"]; d = grid["delta"]; n = grid["counts"]
    dat = grid["data"]
    dx = d[0][0]; dy = d[1][1]; dz = d[2][2]
    fx = (x - o[0]) / dx
    fy = (y - o[1]) / dy
    fz = (z - o[2]) / dz
    ix, iy, iz = int(fx), int(fy), int(fz)
    nx, ny, nz = n
    if ix < 0 or iy < 0 or iz < 0 or ix >= nx-1 or iy >= ny-1 or iz >= nz-1:
        return 0.0
    tx, ty, tz = fx-ix, fy-iy, fz-iz
    def idx(a,b,c): return a*ny*nz + b*nz + c
    return float(
        dat[idx(ix,  iy,  iz  )]*(1-tx)*(1-ty)*(1-tz) +
        dat[idx(ix+1,iy,  iz  )]*tx    *(1-ty)*(1-tz) +
        dat[idx(ix,  iy+1,iz  )]*(1-tx)*ty    *(1-tz) +
        dat[idx(ix,  iy,  iz+1)]*(1-tx)*(1-ty)*tz     +
        dat[idx(ix+1,iy+1,iz  )]*tx    *ty    *(1-tz) +
        dat[idx(ix+1,iy,  iz+1)]*tx    *(1-ty)*tz     +
        dat[idx(ix,  iy+1,iz+1)]*(1-tx)*ty    *tz     +
        dat[idx(ix+1,iy+1,iz+1)]*tx    *ty    *tz
    )


def dG_PB(patch_residues: list, atoms_by_res: dict, grid: dict, T: float) -> dict:
    """
    ΔG_elec = Σᵢ qᵢ·φ(rᵢ)/2   [kT → kcal/mol]
    Factor /2: avoids double counting of self-energy.
    """
    kT_kcal = KB * T * NA / 4184.0
    total = 0.0
    net_q  = 0.0
    for res in patch_residues:
        atoms = atoms_by_res.get(res["key"], [])
        for a in atoms:
            q   = a.get("charge_pqr", 0.0)
            phi = phi_at(grid, a["x"], a["y"], a["z"])
            total += q * phi
            net_q += q
    dG = total * kT_kcal / 2.0
    return {"dG_kcal": round(dG, 5), "net_charge_pqr": round(net_q, 4)}


def load_pqr_charges(pqr: Path, atoms_by_res: dict) -> None:
    for line in pqr.read_text().splitlines():
        if not line.startswith(("ATOM","HETATM")): continue
        try:
            aname = line[12:16].strip()
            cid   = line[21:22].strip() or "A"
            rid   = int(line[22:26].strip())
            q     = float(line[54:62].strip())
            key   = f"{cid}_{rid}"
            for a in atoms_by_res.get(key, []):
                if a["atomName"] == aname:
                    a["charge_pqr"] = q; break
        except (ValueError, IndexError):
            continue

# ════════════════════════════════════════════════════════════════════════════
# LAYER 4 — Eisenberg-Weiss ASP hydrophobic ΔG
# ════════════════════════════════════════════════════════════════════════════
def dG_hydrophobic(residues: list, mat: MaterialParams) -> float:
    mh = 1.0 - mat.hydrophilicity
    total = 0.0
    for r in residues:
        aa = AA.get(r["res_name"])
        if not aa: continue
        sig = 0.50*aa[5] + 0.20*aa[6] + 0.20*aa[7] + 0.10*aa[8]
        total += sig * r.get("sasa", 0.0) * mh
    return round(total / 1000.0, 5)   # cal → kcal

# ════════════════════════════════════════════════════════════════════════════
# LAYER 5 — H-bond + π-π
# ════════════════════════════════════════════════════════════════════════════
def parse_func_groups(s: str) -> dict:
    u = s.upper()
    d = a = pi = 0
    if "-OH"   in u: d+=1; a+=1
    if "NH2"   in u: d+=2
    if "NH"    in u and "NH2" not in u: d+=1
    if "COOH"  in u: d+=1; a+=2
    if "-SH"   in u: d+=1; a+=1
    if "PO4"   in u or "-PO" in u: a+=3
    if "π"     in u or "PI" in u:  pi+=1
    if "TIO"   in u or "ZRO" in u: a+=2
    if "-AU"   in u: a+=1
    return {"donor": d, "acceptor": a, "pi": pi}

def dG_hbond(hbD: int, hbA: int, ar_frac: float, mat: MaterialParams) -> float:
    mg  = parse_func_groups(mat.functional)
    hb  = min(hbD, mg["acceptor"]*2) + min(hbA, mg["donor"]*2)
    pip = ar_frac * mg["pi"] * -3.0
    return round(hb * (-0.8) + pip, 5)

# ════════════════════════════════════════════════════════════════════════════
# LAYER 6 — Hamaker vdW
# ════════════════════════════════════════════════════════════════════════════
def dG_vdw(area_A2: float, mat: MaterialParams) -> float:
    A   = HAMAKER.get(mat.preset_key.lower(), 10) * 1e-21   # J
    D   = 3.0e-10           # 3 Å contact distance
    ar  = area_A2 * 1e-20   # Å² → m²
    dG  = (-A / (6*math.pi*D*D)) * ar / 4184.0
    return round(dG, 5)

# ════════════════════════════════════════════════════════════════════════════
# LAYER 7 — Solvation entropy (water release)
# ════════════════════════════════════════════════════════════════════════════
def dG_solvation(area_A2: float, mat: MaterialParams, T: float) -> float:
    n_w = area_A2 / 12.0
    mh  = 1.0 - mat.hydrophilicity
    dS  = mh*0.35 + (1-mh)*0.10   # cal/mol/K per water
    return round(-(T * n_w * dS) / 1000.0, 5)

# ════════════════════════════════════════════════════════════════════════════
# LAYER 8 — Counterion release entropy
# ════════════════════════════════════════════════════════════════════════════
def dG_counterion(net_q: float, ionic_mM: float, T: float) -> float:
    n = abs(net_q)
    if n < 0.05: return 0.0
    V = max(1.0, 1000.0/ionic_mM)
    return round(-(n * KB * T * NA * math.log(V)) / 4184.0, 5)

# ════════════════════════════════════════════════════════════════════════════
# LAYER 9 — Denaturation risk
# ════════════════════════════════════════════════════════════════════════════
def dG_denaturation(ss: dict, n_res: int, mat: MaterialParams) -> float:
    mh   = 1.0 - mat.hydrophilicity
    cost = ((ss.get("H",0)*2.0 + ss.get("E",0)*1.5 + ss.get("C",0)*0.3)
            / max(n_res, 1))
    return round(cost * mh * 0.25, 5)

# ════════════════════════════════════════════════════════════════════════════
# LAYER 10 — Curvature (PCA of CA coords)
# ════════════════════════════════════════════════════════════════════════════
def curvature(coords: list) -> dict:
    if len(coords) < 3:
        return {"index": 0, "variance": 30.0}
    arr = np.array(coords)
    var = float(np.mean(np.sum((arr - arr.mean(0))**2, axis=1)))
    idx = 1 if var > 60 else (-1 if var < 25 else 0)
    return {"index": idx, "variance": round(var, 2)}

def curvature_mod(curv: dict, mat: MaterialParams) -> float:
    r = mat.roughness
    if curv["index"] ==  1: return round(1.0 + r/150.0, 4)
    if curv["index"] == -1: return 0.75
    return 1.0

# ════════════════════════════════════════════════════════════════════════════
# DEBYE-HUCKEL FALLBACK (Layer 3 when APBS unavailable)
# ════════════════════════════════════════════════════════════════════════════
def dG_DH(net_charge: float, dipole: float,
          mat: MaterialParams, ionic_mM: float, T: float) -> dict:
    epsr = mat.dielectric
    I_m3 = (ionic_mM/1000) * 1000
    k2   = (2*NA*E*E*I_m3) / (EPS0*epsr*KB*T)
    lD   = (1/math.sqrt(k2)) * 1e10   # Å
    D    = 5.0
    eps  = EPS0 * epsr
    qp   = net_charge * E
    qm   = (mat.charge/1000) * E
    cc   = (qp*qm) / (4*math.pi*eps*D*1e-10) * math.exp(-D/lD) / (NA*4184)
    dp   = dipole * abs(mat.charge/100) * 0.005 * math.exp(-D/lD)
    return {"dG_kcal": round(cc+dp, 5), "debye_len": round(lD,2),
            "method": "DebyeHuckel_fallback"}

# ════════════════════════════════════════════════════════════════════════════
# LAYER 11 — Orientational + Translational Entropy Loss
# ════════════════════════════════════════════════════════════════════════════
def dG_orientational(n_res: int, T: float) -> float:
    """
    When a protein patch adsorbs it loses:
      • 3 translational DoF  → ΔS_trans = −kB·ln(V_site/V_bulk)
      • 3 rotational DoF     → ΔS_rot   = −kB·ln(8π²/Δω³)

    Canonical approximation (Finkelstein & Janin, 1989):
        ΔG_orient = +kB·T·NA · [ln(8π²/Δω³) + ln(V_site/V_bulk)] / 4184
    Parameters (typical):
        Δω = 0.30 rad  (angular tolerance on a surface)
        V_bulk / V_site ≈ 55.5 M / 1 mM active site ≈ 55,500   → translational
    Protein size correction: larger patch = more constrained → steeper loss.
    """
    kT_kcal = KB * T * NA / 4184.0
    delta_omega = 0.30   # rad, angular freedom on surface (empirical)
    # Rotational: 3 angles each constrained to Δω
    dG_rot = kT_kcal * math.log(8 * math.pi**2 / delta_omega**3)
    # Translational: protein confined from bulk (55 M) to surface layer (~1 nm thick)
    # ΔG_trans ≈ kT·ln(c_bulk_sites/c_surface) — simplified to fixed term
    dG_trans = kT_kcal * math.log(55.5e3 / 1.0)   # ~5.5 kcal/mol baseline
    # Size correction: larger patch loses more orientational freedom
    # Scale ~ ln(n_res) because each residue-surface contact restricts one DoF
    size_factor = 1.0 + 0.05 * math.log(max(n_res, 1))
    return round((dG_rot + dG_trans) * size_factor * 0.10, 5)
    # ×0.10: partial loss factor (only fraction of full DoF lost on adsorption)


# ════════════════════════════════════════════════════════════════════════════
# LAYER 12 — Isoelectric point (pI) helper + pI-surface charge mismatch
# ════════════════════════════════════════════════════════════════════════════
def _compute_pI(res_names: list) -> float:
    """Binary search for pH where Σ q_i(pH) = 0."""
    lo, hi = 0.0, 14.0
    for _ in range(60):   # 60 iterations → ±0.0001 pH accuracy
        mid = (lo + hi) / 2.0
        q   = sum(residue_charge(r, mid) for r in res_names)
        if q > 0: lo = mid
        else:     hi = mid
    return (lo + hi) / 2.0


def dG_pI_mismatch(pI: float, pH: float, mat_charge_mV: float,
                   net_charge: float, T: float) -> float:
    """
    The further pH is from pI, the stronger the net protein charge,
    and the stronger the electrostatic attraction or repulsion with
    the surface. This is already partially captured in L3, but pI
    gives an explicit "charge-matching bonus/penalty" independent of
    the PB grid.

    Sign logic:
      • protein net (+), mat_charge (−) → attractive → negative ΔG_pI
      • protein net (+), mat_charge (+) → repulsive  → positive ΔG_pI

    Magnitude ∝ |pH − pI| (how charged the protein is)
              × |mat_charge_mV| / 100 (how charged the surface is)
    """
    kT_kcal = KB * T * NA / 4184.0
    charge_strength  = abs(pH - pI) * 0.20   # 0…~2.8 on scale
    surface_strength = abs(mat_charge_mV) / 50.0
    # Determine sign: attractive (−) if opposite charges
    protein_sign = 1.0 if net_charge > 0 else -1.0
    surface_sign = -1.0 if mat_charge_mV < 0 else 1.0
    sign = protein_sign * surface_sign   # −1 = attractive, +1 = repulsive
    dG = sign * charge_strength * surface_strength * kT_kcal * 0.5
    return round(dG, 5)


# ════════════════════════════════════════════════════════════════════════════
# LAYER 13 — Image Charge / Born Energy (dielectric discontinuity)
# ════════════════════════════════════════════════════════════════════════════
def dG_image_charge(net_charge: float, area_A2: float,
                    mat_dielectric: float, D_Å: float = 5.0) -> float:
    """
    At the protein–water–material interface the dielectric discontinuity
    (ε_protein ≈ 2, ε_water ≈ 80) induces image charges that repel the
    approaching protein.  For a charge q at distance D from a planar
    interface between ε₁ (protein) and ε₂ (material):

        ΔG_image = q² / (4πε₀) × (ε₁ − ε₂) / (ε₁ + ε₂) × 1/(4D)
                                    ↑ reflection coefficient k

    Here we approximate over the whole patch using net_charge.
    Material dielectric enters via the protein-material contrast term.
    Result in kcal/mol.

    Sign: if ε_mat > ε_protein (e.g. water-like material, ε>2) → k negative
          → image charge ATTRACTS → favourable (negative ΔG).
          Dry hydrophobic materials (ε≈2) → k≈0, no image effect.
    """
    if abs(net_charge) < 0.05: return 0.0
    eps1 = EPS_PROTEIN          # protein interior
    eps2 = max(mat_dielectric, 1.0)
    k    = (eps1 - eps2) / (eps1 + eps2)   # reflection coefficient (−1…+1)
    D_m  = D_Å * 1e-10
    q    = net_charge * E
    # Energy of one charge q near the interface
    dG_J = (q**2 / (4 * math.pi * EPS0)) * k / (4 * D_m)
    dG_per_mol = dG_J * NA / 4184.0  # kcal/mol
    # Scale by sqrt(area) as proxy for number of charges in contact zone
    area_factor = math.sqrt(area_A2 / 100.0)   # normalised to 10 Å radius patch
    return round(dG_per_mol * area_factor, 5)


# ════════════════════════════════════════════════════════════════════════════
# LAYER 14 — Steric Barrier (Alexander–de Gennes PEG brush)
# ════════════════════════════════════════════════════════════════════════════
def dG_steric_peg(mat: MaterialParams, T: float) -> float:
    """
    Alexander–de Gennes brush free energy for protein approaching
    a polymer-grafted surface:

        G(D) = kT·s⁻³ · L⁵/³·D⁻⁵/³ / (5/3)    D < L  (compression)
             + kT·s⁻³ · L⁸/³·D⁻²/³ / (8/3)      ← osmotic term

    Here we use the total free energy at contact (D → protein radius ≈ 20 Å):

        ΔG_steric ≈ kT × (L/D)^(5/3) / (sigma × s²)

    Parameters:
        L     = brush height (Å) from PEG_BRUSH dict
        sigma = grafting density (chains/nm²)
        s     = mean distance between grafting points = 1/√sigma  (nm → Å)
        D     = protein approach distance ≈ 20 Å (small globular protein radius)

    If the material has no PEG brush, ΔG_steric = 0.
    """
    key = mat.preset_key.lower()
    L_Å, sigma = PEG_BRUSH.get(key, (0.0, 0.0))
    if L_Å < 1.0 or sigma < 1e-6:
        return 0.0
    kT_kcal = KB * T * NA / 4184.0
    s_nm  = 1.0 / math.sqrt(sigma)          # grafting spacing (nm)
    s_Å   = s_nm * 10.0                     # convert to Å
    D_Å   = 20.0                            # protein contact distance
    if D_Å >= L_Å:
        return 0.0   # protein doesn't compress brush
    # Compression free energy per chain (de Gennes)
    ratio = L_Å / D_Å
    dG_chain = kT_kcal * ((ratio**(5.0/3.0))/(5.0/3.0) +
                           (ratio**(8.0/3.0))/(8.0/3.0))
    # Chains per patch contact area
    n_chains = (100.0 / (s_Å**2))           # per 100 Å² contact
    return round(+dG_chain * n_chains, 5)   # positive = repulsive


# ════════════════════════════════════════════════════════════════════════════
# LAYER 15 — B-factor weighted denaturation (local flexibility)
# ════════════════════════════════════════════════════════════════════════════
def dG_bfactor_flex(residues: list, mat: MaterialParams) -> float:
    """
    PDB isotropic B-factor (Å²) is proportional to mean-square displacement:
        <u²> = 3·B / (8π²)

    High B-factor residues are thermally disordered and lose more
    conformational entropy upon adsorption. We weight the L9 denaturation
    cost by normalised B-factor:

        ΔG_flex = BFAC_SCALE × Σᵢ (B_i / B_ref) × ΔG_denat_i

    B_ref = 30 Å² (typical folded protein surface)
    ΔG_denat_i per-residue: helix=2.0, sheet=1.5, loop=0.3 kcal/mol
    Mat hydrophobicity multiplier same as L9.
    """
    B_REF  = 30.0   # Å², typical folded surface residue
    SS_COST = {"H": 2.0, "E": 1.5, "C": 0.3}
    mh = 1.0 - mat.hydrophilicity
    total = 0.0
    for r in residues:
        B     = r.get("bfac", B_REF)
        w     = min(B / B_REF, 3.0)   # cap at 3× to avoid outlier dominance
        cost  = SS_COST.get(r.get("ss","C"), 0.3)
        total += w * cost
    # average × mat factor × scale constant
    n   = max(len(residues), 1)
    dG  = (total / n) * mh * BFAC_SCALE
    return round(dG, 5)


# ════════════════════════════════════════════════════════════════════════════
# FUNCTIONAL GROUP ↔ RESIDUE INTERACTION DATABASE
# ════════════════════════════════════════════════════════════════════════════
# Material-independent. Each surface functional group interacts with specific
# residue types through a defined physical mechanism. Energies in kcal/mol per
# contacting residue (literature-derived; negative = attractive).
#
# ════════════════════════════════════════════════════════════════════════════
# LITERATURE REFERENCES for the interaction energies
# ════════════════════════════════════════════════════════════════════════════
# Each per-residue energy in FG_INTERACTIONS is grounded in the peer-reviewed
# literature. The citation keys below are surfaced in the UI (Theory page) and
# in generated reports so every number is traceable.
REFERENCES = {
    "AuS": ("Au–S dative/semi-covalent bond ≈ 40–50 kcal/mol. "
            "Hakkinen, Nat. Chem. 2012, 4, 443; "
            "Vericat et al., Chem. Soc. Rev. 2010, 39, 1805; "
            "Pensa et al., Acc. Chem. Res. 2012, 45, 1183."),
    "cation_pi": ("Cation–π binding ≈ 4–7 kcal/mol (up to ~9 kcal/mol in proteins). "
                  "Dougherty, Science 1996, 271, 163; "
                  "Gallivan & Dougherty, PNAS 1999, 96, 9459; "
                  "Dougherty, Chem. Rev. 2025 (cation-π review)."),
    "pi_pi": ("π–π stacking: Phe–Phe ≈ −3.3, Phe–Trp ≈ −4.2, Trp–Trp ≈ −5.2 kcal/mol (QM, in water). "
              "Aromatic residues re-evaluation, J. Phys. Chem. B 2024, 128, 9509 (PMC11403661)."),
    "Ca_carboxylate": ("Carboxylate–Ca²⁺ coordination at hydroxyapatite; Asp/Glu-rich and "
                       "phosphorylated motifs dominate binding. "
                       "Kawasaki, J. Chromatogr. 1991; Goobes et al., PNAS 2006; "
                       "Capriotti et al. statherin–HAP studies."),
    "metal_oxide": ("Carboxylate coordination to TiO₂/ZrO₂ surfaces. "
                    "Monti et al., J. Phys. Chem. C 2012; "
                    "Skelton et al., ACS Appl. Mater. Interfaces 2009 (titania-peptide)."),
    "metal_coord": ("Transition-metal (Zn/Ni/Cu/Fe) coordination by His, Cys, Asp, Glu. "
                    "ZincBind database, Nucleic Acids/Database 2019, baz006; "
                    "Alberts et al., J. Mol. Biol. 1998; Berg, J. Biol. Chem. 1990."),
    "salt_bridge": ("Salt bridges contribute ≈ 1–5 kcal/mol when solvent-exposed "
                    "(larger when buried). Kumar & Nussinov, "
                    "ChemBioChem 2002; Bosshard et al., J. Mol. Recognit. 2004."),
    "hbond": ("Surface H-bond ≈ 0.5–1.5 kcal/mol each. "
              "Fersht, Structure and Mechanism in Protein Science, 1999; "
              "Sheu et al., PNAS 2003."),
    "hydrophobic": ("Hydrophobic/CH–π contact ≈ 1–2 kcal/mol per nonpolar contact. "
                    "Eisenberg & McLachlan, Nature 1986, 319, 199 (atomic solvation parameters)."),
    "phosphate": ("Phosphate–Arg/Lys electrostatic & bidentate contacts. "
                  "Woods & Ferre, J. Proteome Res. 2005; Hendsch & Tidor, Protein Sci. 1994."),
}

# Map each functional group to the reference key(s) that justify its energies
FG_REFERENCES = {
    "-SH":        ["AuS"],
    "Ca2+":       ["Ca_carboxylate", "cation_pi"],
    "Ti-O":       ["metal_oxide"],
    "Zr-O":       ["metal_oxide"],
    "pi":         ["pi_pi", "cation_pi"],
    "-NH2":       ["salt_bridge", "hbond"],
    "-COOH":      ["salt_bridge", "hbond"],
    "-OH":        ["hbond"],
    "-CH3":       ["hydrophobic"],
    "PO4":        ["phosphate"],
    "-C=O":       ["hbond"],
    "Metal":      ["metal_coord"],
    "SaltBridge": ["salt_bridge"],
}

# Each entry: group -> {residue: (E_kcal_per_residue, mechanism_label)}
FG_INTERACTIONS: Dict[str, Dict[str, tuple]] = {
    "-SH": {                       # thiol / gold surface
        "CYS": (-45.0, "Au-S semi-covalent"),
        "MET": (-5.0,  "weak S coordination"),
    },
    "Ca2+": {                      # calcium sites (hydroxyapatite, calcium phosphate)
        "ASP": (-10.0, "carboxylate-Ca coordination"),
        "GLU": (-10.0, "carboxylate-Ca coordination"),
        "SER": (-6.0,  "hydroxyl/phospho-Ser binding"),
        "THR": (-4.0,  "hydroxyl binding"),
        "PHE": (-4.0,  "cation-π"),
        "TRP": (-6.0,  "cation-π"),
        "TYR": (-5.0,  "cation-π"),
    },
    "Ti-O": {                      # titanium / zirconium oxide
        "ASP": (-8.0,  "carboxylate-metal coordination"),
        "GLU": (-8.0,  "carboxylate-metal coordination"),
        "SER": (-3.0,  "hydroxyl-oxide H-bond"),
        "THR": (-3.0,  "hydroxyl-oxide H-bond"),
        "TYR": (-3.0,  "phenol-oxide H-bond"),
    },
    "Zr-O": {
        "ASP": (-8.0,  "carboxylate-metal coordination"),
        "GLU": (-8.0,  "carboxylate-metal coordination"),
        "SER": (-3.0,  "hydroxyl-oxide H-bond"),
        "THR": (-3.0,  "hydroxyl-oxide H-bond"),
    },
    "pi": {                        # aromatic / graphitic surface (π system)
        "PHE": (-4.0,  "π-π stacking"),
        "TRP": (-5.0,  "π-π stacking"),
        "TYR": (-4.0,  "π-π stacking"),
        "ARG": (-3.0,  "cation-π (guanidinium)"),
        "LYS": (-2.0,  "cation-π"),
        "HIS": (-2.5,  "π-π / cation-π"),
    },
    "-NH2": {                      # amine surface (cationic)
        "ASP": (-4.0,  "electrostatic + H-bond"),
        "GLU": (-4.0,  "electrostatic + H-bond"),
        "SER": (-0.8,  "H-bond acceptor"),
        "THR": (-0.8,  "H-bond acceptor"),
        "ASN": (-0.8,  "H-bond"),
        "GLN": (-0.8,  "H-bond"),
        "TYR": (-1.0,  "H-bond"),
    },
    "-COOH": {                     # carboxyl surface (anionic at neutral pH)
        "LYS": (-4.0,  "electrostatic + H-bond"),
        "ARG": (-4.0,  "electrostatic + H-bond"),
        "HIS": (-2.0,  "electrostatic (if protonated)"),
        "SER": (-0.8,  "H-bond donor"),
        "THR": (-0.8,  "H-bond donor"),
        "ASP": (+2.0,  "electrostatic repulsion"),
        "GLU": (+2.0,  "electrostatic repulsion"),
    },
    "-OH": {                       # hydroxyl surface (H-bonding, neutral)
        "SER": (-0.8,  "H-bond"),
        "THR": (-0.8,  "H-bond"),
        "TYR": (-1.0,  "H-bond"),
        "ASN": (-0.8,  "H-bond"),
        "GLN": (-0.8,  "H-bond"),
        "ASP": (-1.0,  "H-bond acceptor"),
        "GLU": (-1.0,  "H-bond acceptor"),
        "LYS": (-0.8,  "H-bond"),
        "ARG": (-0.8,  "H-bond"),
    },
    "-CH3": {                      # methyl / hydrophobic surface
        "ALA": (-1.0,  "hydrophobic contact"),
        "VAL": (-1.8,  "hydrophobic contact"),
        "LEU": (-2.2,  "hydrophobic contact"),
        "ILE": (-2.2,  "hydrophobic contact"),
        "MET": (-1.8,  "hydrophobic contact"),
        "PHE": (-2.0,  "hydrophobic + CH-π"),
        "TRP": (-2.2,  "hydrophobic + CH-π"),
        "PRO": (-1.0,  "hydrophobic contact"),
    },
    "PO4": {                       # phosphate surface
        "LYS": (-6.0,  "electrostatic (phosphate-ammonium)"),
        "ARG": (-7.0,  "bidentate (phosphate-guanidinium)"),
        "HIS": (-3.0,  "electrostatic"),
        "SER": (-1.0,  "H-bond"),
    },
    "-C=O": {                      # carbonyl surface (H-bond acceptor)
        "LYS": (-1.2,  "H-bond"),
        "ARG": (-1.2,  "H-bond"),
        "SER": (-1.0,  "H-bond donor"),
        "THR": (-1.0,  "H-bond donor"),
        "TRP": (-1.0,  "H-bond (indole NH)"),
    },
    "Metal": {                     # transition-metal ion surface (Zn/Ni/Cu/Fe)
        "HIS": (-12.0, "imidazole-metal coordination"),
        "CYS": (-10.0, "thiolate-metal coordination"),
        "ASP": (-6.0,  "carboxylate-metal coordination"),
        "GLU": (-6.0,  "carboxylate-metal coordination"),
        "MET": (-4.0,  "thioether-metal (soft)"),
    },
    "SaltBridge": {                # charged surface forming ionic pairs
        "LYS": (-5.0,  "ammonium-anion salt bridge"),
        "ARG": (-5.0,  "guanidinium-anion salt bridge"),
        "ASP": (-5.0,  "carboxylate-cation salt bridge"),
        "GLU": (-5.0,  "carboxylate-cation salt bridge"),
        "HIS": (-2.5,  "imidazolium salt bridge (if protonated)"),
    },
}

# Canonical list of mappable functional groups (for frontend display)
FG_LIST = list(FG_INTERACTIONS.keys())

# Reference SASA for normalising residue exposure (Å²)
SASA_REF = 50.0


# ════════════════════════════════════════════════════════════════════════════
# SURFACE GROUP DEFINITIONS  (Mode B: direct surface chemistry highlighting)
# ════════════════════════════════════════════════════════════════════════════
# Unlike the affinity engine (which scores how a patch binds an external
# surface), this maps the chemical groups the PROTEIN ITSELF presents on its
# solvent-accessible surface. Each entry lists residues bearing that group and
# the atoms that carry it (for optional atom-level display).
SURFACE_GROUP_DEFS = {
    "-COOH (carboxyl)":      {"ASP": ["OD1","OD2","CG"], "GLU": ["OE1","OE2","CD"]},
    "-OH (hydroxyl)":        {"SER": ["OG"], "THR": ["OG1"], "TYR": ["OH"]},
    "-SH (thiol)":           {"CYS": ["SG"]},
    "Aromatic ring":         {"PHE": ["CG","CD1","CD2","CE1","CE2","CZ"],
                              "TRP": ["CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2"],
                              "TYR": ["CG","CD1","CD2","CE1","CE2","CZ"],
                              "HIS": ["CG","ND1","CD2","CE1","NE2"]},
    "-NH2 / amine":          {"LYS": ["NZ"]},
    "Guanidinium":           {"ARG": ["NE","NH1","NH2","CZ"]},
    "Imidazole":             {"HIS": ["ND1","CE1","NE2","CD2","CG"]},
    "Metal-binding site":    {"HIS": ["ND1","NE2"], "CYS": ["SG"],
                              "ASP": ["OD1","OD2"], "GLU": ["OE1","OE2"]},
    "H-bond donor":          {"SER": ["OG"], "THR": ["OG1"], "TYR": ["OH"],
                              "ASN": ["ND2"], "GLN": ["NE2"], "TRP": ["NE1"],
                              "LYS": ["NZ"], "ARG": ["NE","NH1","NH2"], "HIS": ["ND1","NE2"]},
    "H-bond acceptor":       {"ASP": ["OD1","OD2"], "GLU": ["OE1","OE2"],
                              "ASN": ["OD1"], "GLN": ["OE1"], "SER": ["OG"],
                              "THR": ["OG1"], "TYR": ["OH"], "HIS": ["ND1","NE2"]},
    "Hydrophobic":           {"ALA": ["CB"], "VAL": ["CB","CG1","CG2"],
                              "LEU": ["CB","CG","CD1","CD2"], "ILE": ["CB","CG1","CG2","CD1"],
                              "MET": ["CB","CG","SD","CE"], "PHE": ["CG","CZ"],
                              "PRO": ["CB","CG","CD"]},
}
SURFACE_GROUP_LIST = list(SURFACE_GROUP_DEFS.keys())


# ════════════════════════════════════════════════════════════════════════════
# FEATURE AFFINITY ENGINE  (material-independent)
# ════════════════════════════════════════════════════════════════════════════
def functional_group_affinity(residues: list) -> dict:
    """
    For each functional group, compute how strongly THIS patch would bind a
    surface presenting that group. Sums per-residue interaction energies,
    weighted by solvent accessibility (only exposed residues can interact).

    Returns dict: group -> {dG, n_contacts, contributing_residues, mechanisms}
    """
    out = {}
    for group, residue_map in FG_INTERACTIONS.items():
        dG = 0.0
        contacts = []
        mechanisms = set()
        for r in residues:
            rn = r["res_name"]
            if rn in residue_map:
                E, mech = residue_map[rn]
                # weight by accessibility (buried residues can't reach surface)
                w = min(r.get("sasa", 0.0) / SASA_REF, 1.0)
                contribution = E * w
                dG += contribution
                if abs(contribution) > 0.05:
                    contacts.append({
                        "res": f"{rn}{r['res_seq']}{r['chain'] if r['chain']!='A' else ''}",
                        "dG": round(contribution, 2),
                        "mechanism": mech,
                    })
                    mechanisms.add(mech)
        # sort contributing residues by strength
        contacts.sort(key=lambda c: c["dG"])
        out[group] = {
            "dG": round(dG, 3),
            "n_contacts": len(contacts),
            "residues": contacts[:8],          # top 8 strongest
            "mechanisms": sorted(mechanisms),
        }
    return out


def surface_property_affinity(pm: dict, pb: Optional[dict],
                              env: EnvParams) -> dict:
    """
    Material-independent surface-property profile for this patch.
    Describes what KIND of surface this patch prefers, on physical axes
    that any material can be characterised by.

    Axes:
      charge_preference    : which surface charge sign is attractive
      hydrophobicity_pref  : hydrophilic vs hydrophobic surface
      hbond_capacity       : donor/acceptor potential
      electrostatic_strength : magnitude of electrostatic response (APBS)
      dielectric_sensitivity : how much it cares about ε (image charge)
    """
    nc = pm["net_charge"]

    # Charge preference: a negatively-charged patch prefers a positive surface
    if nc > 0.5:
        charge_pref = "negative surface (attracts +patch)"
        charge_sign = -1
    elif nc < -0.5:
        charge_pref = "positive surface (attracts −patch)"
        charge_sign = +1
    else:
        charge_pref = "neutral / weakly charged surface"
        charge_sign = 0

    # Hydrophobicity preference
    hydro = pm["hydro"]
    if hydro > 0.25:
        hydro_pref = "hydrophobic surface"
    elif hydro < -0.25:
        hydro_pref = "hydrophilic surface"
    else:
        hydro_pref = "amphipathic / mixed"

    # APBS electrostatic strength (if available)
    if pb:
        elec_strength = abs(pb.get("dG_kcal", 0.0))
        elec_method = "APBS_LinearPB"
    else:
        # Debye-Huckel proxy from net charge
        elec_strength = abs(nc) * 1.5
        elec_method = "DebyeHuckel_fallback"

    return {
        "net_charge":        round(nc, 2),
        "charge_preference": charge_pref,
        "charge_sign":       charge_sign,
        "hydrophobicity":    round(hydro, 3),
        "hydrophobicity_preference": hydro_pref,
        "hbond_donors":      pm["hb_donors"],
        "hbond_acceptors":   pm["hb_acceptors"],
        "aromatic_fraction": pm["ar_frac"],
        "pI":                pm["pI"],
        "dipole":            pm["dipole"],
        "electrostatic_strength": round(elec_strength, 3),
        "electrostatic_method":   elec_method,
        "isoelectric_note": (
            f"At pH {env.pH}, patch is "
            + ("positively" if nc > 0.5 else "negatively" if nc < -0.5 else "weakly")
            + " charged"
        ),
    }


def patch_intrinsic_penalties(pm: dict, env: EnvParams) -> dict:
    """
    Material-independent costs that always oppose adsorption, regardless of
    surface. These set a baseline that any attractive interaction must overcome.
      - orientational entropy loss (L11)
      - denaturation risk from secondary structure (L9, B-factor weighted L15)
    """
    dG_orient = dG_orientational(pm["n"], env.temp)

    # Denaturation: SS cost weighted by B-factor (material hydrophobicity removed)
    SS_COST = {"H": 2.0, "E": 1.5, "C": 0.3}
    B_REF = 30.0
    total = 0.0
    for r in pm["residues"]:
        B = r.get("bfac", B_REF)
        w = min(B / B_REF, 3.0)
        cost = SS_COST.get(r.get("ss", "C"), 0.3)
        total += w * cost
    dG_denat = (total / max(pm["n"], 1)) * 0.25

    return {
        "orientational_entropy": round(dG_orient, 4),
        "denaturation_risk":     round(dG_denat, 4),
        "total_penalty":         round(dG_orient + dG_denat, 4),
    }


def build_affinity_profile(pm: dict, pb: Optional[dict],
                           env: EnvParams) -> dict:
    """
    Complete material-independent affinity profile for a patch.
    Combines functional-group affinity + surface-property profile +
    intrinsic penalties, then produces a ranked recommendation list.
    """
    fg = functional_group_affinity(pm["residues"])
    sp = surface_property_affinity(pm, pb, env)
    pen = patch_intrinsic_penalties(pm, env)

    # Ranked functional-group recommendations (most attractive first)
    ranked_fg = sorted(
        [{"group": g, **v} for g, v in fg.items()],
        key=lambda x: x["dG"]
    )

    # Best functional group (net of intrinsic penalty)
    best = ranked_fg[0] if ranked_fg else None
    net_best = (best["dG"] + pen["total_penalty"]) if best else 0.0

    # Affinity score 0-100: how "useful" this patch is as a binding site
    # Based on the strongest available interaction overcoming penalties
    affinity_score = round(max(0.0, min(100.0, 50.0 - net_best * 2.0)), 1)

    return {
        "affinity_score":      affinity_score,
        "functional_groups":   fg,           # heat-map data: all groups
        "ranked_groups":       ranked_fg,    # ranked recommendation list
        "surface_properties":  sp,
        "penalties":           pen,
        "best_group":          best["group"] if best else None,
        "best_group_dG":       round(net_best, 3) if best else 0.0,
    }

# ════════════════════════════════════════════════════════════════════════════
# BIOPYTHON PIPELINE
# ════════════════════════════════════════════════════════════════════════════
def parse_struct(pdb: Path):
    return PDBParser(QUIET=True).get_structure("p", str(pdb))

def compute_sasa(struct) -> dict:
    sr = ShrakeRupley(probe_radius=1.40, n_points=100)
    sr.compute(struct, level="R")
    out = {}
    for model in struct:
        for chain in model:
            for res in chain:
                out[f"{chain.id}_{res.id[1]}"] = round(float(res.sasa), 2)
    return out

def compute_dssp(struct, pdb: Path) -> dict:
    """DSSP secondary structure. Backbone-geometry fallback if binary missing."""
    out = {}
    if MKDSSP:
        try:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                dssp = DSSP(struct[0], str(pdb), dssp=MKDSSP)
            for key in dssp:
                chain, res_id = key[0], key[1]
                code = dssp[key][2]
                ss   = "H" if code in ("H","G","I") else ("E" if code in ("B","E") else "C")
                out[f"{chain}_{res_id[1]}"] = ss
            return out
        except Exception as ex:
            log.warning(f"DSSP failed ({ex}), falling back to geometry")
    # Backbone geometry: d(CA_i-1, CA_i+1) < 6 → helix, > 7 → sheet
    for model in struct:
        for chain in model:
            cas = [r for r in chain if r.has_id("CA")]
            for i, res in enumerate(cas):
                key = f"{chain.id}_{res.id[1]}"
                if i == 0 or i == len(cas)-1:
                    out[key] = "C"; continue
                try:
                    d = float(cas[i-1]["CA"] - cas[i+1]["CA"])
                except: d = 10.0
                out[key] = "H" if d < 6.0 else ("E" if d > 7.0 else "C")
    return out

def build_residue_list(struct, sasa: dict, ss: dict,
                       pH: float, thr: float = 5.0):
    """Returns (surface_residues, atoms_by_res)."""
    surface = []; atoms_by_res = {}
    for model in struct:
        for chain in model:
            for res in chain:
                rn  = res.get_resname().strip()
                rid = res.id[1]; cid = chain.id
                key = f"{cid}_{rid}"
                if sasa.get(key, 0) < thr: continue
                if not res.has_id("CA"):   continue
                ca = res["CA"]
                surface.append({
                    "key": key, "res_name": rn, "res_seq": rid, "chain": cid,
                    "x": float(ca.coord[0]), "y": float(ca.coord[1]), "z": float(ca.coord[2]),
                    "charge_ph": residue_charge(rn, pH),
                    "sasa": sasa.get(key, 0.0),
                    "ss":   ss.get(key, "C"),
                    "bfac": round(float(ca.bfactor), 2),
                })
                atoms_by_res[key] = [
                    {"atomName": a.name,
                     "x": float(a.coord[0]), "y": float(a.coord[1]), "z": float(a.coord[2]),
                     "charge_pqr": 0.0}
                    for a in res.get_atoms() if a.element != "H"
                ]
    return surface, atoms_by_res

# ════════════════════════════════════════════════════════════════════════════
# PATCH CLUSTERING
# ════════════════════════════════════════════════════════════════════════════
def cluster(surface: list, pH: float, max_patches: int = 14,
            radius: float = 8.0, max_size: int = 30) -> list:
    """
    Build local surface patches.

    Previous version used a 12 Å BFS with no size cap, which on a large
    protein merged the entire surface into one giant connected component.
    This version grows compact, bounded patches:

      - each unvisited surface residue can seed a patch
      - a patch grows by adding nearby residues within `radius` (Å, Cα-Cα)
        of the SEED (not chained), keeping patches local and roughly disk-sized
      - a chemical-similarity bonus slightly extends the radius for like residues
      - each patch is capped at `max_size` residues
      - patches are ranked by size; the largest `max_patches` are returned
    """
    DIST = radius
    CHEM = 2.0
    visited = set()
    patches = []

    # seed order: most-buried-first tends to give better-centered patches,
    # but simple input order is fine and deterministic
    for seed in surface:
        if seed["key"] in visited:
            continue
        visited.add(seed["key"])
        patch = [seed]
        cat_seed = (AA.get(seed["res_name"]) or [None]*10)[9]

        # gather neighbours within radius OF THE SEED (local disk, not chained)
        candidates = []
        for other in surface:
            if other["key"] in visited:
                continue
            dx = seed["x"]-other["x"]; dy = seed["y"]-other["y"]; dz = seed["z"]-other["z"]
            d = math.sqrt(dx*dx + dy*dy + dz*dz)
            cat_o = (AA.get(other["res_name"]) or [None]*10)[9]
            thresh = DIST + (CHEM if cat_o == cat_seed else 0.0)
            if d < thresh:
                candidates.append((d, other))

        # nearest first, capped at max_size
        candidates.sort(key=lambda t: t[0])
        for _, other in candidates[:max_size-1]:
            visited.add(other["key"])
            patch.append(other)

        if len(patch) >= 3:        # ignore tiny 1-2 residue specks
            patches.append(patch)

    patches.sort(key=len, reverse=True)
    return patches[:max_patches]

# ════════════════════════════════════════════════════════════════════════════
# PATCH METRICS
# ════════════════════════════════════════════════════════════════════════════
def patch_metrics(patch: list, pH: float, atoms_by_res: dict,
                  grid: Optional[dict]) -> dict:
    n = len(patch)
    hy = nc = hbD = hbA = ar = 0.0
    ss_cnt = {"H": 0, "E": 0, "C": 0}
    residues = []; coords = []

    for r in patch:
        aa   = AA.get(r["res_name"])
        q    = r.get("charge_ph", residue_charge(r["res_name"], pH))
        nc  += q
        if aa:
            hy  += aa[0]; hbD += aa[3]; hbA += aa[4]
            if aa[9] == "aromatic": ar += 1
        ss = r.get("ss","C")
        ss_cnt[ss if ss in ss_cnt else "C"] += 1
        coords.append((r["x"], r["y"], r["z"]))
        residues.append({
            "res_name": r["res_name"], "res_seq": r["res_seq"],
            "chain": r["chain"], "key": r["key"],
            "charge": round(q, 3), "sasa": round(r.get("sasa",0.0), 2),
            "ss": ss, "category": aa[9] if aa else "other",
            "bfac": round(r.get("bfac", 0.0), 2),
        })

    cx = sum(c[0] for c in coords)/n
    cy = sum(c[1] for c in coords)/n
    cz = sum(c[2] for c in coords)/n
    maxR = max(math.sqrt((c[0]-cx)**2+(c[1]-cy)**2+(c[2]-cz)**2) for c in coords)
    area = math.pi * maxR * maxR

    # Dipole
    dpx = dpy = dpz = 0.0
    for r, c in zip(residues, coords):
        q = r["charge"]
        dpx += q*(c[0]-cx); dpy += q*(c[1]-cy); dpz += q*(c[2]-cz)
    dipole = math.sqrt(dpx**2 + dpy**2 + dpz**2)

    # PB electrostatics
    pb_elec = None
    if grid:
        T = 310.0  # will be overridden at scoring
        pb_elec = dG_PB(residues, atoms_by_res, grid, T)

    curv = curvature(coords)

    # Mean B-factor of patch (for L15)
    bfacs = [r["bfac"] for r in residues]
    mean_bfac = round(sum(bfacs)/len(bfacs), 2) if bfacs else 0.0
    max_bfac  = max(bfacs) if bfacs else 1.0

    # Isoelectric point of patch (for L12) — binary search
    pI_patch = _compute_pI([r["res_name"] for r in residues])

    return {
        "n": n,
        "hydro":       round(hy/n, 4),
        "net_charge":  round(nc, 3),
        "hb_donors":   int(hbD),
        "hb_acceptors":int(hbA),
        "ar_frac":     round(ar/n, 3),
        "area":        round(area, 1),
        "total_sasa":  round(sum(r["sasa"] for r in residues), 1),
        "dipole":      round(dipole, 3),
        "ss":          {k: round(v/n,3) for k,v in ss_cnt.items()},
        "curvature":   curv,
        "centroid":    {"x": round(cx,2), "y": round(cy,2), "z": round(cz,2)},
        "residues":    residues,
        "pb_elec":     pb_elec,
        "mean_bfac":   mean_bfac,
        "max_bfac":    round(max_bfac, 2),
        "pI":          round(pI_patch, 2),
    }

# ════════════════════════════════════════════════════════════════════════════
# ENDPOINTS
# ════════════════════════════════════════════════════════════════════════════
@app.get("/health")
def health():
    return {
        "status":      "ok",
        "version":     "4.0",
        "mode":        "feature-affinity-mapping",
        "apbs":        APBS    is not None,
        "pdb2pqr":     PDB2PQR is not None,
        "dssp":        MKDSSP  is not None,
        "functional_groups": FG_LIST,
        "surface_groups": SURFACE_GROUP_LIST,
        "surface_features": FEATURE_LIST,
        "compatibility_types": COMPATIBILITY_LIST,
        "references": REFERENCES,
        "fg_references": FG_REFERENCES,
    }


# ════════════════════════════════════════════════════════════════════════════
# RESIDUE-LEVEL SURFACE FEATURE MAPPING  (no patch clustering)
# ════════════════════════════════════════════════════════════════════════════
# These two endpoints give a direct, per-residue picture of (a) the chemically
# relevant features each solvent-accessible residue presents, and (b) how
# compatible each residue is with a chosen generic surface chemistry. No material
# is selected; nothing is clustered. Scoring is deliberately simple/transparent:
#
#     exposure_weight   = min(SASA / SASA_REF, 1.0)        # 0..1, how exposed
#     feature_intensity = exposure_weight                  # 0..1, Mode A
#     raw_dG            = exposure_weight * E(group,res)    # kcal/mol-like, Mode B

# Which residues present which residue-level feature (Mode A — Surface Features)
FEATURE_RESIDUES = {
    "charge_pos":     {"LYS", "ARG", "HIS"},
    "charge_neg":     {"ASP", "GLU"},
    "hbond_donor":    {"SER","THR","TYR","ASN","GLN","TRP","LYS","ARG","HIS"},
    "hbond_acceptor": {"ASP","GLU","ASN","GLN","SER","THR","TYR","HIS"},
    "hydrophobic":    {"ALA","VAL","LEU","ILE","MET","PHE","PRO","TRP"},
    "aromatic":       {"PHE","TRP","TYR","HIS"},
    "metal_binding":  {"HIS","CYS","ASP","GLU","MET"},
    "thiol":          {"CYS"},
    "carboxyl":       {"ASP","GLU"},
    "amine":          {"LYS","ARG","HIS"},
}

FEATURE_INFO = {
    "charge_pos":     ("Positive charge", "Cationic side chains (Lys, Arg, His+) — ionic & H-bond interactions."),
    "charge_neg":     ("Negative charge", "Anionic side chains (Asp, Glu) — ionic interactions, metal/Ca coordination."),
    "hbond_donor":    ("H-bond donor", "Donates hydrogen bonds (OH, NH groups)."),
    "hbond_acceptor": ("H-bond acceptor", "Accepts hydrogen bonds (carbonyl/hydroxyl O, ring N)."),
    "hydrophobic":    ("Hydrophobic", "Nonpolar side chains — hydrophobic & CH-π contacts."),
    "aromatic":       ("Aromatic", "Aromatic rings — π-π stacking and cation-π."),
    "metal_binding":  ("Metal-binding", "Coordinates transition metals (His, Cys, Asp, Glu, Met)."),
    "thiol":          ("Thiol (Cys)", "Free cysteine thiol — strong affinity for gold/soft metals."),
    "carboxyl":       ("Carboxyl (Asp/Glu)", "Carboxylate groups — Ca/oxide coordination, anionic."),
    "amine":          ("Amine (Lys/Arg/His)", "Basic side chains — cationic surface compatibility."),
}
FEATURE_LIST = list(FEATURE_RESIDUES.keys())

# Mode B — generic surface chemistry types mapped onto FG_INTERACTIONS.
# Each type names the FG_INTERACTIONS group it draws per-residue energies from.
COMPATIBILITY_TYPES = {
    "cationic":      ("-COOH",      "Cationic surface compatibility",
                      "Positively charged surface; binds anionic residues (Asp, Glu)."),
    "anionic":       ("-NH2",       "Anionic surface compatibility",
                      "Negatively charged surface; binds cationic residues (Lys, Arg, His)."),
    "hbond_donor":   ("-C=O",       "H-bond donor surface",
                      "Surface donates H-bonds; pairs with acceptor residues."),
    "hbond_acceptor":("-OH",        "H-bond acceptor surface",
                      "Surface accepts H-bonds; pairs with donor residues."),
    "pi_carbon":     ("pi",         "π / carbon-like surface",
                      "Graphitic/aromatic surface; π-π and cation-π with aromatic & cationic residues."),
    "hydrophobic":   ("-CH3",       "Hydrophobic surface",
                      "Nonpolar surface; hydrophobic contacts with aliphatic/aromatic residues."),
    "oxide":         ("Ti-O",       "Oxide surface compatibility",
                      "Metal-oxide (TiO2/ZrO2); carboxylate and hydroxyl coordination."),
    "hydroxyapatite":("Ca2+",       "Hydroxyapatite / Ca²⁺ compatibility",
                      "Calcium-rich surface; carboxylate-Ca coordination, cation-π."),
    "metal_coord":   ("Metal",      "Metal coordination compatibility",
                      "Transition-metal sites; His/Cys/Asp/Glu coordination."),
    "gold":          ("-SH",        "Gold affinity",
                      "Au surface; strong Cys thiol-gold bonding."),
    "phosphate":     ("PO4",        "Phosphate surface compatibility",
                      "Phosphate groups; electrostatic with Lys/Arg."),
}
COMPATIBILITY_LIST = list(COMPATIBILITY_TYPES.keys())


def _prepare_surface(req: "AnalyzeRequest", workdir: Path):
    """
    Shared setup with SCIENTIFIC DEPTH:
      - parse structure, SASA, DSSP, surface residue list (as before)
      - run pdb2pqr + APBS and attach the real electrostatic potential phi
        (kT/e) to every surface residue (graceful fallback if APBS missing)
      - compute a 3-D neighbourhood for each residue: the surface residues
        within NEIGHBOR_RADIUS Å, used to build a local-environment context so
        scoring reflects clustered chemistry rather than isolated residues.
    Returns (surface, stats) where stats carries electrostatics info.
    """
    NEIGHBOR_RADIUS = 8.0   # Å, Cα-Cα — local environment shell

    pdb = workdir / "input.pdb"
    if req.pdb_text:
        pdb.write_text(req.pdb_text)
    elif req.pdb_id:
        import urllib.request
        url = f"https://files.rcsb.org/download/{req.pdb_id.upper()}.pdb"
        urllib.request.urlretrieve(url, pdb)
    else:
        raise HTTPException(400, "Provide pdb_id or pdb_text")

    struct = parse_struct(pdb)
    sasa   = compute_sasa(struct)
    ss     = compute_dssp(struct, pdb)
    surface, _ = build_residue_list(struct, sasa, ss, req.env.pH)

    # ── electrostatics: pdb2pqr → APBS → potential at each residue Cα ──
    electrostatics = "none"
    grid = None
    if PDB2PQR:
        try:
            pqr = workdir / "mol.pqr"
            run_pdb2pqr(pdb, pqr, req.env.pH)
            if APBS:
                try:
                    stem = workdir / "pot"
                    apin = build_apbs_input(pqr, stem, req.env.ionic/1000.0, req.env.temp)
                    dx   = run_apbs(apin, workdir)
                    grid = parse_dx(dx)
                    electrostatics = "APBS_LinearPB"
                except Exception as e:
                    log.warning(f"APBS step failed, continuing without phi: {e}")
                    electrostatics = "pdb2pqr_only"
            else:
                electrostatics = "pdb2pqr_only"
        except Exception as e:
            log.warning(f"pdb2pqr failed: {e}")

    # attach phi (electrostatic potential, kT/e) to each surface residue
    for r in surface:
        if grid is not None:
            try:
                r["phi"] = round(phi_at(grid, r["x"], r["y"], r["z"]), 3)
            except Exception:
                r["phi"] = 0.0
        else:
            r["phi"] = None

    # ── 3-D neighbourhood: surface residues within NEIGHBOR_RADIUS of each ──
    # Two shells are built:
    #   _neighbors      : within NEIGHBOR_RADIUS (8 Å) — used for the local
    #                     'context' bonus in the propensity composite.
    #   _patch_neighbors: within patch_radius (user-set, default 12 Å) — used
    #                     for the SAP-style patch-density score (Chennamsetty
    #                     et al. 2009, generalised to surface chemistry).
    R_patch = float(getattr(req.env, "patch_radius", 12.0) or 12.0)
    n = len(surface)
    for i in range(n):
        ri = surface[i]
        neigh = []; patch = []
        for j in range(n):
            if i == j: continue
            rj = surface[j]
            dx = ri["x"]-rj["x"]; dy = ri["y"]-rj["y"]; dz = ri["z"]-rj["z"]
            d2 = dx*dx + dy*dy + dz*dz
            if d2 <= NEIGHBOR_RADIUS*NEIGHBOR_RADIUS:
                neigh.append(rj)
            if d2 <= R_patch*R_patch:
                patch.append(rj)
        ri["_neighbors"] = neigh          # internal use (not serialised raw)
        ri["n_neighbors"] = len(neigh)
        ri["_patch_neighbors"] = patch    # internal use for SAP-style density

    n_atoms = sum(1 for m in struct for c in m for r in c for _ in r.get_atoms())
    n_res   = sum(1 for m in struct for c in m for _ in c)
    stats = {"n_atoms": n_atoms, "n_residues": n_res,
             "n_surface_res": len(surface), "electrostatics": electrostatics,
             "patch_radius": R_patch}
    return surface, stats


@app.post("/surface_feature_map")
async def surface_feature_map(req: AnalyzeRequest):
    """
    Mode A — residue-level surface feature map (material-independent).
    For every solvent-accessible residue, report which chemical features it
    presents and a 0..1 intensity that combines solvent exposure with the
    feature. No patches, no material.
    """
    workdir = Path(tempfile.mkdtemp(prefix="ppia_feat_"))
    try:
        surface, stats = _prepare_surface(req, workdir)
        pH = req.env.pH

        residues = []
        for r in surface:
            rn = r["res_name"]
            exposure = min(r["sasa"] / SASA_REF, 1.0)
            q = residue_charge(rn, pH)
            labels = [f for f in FEATURE_LIST if rn in FEATURE_RESIDUES[f]]
            # charge labels are pH-aware: only flag if actually charged at this pH
            if "charge_pos" in labels and q <= 0.05: labels.remove("charge_pos")
            if "charge_neg" in labels and q >= -0.05: labels.remove("charge_neg")
            residues.append({
                "res_name": rn, "res_seq": r["res_seq"], "chain": r["chain"],
                "sasa": round(r["sasa"], 1),
                "exposure": round(exposure, 3),
                "charge": round(q, 2),
                "ss": r["ss"],
                "features": labels,
                # per-feature intensity = exposure (feature is binary-present here)
                "intensity": round(exposure, 3),
            })

        # per-feature summary: which residues, counts
        feature_summary = {}
        for f in FEATURE_LIST:
            members = [x for x in residues if f in x["features"]]
            members.sort(key=lambda m: m["intensity"], reverse=True)
            label, mech = FEATURE_INFO[f]
            feature_summary[f] = {
                "label": label, "mechanism": mech,
                "n": len(members),
                "mean_intensity": round(
                    sum(m["intensity"] for m in members)/len(members), 3) if members else 0.0,
                "residues": members,
            }

        return {
            "status": "ok",
            "mode": "surface-features",
            "stats": stats,
            "feature_list": FEATURE_LIST,
            "feature_info": {f: {"label": FEATURE_INFO[f][0],
                                 "mechanism": FEATURE_INFO[f][1]} for f in FEATURE_LIST},
            "features": feature_summary,
            "residues": residues,
            "note": ("Solvent-accessible residue-level surface feature mapping. "
                     "Intensity reflects solvent exposure, not a binding prediction."),
        }
    except HTTPException:
        raise
    except Exception as ex:
        log.error(f"surface_feature_map: {ex}", exc_info=True)
        raise HTTPException(500, str(ex))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


@app.post("/compatibility_map")
async def compatibility_map(req: AnalyzeRequest):
    """
    Mode B — residue-level functional-group compatibility map.
    The user picks a generic surface chemistry (cationic, gold, oxide, ...).
    For each solvent-accessible residue we look up its per-residue interaction
    energy in FG_INTERACTIONS for the corresponding group, weight it by solvent
    exposure, and return both the raw ΔG-like contribution and a 0..1 normalized
    compatibility score for colouring. No material selection, no clustering.
    """
    workdir = Path(tempfile.mkdtemp(prefix="ppia_compat_"))
    try:
        surface, stats = _prepare_surface(req, workdir)

        out = {}
        for ctype, (fg_key, label, mech) in COMPATIBILITY_TYPES.items():
            energies = FG_INTERACTIONS.get(fg_key, {})
            members = []
            raw_sum = 0.0
            for r in surface:
                rn = r["res_name"]
                if rn not in energies:
                    continue
                E, mechanism = energies[rn]
                exposure = min(r["sasa"] / SASA_REF, 1.0)
                raw = exposure * E          # kcal/mol-like contribution
                raw_sum += raw
                members.append({
                    "res_name": rn, "res_seq": r["res_seq"], "chain": r["chain"],
                    "sasa": round(r["sasa"], 1),
                    "exposure": round(exposure, 3),
                    "raw_dG": round(raw, 2),
                    "mechanism": mechanism,
                })
            # normalize: most attractive (most negative raw) -> 1.0
            if members:
                most_neg = min(m["raw_dG"] for m in members)
                most_pos = max(m["raw_dG"] for m in members)
                span = most_pos - most_neg
                for m in members:
                    if span < 1e-6:
                        # all equal (incl. single residue): attractive -> 1, else 0.5
                        m["score"] = 1.0 if m["raw_dG"] < 0 else 0.5
                    else:
                        # 1.0 at most attractive, 0.0 at least attractive/repulsive
                        m["score"] = round((most_pos - m["raw_dG"]) / span, 3)
                members.sort(key=lambda m: m["raw_dG"])   # strongest binders first
            out[ctype] = {
                "label": label, "mechanism": mech, "fg_group": fg_key,
                "n": len(members),
                "raw_sum": round(raw_sum, 1),
                "mean_score": round(sum(m["score"] for m in members)/len(members), 3) if members else 0.0,
                "residues": members,
            }

        return {
            "status": "ok",
            "mode": "compatibility",
            "stats": stats,
            "compatibility_list": COMPATIBILITY_LIST,
            "compatibility": out,
            "note": ("Generic functional-group compatibility mapping. Scores are "
                     "literature-energy-weighted estimates for ranking residues, "
                     "not predictions of exact adsorption sites."),
        }
    except HTTPException:
        raise
    except Exception as ex:
        log.error(f"compatibility_map: {ex}", exc_info=True)
        raise HTTPException(500, str(ex))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


@app.post("/analyze_surface")
async def analyze_surface(req: AnalyzeRequest):
    """
    Unified single-call analysis. Runs the shared pipeline ONCE and returns
    BOTH residue-level surface features (the protein's own chemistry) AND
    functional-group compatibility (fit to generic surface chemistries), so the
    frontend can offer every filter in one list without re-running anything.
    """
    workdir = Path(tempfile.mkdtemp(prefix="ppia_surf_"))
    try:
        # _prepare_surface runs blocking work (pdb2pqr, APBS). Offload to a
        # worker thread so the event loop stays responsive during computation.
        import anyio
        surface, stats = await anyio.to_thread.run_sync(_prepare_surface, req, workdir)
        pH = req.env.pH

        # ── per-residue feature labels (Mode A logic) ──
        residues = []
        for r in surface:
            rn = r["res_name"]
            exposure = min(r["sasa"] / SASA_REF, 1.0)
            q = residue_charge(rn, pH)
            labels = [f for f in FEATURE_LIST if rn in FEATURE_RESIDUES[f]]
            if "charge_pos" in labels and q <= 0.05: labels.remove("charge_pos")
            if "charge_neg" in labels and q >= -0.05: labels.remove("charge_neg")
            residues.append({
                "res_name": rn, "res_seq": r["res_seq"], "chain": r["chain"],
                "sasa": round(r["sasa"], 1), "exposure": round(exposure, 3),
                "charge": round(q, 2), "ss": r["ss"],
                "phi": r.get("phi"), "n_neighbors": r.get("n_neighbors", 0),
                "features": labels, "intensity": round(exposure, 3),
            })

        features = {}
        for f in FEATURE_LIST:
            members = [x for x in residues if f in x["features"]]
            members.sort(key=lambda m: m["intensity"], reverse=True)
            label, mech = FEATURE_INFO[f]
            features[f] = {
                "kind": "feature", "label": label, "mechanism": mech,
                "n": len(members),
                "mean_intensity": round(sum(m["intensity"] for m in members)/len(members), 3) if members else 0.0,
                "residues": members,
            }

        # ── per-residue compatibility scores (with 3D context + APBS) ──
        # Electrostatic surface types are weighted by the real APBS potential;
        # all types get a local-environment bonus when neighbouring surface
        # residues share the same compatible chemistry (clustered chemistry
        # binds better than an isolated residue).
        ELECTROSTATIC_TYPES = {"cationic": -1, "anionic": +1, "phosphate": +1}
        #   sign = which potential sign REINFORCES binding:
        #   cationic surface (negative) binds residues in positive potential? No —
        #   a cationic surface binds ANIONIC residues, which sit in negative phi.
        #   We use: bonus if residue phi has the sign that attracts the surface.
        def _phi_factor(ctype, phi):
            if phi is None or ctype not in ELECTROSTATIC_TYPES:
                return 1.0
            # cationic surface attracts residues sitting in NEGATIVE potential
            # anionic / phosphate surfaces attract residues in POSITIVE potential
            want_negative = (ctype == "cationic")
            aligned = (phi < 0) if want_negative else (phi > 0)
            mag = min(abs(phi) / 5.0, 1.0)          # saturate at |phi|=5 kT/e
            return 1.0 + (0.5 * mag if aligned else -0.3 * mag)

        compatibility = {}
        for ctype, (fg_key, label, mech) in COMPATIBILITY_TYPES.items():
            energies = FG_INTERACTIONS.get(fg_key, {})

            # ── SAP-style patch density (Chennamsetty et al. 2009, PNAS),
            # generalised from hydrophobicity to this surface chemistry.
            # For each surface residue i, sum over residues j within the patch
            # radius (i included) of: exposure(j) × |base_energy(j)| if j is
            # compatible with this chemistry, else 0. High values mark residues
            # sitting in the middle of a cluster of compatible, exposed residues.
            # The value is assigned to the central residue, exactly as in SAP.
            patch_raw = {}   # "chain:res_seq" -> raw patch density
            for r in surface:
                key = f"{r['chain']}:{r['res_seq']}"
                # central residue contribution + its patch neighbours
                shell = [r] + r.get("_patch_neighbors", [])
                s = 0.0
                for nb in shell:
                    if nb["res_name"] in energies:
                        E_nb, mech_nb = energies[nb["res_name"]]
                        exp_nb = min(nb["sasa"] / SASA_REF, 1.0)
                        # pH-dependent protonation scaling per neighbour, so the
                        # SAP-style patch density itself shifts with pH
                        prot_nb = _protonation_factor(nb["res_name"], mech_nb, pH)
                        s += exp_nb * abs(E_nb) * prot_nb
                patch_raw[key] = s
            # normalise patch density 0..100 across the protein for this chemistry
            pvals = [v for v in patch_raw.values() if v > 0]
            pmax = max(pvals) if pvals else 1.0
            pmax = pmax or 1.0

            members = []
            for r in surface:
                rn = r["res_name"]
                if rn not in energies: continue
                E, mechanism = energies[rn]      # literature energy (kcal/mol)
                exposure = min(r["sasa"] / SASA_REF, 1.0)

                # 3-D neighbourhood bonus: fraction of neighbours that also
                # carry this compatible chemistry (clustered chemistry wins)
                neigh = r.get("_neighbors", [])
                if neigh:
                    share = sum(1 for nb in neigh if nb["res_name"] in energies) / len(neigh)
                else:
                    share = 0.0
                context = 1.0 + 0.5 * share          # up to +50% for a full cluster

                # APBS electrostatic factor for charged surfaces
                phi_f = _phi_factor(ctype, r.get("phi"))

                # pH-dependent protonation factor (literature-based): scales the
                # residue's contribution by its protonation state for this
                # interaction mechanism (e.g. His+ cation-π gains at low pH;
                # Asp/Glu carboxylate-oxide coordination weakens at low pH).
                prot_f = _protonation_factor(rn, mechanism, pH)

                # composite (unitless) — NOT a free energy. This blends the
                # literature base energy with our exposure/context/phi weights,
                # so it is reported as a dimensionless "interaction propensity".
                composite = exposure * E * context * phi_f * prot_f
                key = f"{r['chain']}:{r['res_seq']}"
                patch_density = round(patch_raw.get(key, 0.0) / pmax * 100, 1)
                members.append({
                    "res_name": rn, "res_seq": r["res_seq"], "chain": r["chain"],
                    "sasa": round(r["sasa"], 1), "exposure": round(exposure, 3),
                    "phi": r.get("phi"),
                    "n_like_neighbors": int(round(share * len(neigh))) if neigh else 0,
                    "context": round(context, 2),
                    "patch_density": patch_density,   # SAP-style 0-100
                    "base_energy": E,            # literature value (kcal/mol)
                    "prot_factor": round(prot_f, 3),  # pH protonation scaling
                    "_composite": composite,     # internal, used for normalisation
                    "mechanism": mechanism,
                })
            if members:
                mn = min(m["_composite"] for m in members)
                mx = max(m["_composite"] for m in members)
                span = mx - mn
                for m in members:
                    # propensity: 0-100, 100 = strongest (most negative composite)
                    if span < 1e-6:
                        m["propensity"] = 100.0 if m["_composite"] < 0 else 50.0
                    else:
                        m["propensity"] = round((mx - m["_composite"]) / span * 100, 1)
                    del m["_composite"]
                members.sort(key=lambda m: m["propensity"], reverse=True)
            mean_prop = round(sum(m["propensity"] for m in members)/len(members), 1) if members else 0.0
            max_patch = max((m["patch_density"] for m in members), default=0.0)
            # the residue at the centre of the strongest patch
            top_patch_res = max(members, key=lambda m: m["patch_density"], default=None) if members else None
            compatibility[ctype] = {
                "kind": "compatibility", "label": label, "mechanism": mech,
                "fg_group": fg_key, "n": len(members),
                "mean_propensity": mean_prop,
                "max_patch_density": max_patch,
                "top_patch_residue": (f"{top_patch_res['res_name']}{top_patch_res['res_seq']} {top_patch_res['chain']}"
                                      if top_patch_res else None),
                "residues": members,
            }

        feat_list = list(FEATURE_LIST)

        return {
            "status": "ok",
            "mode": "surface",
            "stats": stats,
            "feature_list": feat_list,
            "compatibility_list": COMPATIBILITY_LIST,
            "feature_info": {f: {"label": FEATURE_INFO[f][0],
                                 "mechanism": FEATURE_INFO[f][1]} for f in FEATURE_LIST},
            "features": features,
            "compatibility": compatibility,
            "residues": residues,
            "note": ("Residue-level surface feature and generic functional-group "
                     "compatibility mapping. Scores rank residues on a consistent "
                     "scale; they are not predictions of exact adsorption sites."),
        }
    except HTTPException:
        raise
    except Exception as ex:
        log.error(f"analyze_surface: {ex}", exc_info=True)
        raise HTTPException(500, str(ex))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ════════════════════════════════════════════════════════════════════════════
# SERVE FRONTEND  (so the whole app is one URL: http://localhost:8000)
# ════════════════════════════════════════════════════════════════════════════
# The frontend lives in ../frontend/index.html relative to this file.
_FRONTEND_DIR = (Path(__file__).parent.parent / "frontend").resolve()

@app.get("/")
def _serve_index():
    idx = _FRONTEND_DIR / "index.html"
    if idx.exists():
        from fastapi.responses import HTMLResponse
        return HTMLResponse(idx.read_text(encoding="utf-8"))
    from fastapi.responses import PlainTextResponse
    return PlainTextResponse(
        "Frontend not found. Expected at: " + str(idx), status_code=404)

@app.get("/favicon.png")
@app.get("/favicon.ico")
def _serve_favicon():
    from fastapi.responses import FileResponse, Response
    for name in ("favicon.png", "logo.png"):
        p = _FRONTEND_DIR / name
        if p.exists():
            return FileResponse(str(p))
    return Response(status_code=404)

@app.get("/logo.png")
def _serve_logo():
    from fastapi.responses import FileResponse, Response
    p = _FRONTEND_DIR / "logo.png"
    return FileResponse(str(p)) if p.exists() else Response(status_code=404)

# also serve any other static assets in the frontend folder (if added later)
try:
    from fastapi.staticfiles import StaticFiles
    if _FRONTEND_DIR.exists():
        app.mount("/static", StaticFiles(directory=str(_FRONTEND_DIR)), name="static")
except Exception as _e:
    log.warning(f"static mount skipped: {_e}")


def _open_browser_later(url: str, delay: float = 1.5):
    """Open the default browser to the app once the server is up."""
    import threading, webbrowser, time
    def _go():
        time.sleep(delay)
        try: webbrowser.open(url)
        except Exception: pass
    threading.Thread(target=_go, daemon=True).start()


if __name__ == "__main__":
    import uvicorn
    PORT = int(os.environ.get("PPIA_PORT", "8000"))
    log.info(f"Open the app at http://localhost:{PORT}")
    _open_browser_later(f"http://localhost:{PORT}")
    uvicorn.run(app, host="0.0.0.0", port=PORT, reload=False)
