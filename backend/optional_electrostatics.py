"""Opt-in PDB2PQR/PROPKA/APBS workflow. Canonical 11 maps are unchanged.

This module requires locally installed PDB2PQR and (for the potential grid)
the APBS executable. It never synthesizes a successful result when either is
missing. It is protein-only; missing ligands and cofactors are disclosed.
"""
from __future__ import annotations
import io
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from dataclasses import asdict, dataclass
from pathlib import Path
import numpy as np
from prepare import prepare_pdb_text

class ElectrostaticsError(ValueError):
    pass

@dataclass(frozen=True)
class ElectrostaticSettings:
    pH: float = 7.4
    ionic_mM: float = 150.0
    temperature_K: float = 298.15
    electrolyte: str = "NaCl"
    forcefield: str = "PARSE"
    protein_dielectric: float = 2.0
    solvent_dielectric: float = 78.5
    apbs_enabled: bool = True
    propka_enabled: bool = True

    def validate(self):
        if not 0 <= self.pH <= 14: raise ElectrostaticsError("pH must be 0–14")
        if not 0 <= self.ionic_mM <= 1000: raise ElectrostaticsError("Ionic strength must be 0–1000 mM")
        if not 273.15 <= self.temperature_K <= 373.15: raise ElectrostaticsError("Temperature must be 273.15–373.15 K")
        if self.electrolyte not in ("NaCl","KCl"): raise ElectrostaticsError("Only NaCl/KCl supported")
        if self.forcefield not in ("PARSE","AMBER","CHARMM"): raise ElectrostaticsError("Unknown forcefield")
        if not 1 <= self.protein_dielectric <= 30: raise ElectrostaticsError("Protein dielectric must be 1–30")
        if not 20 <= self.solvent_dielectric <= 110: raise ElectrostaticsError("Solvent dielectric must be 20–110")

def executable():
    pqr=shutil.which("pdb2pqr30") or shutil.which("pdb2pqr")
    if pqr: return [pqr]
    try:
        import importlib.util
        if importlib.util.find_spec("pdb2pqr"): return [sys.executable,"-m","pdb2pqr"]
    except (ImportError,ValueError):
        pass
    return None

def status():
    return {"pdb2pqr":bool(executable()),"apbs":bool(shutil.which("apbs")),
            "optional":True,"changes_standard_maps":False}

def run(args,cwd,seconds):
    try:
        r=subprocess.run(args,cwd=cwd,capture_output=True,text=True,timeout=seconds,
                         env={**os.environ,"OPENBLAS_NUM_THREADS":"1"},check=False)
    except subprocess.TimeoutExpired as e:
        raise ElectrostaticsError("External computation timed out") from e
    except OSError as e:
        raise ElectrostaticsError("External executable failed to start: "+str(e)) from e
    log=(r.stdout or "")+"\n"+(r.stderr or "")
    if r.returncode: raise ElectrostaticsError("External calculation failed: "+log[-1500:])
    return log[-10000:]

def atoms_from_pqr(pqr):
    result=[]
    for line in pqr.splitlines():
        if not line.startswith(("ATOM  ","HETATM")): continue
        words=line.split()
        if len(words)<10: raise ElectrostaticsError("Incomplete PQR atom record")
        try:
            xyz=[float(x) for x in words[-5:-2]]
            charge,radius=map(float,words[-2:])
        except ValueError as e: raise ElectrostaticsError("Invalid PQR charge/radius") from e
        if not all(math.isfinite(x) for x in xyz+[charge,radius]) or radius<=0:
            raise ElectrostaticsError("Invalid or missing PQR parameters")
        result.append((xyz,charge,radius))
    if not result:raise ElectrostaticsError("No PQR atoms")
    return result

def input_for_apbs(settings,atoms):
    xyz=np.array([a[0] for a in atoms],dtype=float)
    lo=xyz.min(axis=0);hi=xyz.max(axis=0)
    extent=hi-lo
    fine=extent+24.0
    coarse=np.maximum(1.6*fine,extent+48.0)
    center=(lo+hi)/2
    choices=(65,97,129,161)
    dims=[next((d for d in choices if d>=math.ceil(n)+1),161) for n in fine]
    spacing=fine/(np.array(dims)-1)
    if max(spacing)>1.8: raise ElectrostaticsError("Structure exceeds safe local grid resolution")
    conc=settings.ionic_mM/1000
    ions=(f"ion charge +1 conc {conc:.6f} radius 2.0\n  ion charge -1 conc {conc:.6f} radius 2.0"
          if conc else "")
    inp=f"""read
  mol pqr prepared.pqr
end
ELEC name electrostatics
  mg-auto
  dime {dims[0]} {dims[1]} {dims[2]}
  cglen {coarse[0]:.4f} {coarse[1]:.4f} {coarse[2]:.4f}
  fglen {fine[0]:.4f} {fine[1]:.4f} {fine[2]:.4f}
  cgcent {center[0]:.4f} {center[1]:.4f} {center[2]:.4f}
  fgcent {center[0]:.4f} {center[1]:.4f} {center[2]:.4f}
  mol 1
  lpbe
  bcfl sdh
  pdie {settings.protein_dielectric:.4f}
  sdie {settings.solvent_dielectric:.4f}
  chgm spl2
  srfm smol
  srad 1.4
  swin 0.3
  sdens 10.0
  temp {settings.temperature_K:.3f}
  {ions}
  calcenergy no
  calcforce no
  write pot dx potential
end
quit
"""
    return inp,{"grid_dimensions":dims,"grid_spacing_A":spacing.tolist(),
                "potential_native_units":"kBT/e","electrolyte_concentration_M":conc,
                "solver":"linearized Poisson-Boltzmann"}

def compute(*,pdb_text,chain,settings):
    settings.validate()
    s=status()
    if not s["pdb2pqr"]:raise ElectrostaticsError("PDB2PQR is not installed")
    if settings.apbs_enabled and not s["apbs"]:raise ElectrostaticsError("APBS executable is not installed")
    if not pdb_text or len(pdb_text)>16000000:raise ElectrostaticsError("Invalid or oversized PDB")
    cleaned,prep=prepare_pdb_text(pdb_text,chain)
    lost=sorted({l[17:20].strip() for l in pdb_text.splitlines() if l.startswith("HETATM")})
    with tempfile.TemporaryDirectory(prefix="interfacescout_pb_") as td:
        root=Path(td)
        (root/"input.pdb").write_text(cleaned,encoding="utf-8")
        cmd=[*executable(),"--ff="+settings.forcefield,"--keep-chain"]
        if settings.propka_enabled:
            cmd.extend(["--titration-state-method=propka",f"--with-ph={settings.pH:.3f}"])
        cmd.extend(["input.pdb","prepared.pqr"])
        prep_log=run(cmd,root,150)
        out=root/"prepared.pqr"
        if not out.exists():raise ElectrostaticsError("PDB2PQR output missing")
        pqr=out.read_text(errors="replace")
        atoms=atoms_from_pqr(pqr)
        charge=sum(a[1] for a in atoms)
        warnings=[]
        if lost: warnings.append("Removed HETATM components (not modeled): "+", ".join(lost))
        if not settings.propka_enabled: warnings.append("PROPKA disabled: pH did not determine protonation")
        result={"status":"completed","settings":asdict(settings),"preparation":prep,
                "parameterized_atomic_charge_e":round(charge,5),
                "n_parameterized_atoms":len(atoms),
                "protonation":"PROPKA" if settings.propka_enabled else "PDB2PQR default",
                "warnings":warnings,"standard_maps_modified":False,
                "scope":"Protein-only continuum electrostatics; no adsorption energy or site prediction"}
        apbs_log=""
        if settings.apbs_enabled:
            inp,grid=input_for_apbs(settings,atoms)
            (root/"potential.in").write_text(inp)
            apbs_log=run([shutil.which("apbs"),"potential.in"],root,300)
            dx=root/"potential.dx"
            if not dx.exists() or dx.stat().st_size<80:raise ElectrostaticsError("APBS potential.dx missing/empty")
            result["apbs"]={"grid":grid,"data_file":"data/potential.dx",
                            "mV_per_kBT_e":round(settings.temperature_K*0.08617333262,6)}
        bundle=io.BytesIO()
        with zipfile.ZipFile(bundle,"w",zipfile.ZIP_DEFLATED,compresslevel=6) as z:
            z.writestr("metadata/results.json",json.dumps(result,indent=2))
            z.writestr("data/prepared.pqr",pqr)
            z.writestr("metadata/pdb2pqr.log",prep_log)
            z.writestr("metadata/prepared_protein.pdb",cleaned)
            if settings.apbs_enabled:
                z.writestr("metadata/apbs.in",inp)
                z.writestr("metadata/apbs.log",apbs_log)
                z.write(dx,"data/potential.dx")
        return result,bundle.getvalue()
