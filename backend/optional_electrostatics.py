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
from physics_surface import read_pqr, read_dx, charge_by_residue, sample_exposed_potentials

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

def _is_hydrogen_atom(name):
    import re
    return bool(re.match(r"^(?:[123])?H",str(name).upper()))

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
        atom_name=words[2]
        if not all(math.isfinite(x) for x in xyz+[charge,radius]) or radius<0 or (radius==0 and not _is_hydrogen_atom(atom_name)):
            raise ElectrostaticsError("Invalid or missing PQR parameters in line: "+repr(line[:180]))
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
        original_pqr=out.read_text(errors="replace")
        # Some PDB2PQR forcefields use charged zero-radius hydrogens.
        # These MUST remain in the APBS PQR: partial charges are physical inputs.
        # Only zero-charge/zero-radius hydrogens can be omitted without loss.
        filtered_lines=[]
        neutral_dummy_count=0
        for line in original_pqr.splitlines(keepends=True):
            if line.startswith(("ATOM  ","HETATM")):
                words=line.split()
                if len(words)>=10:
                    q=float(words[-2]);radius=float(words[-1])
                    if radius==0.0 and q==0.0 and _is_hydrogen_atom(words[2]):
                        neutral_dummy_count+=1
                        continue
                    if radius<0 or (radius==0 and not _is_hydrogen_atom(words[2])) or not all(math.isfinite(x) for x in (q,radius)):
                        raise ElectrostaticsError("Invalid atomic radius/charge: "+repr(line[:180]))
            filtered_lines.append(line)
        pqr="".join(filtered_lines)
        out.write_text(pqr,encoding="utf-8")
        atoms=atoms_from_pqr(pqr)
        charge=sum(a[1] for a in atoms)
        detailed_pqr=read_pqr(pqr)
        residue_charge=charge_by_residue(detailed_pqr)
        warnings=[]
        if neutral_dummy_count:
            warnings.append(f"Removed {neutral_dummy_count} zero-charge/zero-radius nonphysical PDB2PQR placeholders from the APBS input; original PQR preserved")
        if lost: warnings.append("Removed HETATM components (not modeled): "+", ".join(lost))
        if not settings.propka_enabled: warnings.append("PROPKA disabled: pH did not determine protonation")
        result={"status":"completed","settings":asdict(settings),"preparation":prep,
                "parameterized_atomic_charge_e":round(charge,5),
                "n_parameterized_atoms":len(atoms),
                "protonation":"PROPKA" if settings.propka_enabled else "PDB2PQR default",
                "warnings":warnings,"standard_maps_modified":False,
                "scope":"Protein-only continuum electrostatics; no adsorption energy or site prediction"}
        apbs_log=""
        per_residue=[]
        per_patch=[]
        sampled={}
        # Canonical maps are read-only references; APBS data are independent descriptors.
        from interface_engine import analyze_all_maps
        standard=analyze_all_maps(pdb_text=cleaned,pH=settings.pH,
             ionic_mM=settings.ionic_mM,temp_K=settings.temperature_K)
        surface_keys=set(standard["diagnostics"]["surface_residue_keys"])
        matched=surface_keys & set(residue_charge)
        if len(matched)<max(15,int(0.7*len(surface_keys))):
            raise ElectrostaticsError("Protein/PQR residue mapping insufficient: "
               +str(len(matched))+"/"+str(len(surface_keys)))
        if settings.apbs_enabled:
            inp,grid=input_for_apbs(settings,atoms)
            (root/"potential.in").write_text(inp)
            apbs_log=run([shutil.which("apbs"),"potential.in"],root,300)
            dx=root/"potential.dx"
            if not dx.is_file():
                # APBS builds differ in whether they suffix the single
                # processor ID; accept only an unambiguous potential grid.
                found=sorted(p for p in root.glob("potential*.dx") if p.is_file())
                if len(found)==1:
                    dx=found[0]
            if not dx.exists() or dx.stat().st_size<80:
                files=", ".join(sorted(p.name for p in root.iterdir()))
                raise ElectrostaticsError("APBS potential DX missing/empty; output files: "+
                    files+"; APBS log: "+apbs_log[-2000:])
            sampled=sample_exposed_potentials(detailed_pqr,read_dx(dx.read_text(errors="replace")),
                         settings.temperature_K,surface_keys)
            if len(sampled)<max(10,int(.25*len(surface_keys))):
                raise ElectrostaticsError("Insufficient finite solvent-side APBS samples: "
                  +str(len(sampled))+"/"+str(len(surface_keys)))
            result["apbs"]={"grid":grid,"data_file":"data/potential.dx",
                            "mV_per_kBT_e":round(settings.temperature_K*0.08617333262,6)}
        physical_view=None
        if settings.apbs_enabled:
            from physics_aware_11maps import build_physical_11
            physical=build_physical_11(
                prepared_protein_pdb=cleaned,pqr_text=pqr,dx_text=dx.read_text(errors="replace"),
                pH=settings.pH,temperature_K=settings.temperature_K,
                baseline_analysis=standard)
            # Preserve the canonical API schema for rendering, but clearly mark
            # this as an experimentally unvalidated electrostatics-aware variant.
            physical_view=dict(standard)
            physical_view["maps"]=physical["maps"]
            physical_view["version"]=physical["model"]
            physical_view["physical_context"]=physical["physics"]
            physical_view["scope"]={**standard["scope"],
                "electrostatics_aware_experimental_variant":True}
            physical_view["method"]={**standard["method"],
                "electrostatics_aware_variant":True,
                "potential_gate":"fraction of solvent-side points of compatible potential sign for explicit ionic mechanisms",
                "ionization_source":"PDB2PQR/PROPKA assigned-state; not a titration ensemble"}
            result["physics_aware_11map_model"]={
                "version":physical["model"],"n_maps":len(physical["maps"]),
                "mapped_surface_residues":physical["physics"]["charged_residues_mapped"],
                "potential_sampled_residues":physical["physics"]["potential_sampled_residues"],
                "candidate_regions_are_not_adsorption_predictions":True}
        import csv
        for profile in standard["residue_profiles"]:
            key=profile["key"]
            per_residue.append({"residue_key":key,"residue_name":profile["res_name"],
               "surface_exposed":bool(profile["surface_exposed"]),
               "assigned_PQR_residue_charge_e":residue_charge.get(key),
               **sampled.get(key,{})})
        for key in standard["map_order"]:
            for p in standard["maps"][key]["patches"]:
                vals=[sampled[k]["potential_mean_mV"] for k in p["members"] if k in sampled]
                per_patch.append({"map_key":key,"patch_rank":p["display_rank"],
                     "center_key":p["center_key"],"member_count":len(p["members"]),
                     "potential_covered_members":len(vals),
                     "mean_member_potential_mV":float(np.mean(vals)) if vals else None,
                     "potential_sign":"positive" if vals and np.mean(vals)>0 else ("negative" if vals else "unavailable"),
                     "unchanged_standard_patch_rank":True})
        result["regional_surface_coverage"]={"surface_residue_count":len(surface_keys),
            "PQR_mapped_count":len(matched),"APBS_potential_sampled_count":len(sampled),
            "method":"PQR heavy sidechain atom probe-expanded solvent-side points; steric masking; equal angular samples; not exact SES"}
        def to_csv(rows):
            fields=sorted({k for r in rows for k in r});out=io.StringIO()
            writer=csv.DictWriter(out,fieldnames=fields);writer.writeheader();writer.writerows(rows)
            return out.getvalue()
        bundle=io.BytesIO()
        with zipfile.ZipFile(bundle,"w",zipfile.ZIP_DEFLATED,compresslevel=6) as z:
            z.writestr("metadata/results.json",json.dumps(result,indent=2))
            if physical_view is not None:
                z.writestr("data/physics_aware_11maps.json",json.dumps(physical_view))
                z.writestr("data/canonical_11maps.json",json.dumps(standard))
            z.writestr("data/prepared.pqr",pqr)
            z.writestr("metadata/pdb2pqr_unfiltered.pqr",original_pqr)
            z.writestr("data/charged_residue_surface.csv",to_csv(per_residue))
            z.writestr("data/standard_region_APBS_overlay.csv",to_csv(per_patch))
            z.writestr("metadata/pdb2pqr.log",prep_log)
            z.writestr("metadata/prepared_protein.pdb",cleaned)
            if settings.apbs_enabled:
                z.writestr("metadata/apbs.in",inp)
                z.writestr("metadata/apbs.log",apbs_log)
                z.write(dx,"data/potential.dx")
        result["canonical_view"]=standard
        if physical_view is not None:
            result["physics_aware_view"]=physical_view
        return result,bundle.getvalue()
