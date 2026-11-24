"""Focused, reproducible structural/physics QA for ISV76 exceptions.

No in-place substitution of structures, fabricated atoms or parameter fitting.
- 1GCT: inspect ALL incomplete backbone residues, PQR failure, optional explicit
  deletion-only diagnostic of severely incomplete residues, and 1YPH alternate.
- 8C50: actual chain charges, PQR geometry and LPBE 0.8/0.6/1.1 Å
  numerical consistency plus an independently solved nonlinear PB field.
All original maps and the official 76-member list remain unchanged.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import sys
import tempfile
import traceback
from collections import defaultdict

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"backend"),str(ROOT/"validation")]
from isv_full_80_qa import prepare, CANON
from electrostatics import _apbs_input, _pqr_atoms, _run
from surface_potential import sample_residue_sas_potential
from interface_engine import analyze_all_maps
from tv48_numerical_convergence import load_dx_streaming
from tv48_apbs_sas_grid_audit import run_apbs as run_grid
OUT=Path("isv76_special_cases")
OUT.mkdir(exist_ok=True)
ACTIVE=list(csv.DictReader((ROOT/"validation"/"ISV_active_candidates_76_v3.csv").open(encoding="utf-8")))
BY={r["proposed_pdb_id"]:r for r in ACTIVE}


def save(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")


def source_characterization(pid, row):
    meta,pdbfile=prepare(row)
    pdb=pdbfile.read_text()
    by=defaultdict(lambda:{"atoms":set(),"n":0})
    for line in pdb.splitlines():
        if line.startswith("ATOM "):
            key=(line[21],int(line[22:26]),line[26].strip(),line[17:20].strip())
            by[key]["atoms"].add(line[12:16].strip())
            by[key]["n"]+=1
    issues=[]
    for (chain,n,icode,res),v in by.items():
        if res not in CANON:continue
        missing=sorted({"N","CA","C","O"}-v["atoms"])
        if missing:
            issues.append({"chain":chain,"resid":n,"icode":icode,"resname":res,
                           "observed_atoms":sorted(v["atoms"]),"missing_backbone_atoms":missing,
                           "n_heavy_atoms":v["n"],"removal_candidate":len(missing)>=2})
    return meta,pdbfile,pdb,issues


def pqr_prepare(wd,pdb):
    (wd/"structure.pdb").write_text(pdb)
    text=_run([shutil.which("pdb2pqr"),"--ff=AMBER","--keep-chain",
               "--titration-state-method=propka","--whitespace",
               "--with-ph=7.4000","structure.pdb","structure.pqr"],wd,timeout=660)
    return _pqr_atoms(wd/"structure.pqr"),text


def gct_case():
    row=BY["1GCT"]
    meta,pdbfile,pdb,issues=source_characterization("1GCT",row)
    data={"pdb":"1GCT","source_protein_annotation_in_RCSB":"gamma-chymotrypsin",
          "active_manifest_name":"Alpha chymotrypsin (scientifically imprecise; not automatically changed here)",
          "original_bioassembly_sha256":meta["assembly_cif_sha256"],
          "original_prepared_pdb_sha256":meta["prepared_pdb_sha256"],
          "n_chains":meta["n_chains"],
          "n_incomplete_backbone_residues":len(issues),"issues":issues,
          "policy":"No atom reconstruction; any deletion is a separately labeled input-perturbation diagnostic."}
    with tempfile.TemporaryDirectory() as td:
        wd=Path(td)
        try:
            atoms,_=pqr_prepare(wd,pdb)
            data["original_pqr"]={"status":"PASS","n_atoms":len(atoms)}
        except Exception as exc:
            data["original_pqr"]={"status":"FAIL","error":str(exc)[-2000:]}
        to_delete={(x["chain"],x["resid"],x["icode"]) for x in issues
                   if x["removal_candidate"]}
        if not to_delete:
            data["deletion_diagnostic"]={"status":"NOT_APPLICABLE"}
        else:
            def match(line):
                return line.startswith("ATOM ") and (
                    (line[21],int(line[22:26]),line[26].strip()) in to_delete
                )
            edited="\n".join(l for l in pdb.splitlines() if not match(l))+"\n"
            editcount=sum(match(l) for l in pdb.splitlines())
            try:
                atoms,log=pqr_prepare(wd,edited)
                new=analyze_all_maps(pdb_text=edited,pH=7.4,ionic_mM=150.,temp_K=298.)
                old=analyze_all_maps(pdb_text=pdb,pH=7.4,ionic_mM=150.,temp_K=298.)
                original_keys={x["key"] for x in old["residue_profiles"]}
                edited_keys={x["key"] for x in new["residue_profiles"]}
                if original_keys-edited_keys != {f"{c}:{n}:{ic}" for c,n,ic in to_delete}:
                    raise AssertionError("Removal altered unanticipated identity")
                grid, m=run_grid(wd,atoms,150.,.8,"trimmed")
                sas,sm=sample_residue_sas_potential(
                    atoms,*grid,{r["key"] for r in new["residue_profiles"] if r["surface_exposed"]}
                )
                del grid
                oldcount=sum(len(v.get("patches",[])) for v in old["maps"].values())
                newcount=sum(len(v.get("patches",[])) for v in new["maps"].values())
                data["deletion_diagnostic"]={
                    "status":"PASS_DIAGNOSTIC_NOT_OFFICIAL",
                    "removed_residue_keys":sorted(f"{c}:{n}:{ic}" for c,n,ic in to_delete),
                    "removed_source_atom_records":editcount,
                    "original_residues":len(old["residue_profiles"]),
                    "remaining_residues":len(new["residue_profiles"]),
                    "old_11map_patches":oldcount,"new_11map_patches":newcount,
                    "n_pqr_atoms_after_deletion":len(atoms),
                    "n_SAS_potentials":len(sas),
                    "apbs_grid":m,
                    "synthetic_atoms_added":0,
                    "interpretation":"diagnostic structural truncation; not validation of original 1GCT",
                }
            except Exception as exc:
                data["deletion_diagnostic"]={
                    "status":"FAIL_DIAGNOSTIC","error":str(exc)[-1700:],
                    "traceback":traceback.format_exc()[-1600:]}
        # 1YPH alpha-chymotrypsin alternate explicitly separate accession
        altrow={**row,"proposed_pdb_id":"1YPH"}
        try:
            altmeta,altpdb,altstr,altissues=source_characterization("1YPH",altrow)
            ata,log=pqr_prepare(wd,altstr)
            alc=analyze_all_maps(pdb_text=altstr,pH=7.4,ionic_mM=150.,temp_K=298.)
            grid, gmeta=run_grid(wd,ata,150.,.8,"alternate")
            sas,smeta=sample_residue_sas_potential(
                ata,*grid,{r["key"] for r in alc["residue_profiles"] if r["surface_exposed"]}
            )
            data["alternate_1YPH"]={
                "status":"PASS_ALTERNATE_NOT_AUTOMATIC_SUBSTITUTION",
                "source_assembly_sha256":altmeta["assembly_cif_sha256"],
                "prepared_pdb_sha256":altmeta["prepared_pdb_sha256"],
                "n_chains":altmeta["n_chains"],"n_residues":len(alc["residue_profiles"]),
                "n_backbone_incomplete":len(altissues),
                "n_11_maps":alc["n_maps"],
                "n_exposed_sas_potentials":len(sas),"grid":gmeta,
                "scope":"1YPH is high-resolution bovine alpha-chymotrypsin, whereas 1GCT is gamma-chymotrypsin; biologically distinct entries.",
            }
        except Exception as exc:
            data["alternate_1YPH"]={
                "status":"FAILED_ALTERNATE","error":str(exc)[-2000:],
                "traceback":traceback.format_exc()[-1300:]}
    save(OUT/"1GCT_coordinates_and_1YPH_alternate.json",data)
    print("ISV76_GCT_QA "+json.dumps({k:v.get("status") for k,v in data.items()
          if isinstance(v,dict) and "status" in v}),flush=True)
    return data


def score_compare(a,b):
    common=sorted(set(a)&set(b))
    if len(common)<30:raise ValueError("Too few matched SAS residues")
    x=np.asarray([a[k]["surface_potential_kT_per_e"] for k in common],float)
    y=np.asarray([b[k]["surface_potential_kT_per_e"] for k in common],float)
    return {
        "n_common":len(common),
        "median_abs_delta_kT_per_e":round(float(np.median(np.abs(x-y))),7),
        "p95_abs_delta_kT_per_e":round(float(np.quantile(np.abs(x-y),.95)),7),
        "max_abs_delta_kT_per_e":round(float(np.max(np.abs(x-y))),7),
        "spearman":round(float(spearmanr(x,y).statistic),7),
        "median_relative_abs_delta_to_reference":round(float(np.median(np.abs(x-y))/(np.median(np.abs(y))+1e-10)),7)
    }


def make_custom_apbs(wd,atoms,label,kind,ionic=150.):
    source=_apbs_input(atoms,ionic,298.)
    if kind=="npbe":
        source=source.replace("    lpbe\n","    npbe\n")
    source=source.replace("write pot dx potential",f"write pot dx {label}")
    fn=wd/(label+".in")
    fn.write_text(source)
    log=_run([shutil.which("apbs"),fn.name],wd,timeout=1500)
    candidates=[wd/(label+".dx")] if (wd/(label+".dx")).exists() else list(wd.glob(label+"-*.dx"))
    if len(candidates)!=1:raise ValueError("APBS yielded missing or ambiguous DX file")
    grid=load_dx_streaming(candidates[0])
    candidates[0].unlink(missing_ok=True)
    return grid,log


def fibril_case():
    row=BY["8C50"]
    meta,pdbfile,pdb,issues=source_characterization("8C50",row)
    analysis=analyze_all_maps(pdb_text=pdb,pH=7.4,ionic_mM=150.,temp_K=298.)
    exposed={r["key"] for r in analysis["residue_profiles"] if r["surface_exposed"]}
    result={
        "pdb":"8C50","source_sha256":meta["assembly_cif_sha256"],
        "n_chains":meta["n_chains"],"n_residues":len(analysis["residue_profiles"]),
        "original_structure_missing_backbone":len(issues),
        "experimental_scope":"7.6 A medium-resolution cryoEM / modeling-assisted CsgA fibril; not atomically determined sidechain conformations",
        "structure_class":"biofilm_fibril","reference":{"pH":7.4,"ionic_mM":150.,
            "model":"PDB2PQR/PROPKA/linear-PB","SAS_directions":42},
        "variants":{},"comparison_to_reference":{}}
    with tempfile.TemporaryDirectory() as td:
        wd=Path(td)
        atoms,plog=pqr_prepare(wd,pdb)
        charges=defaultdict(float);sizes=defaultdict(int)
        for atom in atoms:
            charges[atom["chain"]]+=float(atom["charge"])
            sizes[atom["chain"]]+=1
        result["charges_by_chain_e"]={k:round(v,6) for k,v in charges.items()}
        result["pqr_atoms_by_chain"]=dict(sizes)
        result["net_pqr_charge_e"]=round(sum(charges.values()),6)
        xyz=np.asarray([a["xyz"] for a in atoms if a["radius"]>=.7],float)
        tree=cKDTree(xyz)
        nearest=tree.query(xyz,k=2)[0][:,1]
        result["heavy_atom_very_close_pairs_lt_0p5A"]=int(sum(nearest<.5)/2)
        variants=[
            ("lpbe_standard_08",.8,"lpbe",150.),
            ("lpbe_coarse_11",1.1,"lpbe",150.),
            ("lpbe_finer_06",.6,"lpbe",150.),
            ("lpbe_300mM",.8,"lpbe",300.),
            ("npbe_standard_08",.8,"npbe",150.),
        ]
        sas_by={}
        for label,target,kind,ionic in variants:
            try:
                if target!=.8:
                    grid,md=run_grid(wd,atoms,ionic,target,label)
                    log=""
                else:
                    grid,log=make_custom_apbs(wd,atoms,label,kind,ionic)
                    md={"grid_dims":list(grid[0].shape),
                        "effective_spacing_A":list(map(float,grid[2]))}
                potentials,smeta=sample_residue_sas_potential(
                    atoms,*grid,exposed,n_directions=42
                )
                del grid
                x=np.asarray([v["surface_potential_kT_per_e"]
                         for v in potentials.values()],float)
                perchain=defaultdict(list)
                for k,v in potentials.items():
                    perchain[k.split(":")[0]].append(v["surface_potential_kT_per_e"])
                result["variants"][label]={
                    "status":"PASS","kind":kind,"ionic_mM":ionic,
                    "grid_target_A":target,"grid":md,
                    "n_SAS_residues":len(potentials),
                    "mean_kT_per_e":round(float(x.mean()),5),
                    "median_kT_per_e":round(float(np.median(x)),5),
                    "minimum_kT_per_e":round(float(x.min()),5),
                    "per_chain_median_kT_per_e":{
                        c:round(float(np.median(val)),5) for c,val in perchain.items()
                    },
                    "apbs_solver_warning_lines":[v.strip()[:350] for v in log.splitlines()
                        if re.search("warn|failure|converge|error",v,re.I)][-15:],
                }
                sas_by[label]=potentials
            except Exception as exc:
                result["variants"][label]={
                    "status":"FAIL","error":str(exc)[-1800:],
                    "traceback":traceback.format_exc()[-1000:]}
            print("ISV76_FIBRIL_QA "+json.dumps({
                "label":label, "status":result["variants"][label]["status"],
                "median":result["variants"][label].get("median_kT_per_e")
            }),flush=True)
        ref=sas_by.get("lpbe_standard_08")
        if ref:
            for label,other in sas_by.items():
                if label!="lpbe_standard_08":
                    result["comparison_to_reference"][label]=score_compare(other,ref)
        result["physical_interpretation"]="Extreme potential magnitudes require quarantine irrespective of solver return code; compare model assumptions and numerical stability."
    save(OUT/"8C50_real_PB_variants_and_chain_QA.json",result)
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--case",choices=("1GCT","8C50"),required=True)
    a=parser.parse_args()
    result=gct_case() if a.case=="1GCT" else fibril_case()
    print("ISV76_SPECIAL_QA_SUMMARY "+json.dumps({
        "pdb":a.case,
        "completed":True,"file":str(OUT),
        "overall":"diagnostic evidence, not automatic biological acceptance"
    }),flush=True)


if __name__=="__main__":main()
