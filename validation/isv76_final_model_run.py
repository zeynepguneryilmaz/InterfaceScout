"""ISV76: frozen TV48-calibrated InterfaceScout applied to new, active structures.

Provenance: ISV_active_candidates_76_v3.csv is the ONLY live inclusion list.
Never score the historical 80 as 76; preserve structural/cofactor caveats.

Use the same production analyze_all_maps engine and its 11 maps, real
PDB2PQR/PROPKA/APBS (pH 7.4, 150 mM, 298 K, 42 SAS directions) and the
standard exportable residue/patch fields. APBS failures are retained as
explicit failures without misrepresenting maps that actually succeeded.
Results are descriptive technical surface atlas, not 76 experimentally
verified protein/material adsorption sites.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import traceback

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"validation"),str(ROOT/"backend")]
from isv_full_80_qa import prepare, CANON
from interface_engine import analyze_all_maps
from export_bundle import MAP_NAMES
from model_settings import MODEL_VERSION

ACTIVE=ROOT/"validation/ISV_active_candidates_76_v3.csv"
PANEL=list(csv.DictReader(ACTIVE.open(encoding="utf-8",newline="")))
assert len(PANEL)==76 and len({q["proposed_pdb_id"] for q in PANEL})==76
EXCLUDED={"3GHG","1QJS","4DN0","6LQH"}
assert not ({q["proposed_pdb_id"] for q in PANEL}&EXCLUDED)
assert "1STP" in {q["proposed_pdb_id"] for q in PANEL}
OUT=Path("isv76_final_model")
OUT.mkdir(exist_ok=True)
EXPECTED_MAPS=list(MAP_NAMES)
PH=7.4
IONIC=150.
TEMP=298.
SHARDS=16


def write_json(path,object):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(object,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")


def write_csv(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    if not rows:
        raise ValueError("No rows to export")
    names=list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w",encoding="utf-8",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=names,extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def top5_members(m):
    return {str(key) for region in m.get("patches",[])[:5]
            for key in region.get("members",[])}


def main_one(candidate):
    pid=candidate["proposed_pdb_id"]
    start=time.monotonic()
    structure,pdbpath=prepare(candidate)
    text=pdbpath.read_text(encoding="utf-8")
    source_sha=structure["assembly_cif_sha256"]
    prepared_sha=hashlib.sha256(text.encode()).hexdigest()
    if prepared_sha!=structure["prepared_pdb_sha256"]:
        raise AssertionError(pid+": prepared coordinates SHA mismatch")
    fixed=dict(pdb_text=text,pH=PH,ionic_mM=IONIC,temp_K=TEMP,chain=None)
    base=analyze_all_maps(**fixed)
    if base.get("version")!=MODEL_VERSION or list(base.get("map_order",[]))!=EXPECTED_MAPS:
        raise AssertionError(pid+": changed version/map order from TV48")
    if len(base.get("maps",{}))!=11 or len(base["residue_profiles"])<1:
        raise AssertionError(pid+": incomplete 11 map result")
    profiles=base["residue_profiles"]
    keys={r["key"] for r in profiles}
    if len(keys)!=len(profiles):
        raise AssertionError(pid+": duplicate residue identifiers")
    if any(r["res_name"]=="GLY" and
           (float(r["scrsa"])!=0 or float(r["sidechain_sasa"])!=0)
           for r in profiles):
        raise AssertionError(pid+": glycine sidechain artifact")
    all_map_members={}
    patches=[]
    mapsummary=[]
    for name in EXPECTED_MAPS:
        m=base["maps"][name]
        for patch in m.get("patches",[]):
            patches.append({"pdb_id":pid,"map":name,**{
                k:v for k,v in patch.items()
                if k not in ("center_residue","member_residues","seed_residues")
                and not isinstance(v,(dict,list,tuple))
            }, "members":";".join(map(str,patch.get("members",[])))})
        mapsummary.append({
            "pdb_id":pid, "map":name,"n_patches":len(m.get("patches",[])),
            "n_favorable_residues":len(m.get("residues",[])),
            "n_top5_residues":len(top5_members(m)),
            "n_pareto_front_1":sum(int(p.get("pareto_front",0))==1
                                    for p in m.get("patches",[])),
        })
        all_map_members[name]=top5_members(m)
    # Canonical map / surface data always saved, even if the optional real
    # external solver is not compatible with a candidate's input.
    rows=[]
    for r in profiles:
        key=r["key"]
        rows.append({
            "pdb_id":pid,"residue_key":key,"res_name":r["res_name"],
            "chain":r["chain"],"res_seq":r["res_seq"],"icode":r.get("icode",""),
            "surface_exposed":bool(r["surface_exposed"]),
            "trsa_raw":r["trsa_raw"],"scrsa_raw":r["scrsa_raw"],
            "residue_sasa_A2":r.get("residue_sasa"),
            "sidechain_sasa_A2":r.get("sidechain_sasa"),
        })
    for name in EXPECTED_MAPS:
        m=base["maps"][name]
        favorable={str(q["key"]):q for q in m.get("residues",[])}
        top5=all_map_members[name]
        for rec in rows:
            key=rec["residue_key"]
            v=favorable.get(key)
            rec[name+"_propensity"]=float(v.get("propensity") or 0.) if v else 0.
            rec[name+"_local_score"]=float(v.get("local_score") or 0.) if v else 0.
            rec[name+"_in_top5"]=key in top5

    protein=OUT/"proteins"/pid
    write_csv(protein/"canonical_residue_11map_values.csv",rows)
    write_csv(protein/"canonical_11map_patch_summary.csv",mapsummary)
    if patches:write_csv(protein/"canonical_patches_all_maps.csv",patches)
    write_json(protein/"canonical_maps_complete.json",base)

    status="computed"
    error=None
    el=None
    n_sas=0
    # Run the actual production physics API, NOT the historical TV48 audit
    # implementation. The same protein / same 11 maps must remain identical.
    try:
        enriched=analyze_all_maps(**fixed,electrostatics=True)
        if (base["maps"] != enriched["maps"] or
            base["residue_profiles"] != enriched["residue_profiles"] or
            base["map_order"] != enriched["map_order"]):
            raise AssertionError("APBS changed canonical map/residue computations")
        el=enriched["electrostatics"]
        if el.get("status")!="computed" or el.get("n_matched_residues")!=len(profiles):
            raise AssertionError("Missing real APBS or incomplete residue matching")
        n_sas=sum(bool(q["surface_exposed"]) and q.get("surface_potential_kT_per_e") is not None
                  for q in el["residues"])
        if n_sas!=el["n_surface_potential_residues"]:
            raise AssertionError("SAS potential count inconsistent")
        if n_sas<.9*sum(q["surface_exposed"] for q in profiles):
            raise AssertionError("SAS coverage below 90%")
        write_json(protein/"real_APBS_SAS_electrostatics.json",el)
        write_csv(protein/"real_APBS_SAS_per_residue.csv",[
            {"pdb_id":pid,**q} for q in el["residues"]
        ])
        status="computed"
    except Exception as exc:
        status="failed"
        error={"type":type(exc).__name__,"message":str(exc),
               "traceback_tail":traceback.format_exc()[-4000:]}
        write_json(protein/"APBS_ERROR.json",error)

    meta={
        "candidate_id":candidate["candidate_id"],
        "pdb_id":pid,"protein":candidate["protein_name"],
        "category":candidate["category"],
        "status_11maps":"PASS","status_real_APBS_SAS":status,
        "APBS_error":error,
        "n_maps":11,"n_residues":len(profiles),
        "n_exposed":sum(r["surface_exposed"] for r in profiles),
        "n_surface_APBS_residues":n_sas,
        "n_map_instances_with_patches":sum(r["n_patches"]>0 for r in mapsummary),
        "n_total_patches":sum(r["n_patches"] for r in mapsummary),
        "n_top5_unique_across_11":len(set().union(*all_map_members.values())),
        "n_protein_chains":structure["n_chains"],
        "missing_sidechain_residues":structure["sidechain_incomplete_residues"],
        "missing_backbone_residues":structure["backbone_incomplete_residues"],
        "prepared_source_sha256":source_sha,
        "prepared_PDB_sha256":prepared_sha,
        "model_version":base["version"],
        "pH":PH,"ionic_mM":IONIC,"temperature_K":TEMP,
        "SAS_directions_per_heavy_atom":42,
        "n_original_maps_unchanged_on_APBS":(11 if status=="computed" else None),
        "scientific_acceptance":"CANDIDATE_REVIEW_ONLY",
        "candidate_structural_caveat":candidate.get("structural_caveat",""),
        "ISV76_membership_sha256":hashlib.sha256(ACTIVE.read_bytes()).hexdigest(),
        "note":"No experimental contact-site validation or readjusted map parameters.",
        "elapsed_sec":round(time.monotonic()-start,2),
    }
    write_json(protein/"SUMMARY.json",meta)
    print("ISV76_FINAL_CASE "+json.dumps({
        "pdb_id":pid,"maps":11,"n_residues":len(profiles),"sas":n_sas,
        "APBS":status,"error":(error or {}).get("message","")[:200],
        "elapsed_sec":meta["elapsed_sec"]}),flush=True)
    return meta


def shard(index, shards):
    items=[row for i,row in enumerate(PANEL) if i%shards==index]
    results=[]
    for row in items:
        pid=row["proposed_pdb_id"]
        try:
            result=main_one(row)
        except Exception as e:
            result={
                "pdb_id":pid,"protein":row["protein_name"],"category":row["category"],
                "status_11maps":"FAIL","status_real_APBS_SAS":"not_run",
                "error_type":type(e).__name__,"error":str(e),
                "traceback_tail":traceback.format_exc()[-5000:],
            }
            write_json(OUT/"failures"/f"{pid}.json",result)
        results.append(result)
    write_json(OUT/f"shard_{index}.json",results)
    failmaps=sum(x["status_11maps"]!="PASS" for x in results)
    failapbs=sum(x["status_real_APBS_SAS"]!="computed" for x in results)
    print("ISV76_FINAL_SHARD "+json.dumps({
        "shard":index,"attempted":len(items),"n_11maps_pass":len(items)-failmaps,
        "n_APBS_pass":len(items)-failapbs,"failed_11maps":failmaps,"failed_APBS":failapbs
    }),flush=True)
    # CI is intentionally strict: all 76 full physics runs must be accounted
    # for, and noncompleted calculations are not presented as successes.
    if failmaps or failapbs:
        raise SystemExit(f"ISV76 shard {index}: maps {failmaps} failed; real APBS {failapbs} failed")


def aggregate(shards):
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    records=[]
    missing=[]
    for i in range(shards):
        p=OUT/f"shard_{i}.json"
        if p.is_file():records.extend(json.loads(p.read_text(encoding="utf-8")))
        else:missing.append(i)
    selected={row["proposed_pdb_id"] for row in PANEL}
    returned={r["pdb_id"] for r in records}
    if len(records)!=76 or len(returned)!=76 or returned!=selected or missing:
        raise ValueError("ISV76 missing/double counted structures: "+str(missing))
    records.sort(key=lambda r:next(i for i,v in enumerate(PANEL)
                                    if v["proposed_pdb_id"]==r["pdb_id"]))
    mapsok=[r for r in records if r["status_11maps"]=="PASS"]
    apsok=[r for r in records if r["status_real_APBS_SAS"]=="computed"]
    write_csv(OUT/"ISV76_ALL_STRUCTURE_STATUS.csv",records)
    maprows=[]
    electro_rows=[]
    residue_long=[]
    for p in mapsok:
        root=OUT/"proteins"/p["pdb_id"]
        maps=pd.read_csv(root/"canonical_11map_patch_summary.csv",dtype={"pdb_id":"string"})
        maprows.append(maps)
        df=pd.read_csv(root/"canonical_residue_11map_values.csv",dtype={"pdb_id":"string","residue_key":"string"})
        if df["residue_key"].duplicated().any():
            raise ValueError("Duplicate 11-map residue keys")
        if len(df)!=int(p["n_residues"]):
            raise ValueError("Map residue counts inconsistent")
        residue_long.append(df)
    for p in apsok:
        root=OUT/"proteins"/p["pdb_id"]
        el=pd.read_csv(root/"real_APBS_SAS_per_residue.csv",dtype={"pdb_id":"string","key":"string"})
        if len(el)!=int(p["n_residues"]):
            raise ValueError("Real APBS map/residue mismatch")
        electro_rows.append(el)
    # PDB identifiers are immutable strings: e.g. 1E05 must NOT be
    # parsed as floating-point scientific notation (100000.0).
    for tables in (maprows,residue_long,electro_rows):
        for df in tables:
            if df["pdb_id"].isna().any():
                raise ValueError("Missing PDB ID in per-protein data")
            if any(pid not in selected for pid in df["pdb_id"]):
                raise ValueError("PDB identifier corrupted by CSV parsing")
    if maprows: pd.concat(maprows,ignore_index=True).to_csv(OUT/"ISV76_ALL_11_MAP_METRICS.csv",index=False)
    if residue_long:pd.concat(residue_long,ignore_index=True).to_csv(OUT/"ISV76_ALL_RESIDUE_11_MAPS.csv",index=False)
    if electro_rows:pd.concat(electro_rows,ignore_index=True).to_csv(OUT/"ISV76_ALL_REAL_SAS_APBS_RESIDUES.csv",index=False)
    success=bool(len(mapsok)==76 and len(apsok)==76)
    # Finished aggregate audit *always* preserves failures, and fails CI
    # unless every candidate is computed in both the frozen map model and APBS.
    summary={
        "status":"PASS" if success else "INCOMPLETE",
        "active_manifest":"ISV_active_candidates_76_v3.csv",
        "manifest_sha256":hashlib.sha256(ACTIVE.read_bytes()).hexdigest(),
        "n_expected":76,"n_map_pass":len(mapsok),"n_real_APBS_pass":len(apsok),
        "n_maps_complete":sum(int(x.get("n_maps",0)) for x in mapsok),
        "n_patches":sum(int(x.get("n_total_patches",0)) for x in mapsok),
        "n_exposed":sum(int(x.get("n_exposed",0)) for x in mapsok),
        "n_sas_potentials":sum(int(x.get("n_surface_APBS_residues",0)) for x in apsok),
        "n_missing_sidechain_residues":sum(int(x.get("missing_sidechain_residues",0)) for x in mapsok),
        "n_backbone_incomplete_residues":sum(int(x.get("missing_backbone_residues",0)) for x in mapsok),
        "n_structures_with_any_map_without_patches":
          sum(int(x.get("n_map_instances_with_patches",0))<11 for x in mapsok),
        "n_biofilm_candidates":sum(q["category"].startswith("biofilm_") for q in PANEL),
        "n_scientifically_accepted_experimental_cases":0,
        "unresolved_biological_scope":["1STP biotin-removed representation","1HC1 very large assembly",
                                       "8C50 curli-like fibril; nonglobular scope",
                                       "missing sidechains/cofactors/glycans in protein-only models"],
        "failed_cases":[{
            "pdb_id":q["pdb_id"],"map_status":q["status_11maps"],
            "APBS_status":q["status_real_APBS_SAS"],
            "reason":(q.get("APBS_error") or {}).get("message",q.get("error",""))
        } for q in records if q["status_11maps"]!="PASS" or q["status_real_APBS_SAS"]!="computed"],
        "sensitivity_parameters_fixed":{"pH":PH,"ionic_mM":IONIC,"temperature_K":TEMP,
             "SAS_directions":42,"APBS_target_grid_A":.8},
        "interpretation":"Technical fixed-parameter application to ISV76 candidate set, not experimental adsorption validation."
    }
    write_json(OUT/"SUMMARY.json",summary)
    print("ISV76_FINAL_AGGREGATE "+json.dumps(summary),flush=True)

    if maprows:
        m=pd.concat(maprows,ignore_index=True)
        fig, ax=plt.subplots(figsize=(8.8,4.1))
        data=[m.loc[m["map"]==k,"n_patches"].values for k in EXPECTED_MAPS]
        ax.boxplot(data,tick_labels=EXPECTED_MAPS,showfliers=False)
        ax.tick_params(axis="x",rotation=53)
        ax.set_ylabel("Patches per structure")
        fig.tight_layout()
        root=OUT/"figures";root.mkdir(exist_ok=True)
        fig.savefig(root/"01_patches_by_map.png",dpi=250,bbox_inches="tight")
        fig.savefig(root/"01_patches_by_map.pdf",bbox_inches="tight")
        plt.close(fig)
        df=pd.DataFrame(mapsok)
        fig,ax=plt.subplots(figsize=(7.5,4.2))
        ax.scatter(df["n_residues"],df["n_total_patches"],s=16)
        ax.set_xlabel("Prepared structure residues")
        ax.set_ylabel("Total candidate patches across 11 maps")
        fig.tight_layout()
        fig.savefig(root/"02_structure_size_vs_patches.png",dpi=250,bbox_inches="tight")
        fig.savefig(root/"02_structure_size_vs_patches.pdf",bbox_inches="tight")
        plt.close(fig)
    if electro_rows:
        ae=pd.DataFrame(apsok)
        fig,ax=plt.subplots(figsize=(7.1,4.2))
        ax.scatter(ae["n_exposed"],ae["n_surface_APBS_residues"],s=14)
        bounds=max(float(ae.n_exposed.max()),float(ae.n_surface_APBS_residues.max()))
        ax.plot([0,bounds],[0,bounds],lw=.8,ls="--")
        ax.set_xlabel("tRSA-exposed residues")
        ax.set_ylabel("Residues with SAS-APBS potential")
        fig.tight_layout()
        fig.savefig(OUT/"figures"/"03_SAS_coverage.png",dpi=250,bbox_inches="tight")
        plt.close(fig)
    if not success:raise SystemExit("ISV76 incomplete: consult SUMMARY.json and per-protein failures")


def cli():
    p=argparse.ArgumentParser()
    p.add_argument("--shard",type=int)
    p.add_argument("--shards",type=int,default=SHARDS)
    p.add_argument("--aggregate",action="store_true")
    a=p.parse_args()
    if a.aggregate:aggregate(a.shards)
    elif a.shard is not None and 0<=a.shard<a.shards:shard(a.shard,a.shards)
    else:p.error("Choose --shard n or --aggregate")
if __name__=="__main__":
    cli()
