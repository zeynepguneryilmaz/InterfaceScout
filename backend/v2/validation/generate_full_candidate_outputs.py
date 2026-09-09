"""Generate full, user-facing InterfaceScout outputs for a blinded candidate panel.

This is a batch exporter, not a new scoring model. It uses the publication-frozen
V2 chemistry corrections/settings, saves the complete residue analysis, exports
all chemistry-map views to B-factor PDB files, and adds transparent structural
network/dynamics tables for reproducibility.

Experimental ground-truth labels are not read by this script.
"""
from __future__ import annotations

import csv
import json
import math
import os
import urllib.request
from pathlib import Path
from typing import Any, Dict, Iterable

import numpy as np

from v2.chemistry_freeze import apply_publication_chemistry
from v2.interface_engine import analyze_interface_v2
from v2.prepare import prepare_pdb_text
from v2.gnm import solve_gnm
from v2.rin import build_rin, annotate_rin_percentiles
from v2.model_settings import MULTISCALE_RADII_A, COARSE_PATCH_RADIUS_A

HERE = Path(__file__).resolve().parent
DEFAULT_MANIFEST = HERE / "candidate_prediction_panel_2026.json"
OUTROOT = Path("full_candidate_outputs")
AA = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}


def _fetch_pdb(pdb_id: str) -> str:
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=60) as r:
        return r.read().decode("utf-8", errors="replace")


def _write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows and not fieldnames:
        path.write_text("", encoding="utf-8")
        return
    if fieldnames is None:
        seen=[]
        for row in rows:
            for k in row:
                if k not in seen: seen.append(k)
        fieldnames=seen
    with path.open("w", newline="", encoding="utf-8-sig") as fh:
        w=csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


def _jsonable(x: Any) -> Any:
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, (np.integer,)): return int(x)
    if isinstance(x, (np.floating,)): return float(x)
    if isinstance(x, set): return sorted(x)
    if isinstance(x, dict): return {str(k): _jsonable(v) for k,v in x.items() if k not in {"atoms","center"}}
    if isinstance(x, (list,tuple)): return [_jsonable(v) for v in x]
    return x


def _res_key(r: dict) -> str:
    return str(r.get("key") or f"{r.get('chain','')}:{r.get('res_seq','')}:{r.get('icode','')}")


def _map_by_key(rows: Iterable[dict], keyfield: str="key") -> dict[str,dict]:
    return {str(r.get(keyfield)): r for r in rows if r.get(keyfield) is not None}


def _bfactor_pdb(prepared_pdb: str, values: dict[str,float], remarks: list[str]) -> str:
    out=[f"REMARK 900 {x}" for x in remarks]
    last_chain=None
    for line0 in prepared_pdb.splitlines():
        rec=line0[:6].strip()
        if rec != "ATOM":
            continue
        resn=line0[17:20].strip().upper()
        if resn not in AA: continue
        chain=line0[21] if len(line0)>21 else " "
        try: resseq=int(line0[22:26].strip())
        except Exception: continue
        icode=(line0[26] if len(line0)>26 else " ").strip()
        key=f"{chain}:{resseq}:{icode}"
        v=max(0.0,min(100.0,float(values.get(key,0.0) or 0.0)))
        if last_chain is not None and chain != last_chain: out.append("TER")
        last_chain=chain
        line=line0.ljust(66)
        out.append(line[:60]+f"{v:6.2f}"+line[66:])
    if last_chain is not None: out.append("TER")
    out.append("END")
    return "\n".join(out)+"\n"


def _normalize_0_100(values: dict[str,float]) -> dict[str,float]:
    arr=[abs(float(v)) for v in values.values() if v is not None and math.isfinite(float(v))]
    mx=max(arr) if arr else 0.0
    return {k:(100.0*float(v)/mx if mx>0 else 0.0) for k,v in values.items()}


def _feature_rows(v1res: dict) -> list[dict]:
    allr=v1res.get("all_residues") or v1res.get("surface_residues") or []
    sets={f:set(str(x) for x in (v1res.get("features",{}).get(f,{}) or {}).get("residues",[])) for f in v1res.get("feature_list",[])}
    rows=[]
    for r in allr:
        key=_res_key(r); row={"key":key,"res_name":r.get("res_name"),"res_seq":r.get("res_seq"),"icode":r.get("icode"),"chain":r.get("chain")}
        for f,s in sets.items(): row[f]=1 if key in s else 0
        rows.append(row)
    return rows


def _comprehensive_rows(v1res: dict) -> tuple[list[dict], list[str]]:
    allr=v1res.get("all_residues") or v1res.get("surface_residues") or []
    feature_list=list(v1res.get("feature_list",[]))
    chemistry_list=list(v1res.get("chemistry_list",[]))
    fsets={f:set(str(x) for x in (v1res.get("features",{}).get(f,{}) or {}).get("residues",[])) for f in feature_list}
    inner,outer=MULTISCALE_RADII_A; itag=f"{int(inner)}A"; otag=f"{int(outer)}A"; ctag=f"{int(COARSE_PATCH_RADIUS_A)}A"
    headers=[
        "analysis_version","selected_chain","pH","ionic_strength_mM","temperature_K",
        "residue_key","res_name","res_seq","icode","chain",
        "ca_x_A","ca_y_A","ca_z_A","sidechain_centroid_x_A","sidechain_centroid_y_A","sidechain_centroid_z_A",
        "total_sasa_A2","sidechain_sasa_A2","scrsa_raw","scrsa_clipped","surface_exposed",
        "charge_fraction","charge_descriptor","reference_pKa","ionization_sensitive",
        "secondary_structure","original_bfactor","apbs_phi_aux_kT_per_e","neighbors_8A_aux"
    ]
    headers += [f"feature_{f}" for f in feature_list]
    for c in chemistry_list:
        headers += [
            f"{c}_favorable_member",f"{c}_state_requirement",f"{c}_state_availability",f"{c}_auxiliary_charged_fraction",
            f"{c}_local_score",f"{c}_propensity",f"{c}_patch_{itag}_raw",f"{c}_patch_{otag}_raw",
            f"{c}_patch_{itag}_norm",f"{c}_patch_{otag}_norm",f"{c}_multiscale_persistence",f"{c}_multiscale_geomean",
            f"{c}_compatible_members_{ctag}",f"{c}_electrostatic_relation",f"{c}_ebase_metadata_kcal_mol",f"{c}_mechanism",
            f"{c}_repulsive_member",f"{c}_repulsion_state_requirement",f"{c}_repulsion_state_availability",
            f"{c}_repulsion_score",f"{c}_repulsion_propensity",f"{c}_repulsion_ebase_metadata_kcal_mol",f"{c}_repulsion_mechanism"
        ]
    chem_maps={}
    for c in chemistry_list:
        g=v1res.get("chemistries",{}).get(c,{})
        chem_maps[c]={"fav":_map_by_key(g.get("residues",[])),"rep":_map_by_key(g.get("repulsive_residues",[])),"patch":_map_by_key(g.get("patch_centers",[]),"center_key")}
    settings=v1res.get("settings",{})
    rows=[]
    for r in allr:
        k=_res_key(r)
        row={
            "analysis_version":v1res.get("version"),"selected_chain":settings.get("chain"),"pH":settings.get("pH"),
            "ionic_strength_mM":settings.get("ionic_mM"),"temperature_K":settings.get("temperature_K"),
            "residue_key":k,"res_name":r.get("res_name"),"res_seq":r.get("res_seq"),"icode":r.get("icode",""),"chain":r.get("chain"),
            "ca_x_A":r.get("x"),"ca_y_A":r.get("y"),"ca_z_A":r.get("z"),"sidechain_centroid_x_A":r.get("sc_x"),
            "sidechain_centroid_y_A":r.get("sc_y"),"sidechain_centroid_z_A":r.get("sc_z"),"total_sasa_A2":r.get("total_sasa"),
            "sidechain_sasa_A2":r.get("sidechain_sasa"),"scrsa_raw":r.get("scrsa_raw"),"scrsa_clipped":r.get("scrsa"),
            "surface_exposed":r.get("surface_exposed"),"charge_fraction":r.get("charge_fraction"),"charge_descriptor":r.get("charge_descriptor"),
            "reference_pKa":r.get("pka"),"ionization_sensitive":r.get("ionization_sensitive"),"secondary_structure":r.get("ss"),
            "original_bfactor":r.get("bfactor"),"apbs_phi_aux_kT_per_e":r.get("phi"),"neighbors_8A_aux":r.get("n_neighbors_8A",0)
        }
        for f in feature_list: row[f"feature_{f}"]=1 if k in fsets[f] else 0
        for c in chemistry_list:
            maps=chem_maps[c]; fav=maps["fav"].get(k,{}); rep=maps["rep"].get(k,{}); pat=maps["patch"].get(k,{})
            row.update({
                f"{c}_favorable_member":1 if k in maps["fav"] else 0,f"{c}_state_requirement":fav.get("state_requirement"),
                f"{c}_state_availability":fav.get("state_availability"),f"{c}_auxiliary_charged_fraction":fav.get("auxiliary_charged_fraction"),
                f"{c}_local_score":fav.get("local_score"),f"{c}_propensity":fav.get("propensity"),
                f"{c}_patch_{itag}_raw":pat.get(f"density_{itag}_raw"),f"{c}_patch_{otag}_raw":pat.get(f"density_{otag}_raw"),
                f"{c}_patch_{itag}_norm":pat.get(f"density_{itag}_norm"),f"{c}_patch_{otag}_norm":pat.get(f"density_{otag}_norm"),
                f"{c}_multiscale_persistence":pat.get("multiscale_persistence"),f"{c}_multiscale_geomean":pat.get("multiscale_geomean"),
                f"{c}_compatible_members_{ctag}":";".join(str(x) for x in pat.get(f"compatible_members_{ctag}",[])),
                f"{c}_electrostatic_relation":fav.get("electrostatic_relation"),f"{c}_ebase_metadata_kcal_mol":fav.get("ebase_metadata_kcal_mol"),
                f"{c}_mechanism":fav.get("mechanism"),f"{c}_repulsive_member":1 if k in maps["rep"] else 0,
                f"{c}_repulsion_state_requirement":rep.get("state_requirement"),f"{c}_repulsion_state_availability":rep.get("state_availability"),
                f"{c}_repulsion_score":rep.get("repulsion_score"),f"{c}_repulsion_propensity":rep.get("repulsion_propensity"),
                f"{c}_repulsion_ebase_metadata_kcal_mol":rep.get("ebase_metadata_kcal_mol"),f"{c}_repulsion_mechanism":rep.get("mechanism")
            })
        rows.append(row)
    return rows,headers


def _flatten_v2_patch(p: dict) -> dict:
    center=p.get("center_residue") or {}
    row={k:v for k,v in p.items() if k not in {"center_residue","members","seed_members","rin_context"}}
    row.update({"center_res_name":center.get("res_name"),"center_res_seq":center.get("res_seq"),"center_chain":center.get("chain"),
                "members":";".join(str(x) for x in p.get("members",[])),"seed_members":";".join(str(x) for x in p.get("seed_members",[]))})
    rc=p.get("rin_context") or {}
    for k,v in rc.items():
        if k != "residue_metrics": row[f"rin_{k}"]=v
    return row


def _save_property_track_png(path: Path, title: str, x: list[int], propensity: list[float], persistence: list[float], repulsion: list[float]) -> None:
    try:
        import matplotlib.pyplot as plt
        fig=plt.figure(figsize=(12,4.2)); ax=fig.add_subplot(111)
        ax.plot(x,propensity,label="Residue propensity")
        ax.plot(x,persistence,label="Multiscale persistence")
        if any(abs(v)>1e-12 for v in repulsion): ax.plot(x,repulsion,label="Repulsion")
        ax.set_xlabel("Residue sequence number"); ax.set_ylabel("Relative map value (0-100)"); ax.set_ylim(0,105)
        ax.set_title(title); ax.legend(loc="upper right"); fig.tight_layout(); path.parent.mkdir(parents=True,exist_ok=True); fig.savefig(path,dpi=180); plt.close(fig)
    except Exception as exc:
        path.with_suffix(".txt").write_text(f"PNG generation failed: {exc}\n",encoding="utf-8")


def process_case(case: dict, v1) -> dict:
    cid=case["id"]; root=OUTROOT/cid
    for d in ["PDB/original","PDB/prepared","PDB/Bfactor_maps","CSV","CSV/chemistries","JSON","PNG/property_tracks"]: (root/d).mkdir(parents=True,exist_ok=True)
    raw=_fetch_pdb(case["pdb_id"]); prepared,prep=prepare_pdb_text(raw,chain=case.get("chain"))
    (root/"PDB/original"/f"{case['pdb_id']}.pdb").write_text(raw,encoding="utf-8")
    (root/"PDB/prepared"/f"{case['pdb_id']}_prepared.pdb").write_text(prepared,encoding="utf-8")

    req=v1.AnalyzeRequest(pdb_text=prepared,chain=None,env=v1.EnvParams(pH=float(case["pH"]),ionic=150.0,temp=298.0))
    v1res=v1.analyze(req)
    if hasattr(v1res,"body"): raise RuntimeError(f"V1 analyze failed for {cid}")
    v2res=analyze_interface_v2(surface=case["surface"],pdb_text=raw,chain=case.get("chain"),pH=float(case["pH"]),ionic_mM=150.0,temp_K=298.0)
    (root/"JSON"/"full_residue_analysis.json").write_text(json.dumps(_jsonable(v1res),indent=2),encoding="utf-8")
    (root/"JSON"/"v2_interface_prediction.json").write_text(json.dumps(_jsonable(v2res),indent=2),encoding="utf-8")
    (root/"JSON"/"structure_preparation.json").write_text(json.dumps(_jsonable(prep),indent=2),encoding="utf-8")

    allr=v1res.get("all_residues",[]) or []; surf=v1res.get("surface_residues",[]) or []
    _write_csv(root/"CSV"/"all_residues.csv",allr)
    _write_csv(root/"CSV"/"surface_residues.csv",surf)
    _write_csv(root/"CSV"/"surface_features.csv",_feature_rows(v1res))
    comp,headers=_comprehensive_rows(v1res); _write_csv(root/"CSV"/"all_calculated_properties.csv",comp,headers)

    patch_rows=[_flatten_v2_patch(p) for p in v2res.get("patches",[])]; _write_csv(root/"CSV"/"v2_all_patches.csv",patch_rows)
    primary=[_flatten_v2_patch(p) for p in v2res.get("primary_patches",[])]; _write_csv(root/"CSV"/"v2_pareto_primary_patches.csv",primary)

    all_by_key={_res_key(r):r for r in allr}
    surface_values={k:100.0 if r.get("surface_exposed") else 0.0 for k,r in all_by_key.items()}
    scrsa_values={k:100.0*float(r.get("scrsa",0) or 0) for k,r in all_by_key.items()}
    (root/"PDB/Bfactor_maps"/"surface_exposed_Bfactor.pdb").write_text(_bfactor_pdb(prepared,surface_values,["INTERFACESCOUT SURFACE_EXPOSED BINARY MAP","B-FACTOR STORES 0 OR 100"]),encoding="utf-8")
    (root/"PDB/Bfactor_maps"/"scRSA_Bfactor.pdb").write_text(_bfactor_pdb(prepared,scrsa_values,["INTERFACESCOUT scRSA MAP","B-FACTOR STORES 100 x CLIPPED scRSA"]),encoding="utf-8")

    # Feature B-factor maps.
    frows=_feature_rows(v1res)
    for feature in v1res.get("feature_list",[]):
        vals={r["key"]:100.0*float(r.get(feature,0) or 0) for r in frows}
        (root/"PDB/Bfactor_maps"/f"feature_{feature}_Bfactor.pdb").write_text(_bfactor_pdb(prepared,vals,[f"INTERFACESCOUT FEATURE {feature}","B-FACTOR STORES BINARY 0/100"]),encoding="utf-8")

    inner,outer=MULTISCALE_RADII_A; itag=f"{int(inner)}A"; otag=f"{int(outer)}A"
    chemistry_list=v1res.get("chemistry_list",[])
    for chem in chemistry_list:
        g=v1res.get("chemistries",{}).get(chem,{})
        fav=g.get("residues",[]); rep=g.get("repulsive_residues",[]); pc=g.get("patch_centers",[])
        _write_csv(root/"CSV/chemistries"/f"{chem}_favorable_residues.csv",fav)
        _write_csv(root/"CSV/chemistries"/f"{chem}_repulsive_residues.csv",rep)
        _write_csv(root/"CSV/chemistries"/f"{chem}_patch_centers.csv",pc)
        favmap={_res_key(r):float(r.get("propensity",0) or 0) for r in fav}
        repmap={_res_key(r):float(r.get("repulsion_propensity",0) or 0) for r in rep}
        persmap={str(r.get("center_key")):float(r.get("multiscale_persistence",0) or 0) for r in pc}
        geomap={str(r.get("center_key")):float(r.get("multiscale_geomean",0) or 0) for r in pc}
        for mode,vals,label in [("propensity",favmap,"RESIDUE PROPENSITY"),("persistence",persmap,"MULTISCALE PERSISTENCE"),("repulsion",repmap,"REPULSION PROPENSITY"),("geomean",geomap,"MULTISCALE GEOMEAN")]:
            (root/"PDB/Bfactor_maps"/f"{chem}_{mode}_Bfactor.pdb").write_text(_bfactor_pdb(prepared,vals,[f"INTERFACESCOUT CHEMISTRY {chem}",f"DISPLAY_MODE {mode}",f"B-FACTOR STORES {label} (0-100)",f"MULTISCALE RADII {inner}/{outer} A"]),encoding="utf-8")
        # Sequence property track table and PNG.
        seq=sorted(allr,key=lambda r:(str(r.get("chain")),int(r.get("res_seq",0)),str(r.get("icode",""))))
        track=[]
        for r in seq:
            k=_res_key(r); track.append({"key":k,"chain":r.get("chain"),"res_seq":r.get("res_seq"),"res_name":r.get("res_name"),"propensity":favmap.get(k,0.0),"persistence":persmap.get(k,0.0),"repulsion":repmap.get(k,0.0)})
        _write_csv(root/"CSV/chemistries"/f"{chem}_property_track.csv",track)
        if len({r.get("chain") for r in seq})==1:
            _save_property_track_png(root/"PNG/property_tracks"/f"{chem}_property_track.png",f"{case['protein']} - {chem}",[int(r.get("res_seq",0)) for r in seq],[favmap.get(_res_key(r),0.0) for r in seq],[persmap.get(_res_key(r),0.0) for r in seq],[repmap.get(_res_key(r),0.0) for r in seq])

    # Full descriptive GNM and RIN outputs; these do NOT affect interface ranking.
    old=os.environ.pop("INTERFACESCOUT_VALIDATION_CANONICAL_ONLY",None)
    try:
        gnm=solve_gnm(prepared)
        gr=list(gnm.get("residue_metrics",{}).values()); _write_csv(root/"CSV"/"gnm_residue_metrics.csv",gr)
        contact=[]; nodes=gnm.get("nodes",[]); adj=gnm.get("adjacency")
        if isinstance(adj,np.ndarray):
            for i in range(len(nodes)):
                for j in range(i+1,len(nodes)):
                    if int(adj[i,j])!=0: contact.append({"residue_i":nodes[i]["key"],"residue_j":nodes[j]["key"],"contact":1})
        _write_csv(root/"CSV"/"gnm_contact_matrix_nonzero.csv",contact,["residue_i","residue_j","contact"])
        fluct={r["key"]:float(r.get("normalized_fluctuation",0) or 0) for r in gr}
        (root/"PDB/Bfactor_maps"/"gnm_normalized_fluctuation_Bfactor.pdb").write_text(_bfactor_pdb(prepared,_normalize_0_100(fluct),["INTERFACESCOUT GNM DESCRIPTIVE MAP","B-FACTOR STORES NORMALIZED FLUCTUATION RESCALED TO 0-100","NOT USED IN INTERFACE RANKING"]),encoding="utf-8")
    finally:
        if old is not None: os.environ["INTERFACESCOUT_VALIDATION_CANONICAL_ONLY"]=old

    rin=annotate_rin_percentiles(build_rin(prepared),[_res_key(r) for r in surf])
    rmetrics=list(rin.get("residue_metrics",{}).values()); _write_csv(root/"CSV"/"rin_residue_metrics.csv",rmetrics)
    redges=[]
    for pair,count in rin.get("contact_counts",{}).items():
        a,b=pair.split("|",1); redges.append({"residue_i":a,"residue_j":b,"heavy_atom_contacts":count})
    _write_csv(root/"CSV"/"rin_edges_nonzero.csv",redges,["residue_i","residue_j","heavy_atom_contacts"])
    for field in ["degree_normalized_percentile_all","closeness_percentile_all","betweenness_percentile_all","degree_normalized_percentile_surface","closeness_percentile_surface","betweenness_percentile_surface"]:
        vals={r["key"]:float(r.get(field,0) or 0) for r in rmetrics}
        (root/"PDB/Bfactor_maps"/f"rin_{field}_Bfactor.pdb").write_text(_bfactor_pdb(prepared,vals,[f"INTERFACESCOUT RIN {field}","B-FACTOR STORES PERCENTILE (0-100)","DESCRIPTIVE ONLY; NOT USED IN INTERFACE RANKING"]),encoding="utf-8")

    inventory=[]
    for p in sorted(root.rglob("*")):
        if p.is_file(): inventory.append({"relative_path":str(p.relative_to(root)),"size_bytes":p.stat().st_size})
    _write_csv(root/"file_inventory.csv",inventory,["relative_path","size_bytes"])
    return {"case_id":cid,"protein":case["protein"],"pdb_id":case["pdb_id"],"chain":case.get("chain"),"surface":case["surface"],"pH":case["pH"],"n_files":len(inventory),"n_chemistries":len(chemistry_list),"n_surface_residues":len(surf),"n_v2_patches":len(v2res.get("patches",[]))}


def main(manifest_path: Path=DEFAULT_MANIFEST) -> None:
    manifest=json.loads(manifest_path.read_text(encoding="utf-8"))
    OUTROOT.mkdir(parents=True,exist_ok=True)
    import main as v1
    apply_publication_chemistry(v1)
    summary=[]
    for case in manifest.get("cases",[]):
        print(f"FULL_OUTPUT {case['id']}",flush=True)
        summary.append(process_case(case,v1))
    _write_csv(OUTROOT/"panel_summary.csv",summary)
    (OUTROOT/"panel_manifest_prediction_only.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    (OUTROOT/"README.txt").write_text(
        "InterfaceScout full blinded candidate output package.\n"
        "Experimental ground truth was NOT read by this exporter.\n"
        "Includes comprehensive residue CSVs, chemistry-specific CSVs, all chemistry view B-factor PDBs, V2 patch outputs, original/prepared PDBs, descriptive GNM/RIN outputs, sparse contact tables, and property-track PNGs.\n"
        f"Frozen multiscale radii: {MULTISCALE_RADII_A[0]}/{MULTISCALE_RADII_A[1]} A; coarse patch radius: {COARSE_PATCH_RADIUS_A} A.\n",
        encoding="utf-8")
    print("FULL_OUTPUT_SUMMARY "+json.dumps(summary),flush=True)


if __name__=="__main__": main()
