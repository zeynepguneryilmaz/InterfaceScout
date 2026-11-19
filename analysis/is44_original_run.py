#!/usr/bin/env python3
"""Freeze original InterfaceScout v1 all-11-map outputs for the 44-condition adsorption evidence panel.
One deterministic run per unique PDB/chain/pH, then map results to 44 literature conditions offline.
No experimental residue labels enter the algorithm.
"""
import csv, gzip, hashlib, json, os, sys, time, traceback, urllib.request
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "backend"))
from interface_engine import analyze_all_maps
from model_settings import MODEL_VERSION
root=Path("is44_output")
(root/"maps").mkdir(parents=True,exist_ok=True)
(root/"inputs").mkdir(parents=True,exist_ok=True)
STRUCTURES=[
 ("4F5S","A",[7.4]),("2VUF","A",[7.4]),("2CBA","A",[6.3,8.3]),
 ("1UBQ","A",[7.4]),("1JNJ","A",[7.4]),("3ECA","A",[7.4]),
 ("2HAV","A",[7.4]),("1TGU","A",[7.4]),("3NWV","A",[7.4]),
 ("7RSA","A",[7.4]),("3LN2","A",[7.4]),("1LYZ","A",[7.4]),
 ("1HRC","A",[3.0,7.0,7.4]),("1MBN","A",[7.4]),
 ("2HHB","A",[7.4]),("2PTN","A",[7.4]),("1A5R","A",[7.4]),
 ("2Y0G","A",[7.4]),("256L","A",[7.4]),("1AML","A",[7.4]),
 ("2OED","A",[7.4]),("2GB1","A",[7.4]),("7KWI","A",[4.5,7.4]),
]
def download(pid):
    url=f"https://files.rcsb.org/download/{pid}.pdb"
    errors=[]
    for attempt in range(3):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"InterfaceScout-44-benchmark/1.0 scientific analysis"})
            with urllib.request.urlopen(req,timeout=45) as response: return response.read()
        except Exception as e:
            errors.append(str(e));time.sleep(2*(attempt+1))
    raise RuntimeError("; ".join(errors))
def slim_patch(p):
    keys=("display_rank","pareto_front","center_key","members","seed_members","n_members","n_chemistry_seeds","chemistry_support","mean_accessibility","patch_coherence","orientation_coherence")
    return {k:p.get(k) for k in keys}
rows=[]
for pdb,chain,phs in STRUCTURES:
    try:
        raw=download(pdb)
        sha=hashlib.sha256(raw).hexdigest()
        with gzip.open(root/"inputs"/(pdb+".pdb.gz"),"wb") as f:f.write(raw)
        for ph in phs:
            tic=time.perf_counter();tag=f"{pdb}_{chain}_pH{ph:.1f}"
            try:
                result=analyze_all_maps(pdb_id=pdb,pdb_text=raw.decode("utf-8","replace"),chain=chain,pH=ph,ionic_mM=150.0,temp_K=298.0)
                maps={}
                for name,m in result["maps"].items():
                    patches=m["patches"]
                    unions={}
                    for k in (1,3,5):
                        selected=patches[:k]
                        unions[str(k)]={
                            "centers":[p["center_key"] for p in selected],
                            "residues":sorted({x for p in selected for x in p["members"]}),
                            "seed_residues":sorted({x for p in selected for x in p["seed_members"]}),
                        }
                    maps[name]={"n_patches":len(patches),"top5":[slim_patch(p) for p in patches[:5]],"unions":unions}
                output={"pdb":pdb,"chain":chain,"pH":ph,"ionic_mM":150.0,"temp_K":298.0,"version":result["version"],"source_sha256":sha,"n_maps":result["n_maps"],
                    "structure_preparation":result.get("structure_preparation"),"n_exposed":result["diagnostics"]["n_surface_residues"],
                    "exposed":result["diagnostics"]["surface_residue_keys"],"maps":maps}
                (root/"maps"/(tag+".json")).write_text(json.dumps(output,ensure_ascii=False,separators=(",",":")),encoding="utf-8")
                rows.append(dict(pdb=pdb,chain=chain,pH=ph,version=result["version"],n_maps=result["n_maps"],n_exposed=output["n_exposed"],status="ok",error="",seconds=round(time.perf_counter()-tic,3),sha256=sha))
            except Exception as e:
                rows.append(dict(pdb=pdb,chain=chain,pH=ph,version=MODEL_VERSION,n_maps=0,n_exposed=0,status="error",error=str(e)[:500],seconds=round(time.perf_counter()-tic,3),sha256=sha))
                traceback.print_exc()
    except Exception as e:
        for ph in phs:rows.append(dict(pdb=pdb,chain=chain,pH=ph,version=MODEL_VERSION,n_maps=0,n_exposed=0,status="download_error",error=str(e)[:500],seconds=0,sha256=""))
        traceback.print_exc()
with (root/"run_manifest.csv").open("w",newline="",encoding="utf-8") as f:
    writer=csv.DictWriter(f,fieldnames=["pdb","chain","pH","version","n_maps","n_exposed","status","error","seconds","sha256"])
    writer.writeheader();writer.writerows(rows)
summary={"method":"original-InterfaceScout-v1, fixed weights/no label tuning","model_version":MODEL_VERSION,"parameters":{"SASA_probe_A":1.4,"SASA_points_per_atom":200,"scRSA_cutoff":0.05,"radii_A":[6,9],"patch_radius_A":8,"ionic_mM":150,"temperature_K":298},
"unique_pdb":len(STRUCTURES),"unique_protein_pH":len(rows),"succeeded":sum(r["status"]=="ok" for r in rows),"failed":sum(r["status"]!="ok" for r in rows),
"warning":"Material chemistry channels are qualitative; particle size, actual adsorption orientation, contact affinity and binding energy are not modeled. NMR ensemble PDB uses deterministic source preparation. This is not 44 independent experiments."}
(root/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(json.dumps(summary))
