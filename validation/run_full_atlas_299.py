"""Rebuild the full InterfaceScout Top-k residue atlas from the archived, frozen 300-PDB selection.
No protein structures are resampled; all reported map metrics are checked against the original run.
"""
from __future__ import annotations
import csv, gzip, hashlib, io, json, os, sys, zipfile
from collections import Counter
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"backend"))
from interface_engine import analyze_all_maps, PUBLIC_MAP_SPECS

ARCHIVE=Path(os.environ.get("INPUT_300_ZIP","validation/input_300_proteins.zip"))
OUT=Path("validation/full_atlas_299_results")
OUT.mkdir(parents=True,exist_ok=True)
if not ARCHIVE.is_file():raise FileNotFoundError(ARCHIVE)
EXPECTED_N=299
MAPS=[s[0] for s in PUBLIC_MAP_SPECS]

def readcsv(z,name):return list(csv.DictReader(io.StringIO(z.read(name).decode("utf-8-sig"))))
def writecsv(path,rows):
 if not rows:raise ValueError("Empty output "+str(path))
 with path.open("w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)
def normalized_member_list(patches,k):
 keys=set()
 for p in patches[:k]:keys.update(str(x) for x in p.get("members",[]))
 return sorted(keys)
def res_record(member):
 return {"key":member.get("key"),"chain":member.get("chain"),"res_seq":member.get("res_seq"),"icode":member.get("icode"),"res_name":member.get("res_name")}
manifest_out=[];toprows=[];maprows=[];patchrows=[];errors=[]
with zipfile.ZipFile(ARCHIVE) as z:
 manifest=readcsv(z,"manifest.csv")
 expected={(r["pdb_id"],r["chain"],r["map"]):r for r in readcsv(z,"map_metrics.csv")}
 if len(manifest)!=300:raise ValueError("Expected frozen panel of 300 entries")
 valid=[r for r in manifest if r["pdb_id"]!="3IYK"]
 if len(valid)!=299:raise ValueError("Expected exactly one excluded 3IYK structure")
 for idx,row in enumerate(valid,1):
  pid=row["pdb_id"]; chain=row["chain"]
  path=f"pdb/{pid}_{chain}.pdb"
  data=z.read(path)
  if hashlib.sha256(data).hexdigest()!=row["sha256"]:
   raise RuntimeError(f"PDB SHA mismatch: {pid}/{chain}")
  text=data.decode("utf-8","replace")
  try:
   result=analyze_all_maps(pdb_text=text,chain=chain,pH=7.4,ionic_mM=150,temp_K=298)
   if result["n_maps"]!=11:raise ValueError("Expected 11 maps")
   surf=result["diagnostics"]["n_surface_residues"]
   for chem in MAPS:
    m=result["maps"][chem];patches=m["patches"]
    original=expected[(pid,chain,chem)]
    top5=normalized_member_list(patches,5)
    checks={"n_patches":(len(patches),int(original["n_patches"])),"top5_n_residues":(len(top5),int(original["top5_n_residues"])),"n_surface":(surf,int(original["n_surface"]))}
    mismatches={key:pair for key,pair in checks.items() if pair[0]!=pair[1]}
    if mismatches:raise ValueError(f"{pid} {chem}: frozen map mismatch {mismatches}")
    maprows.append({"pdb_id":pid,"chain":chain,"chemistry":chem,"n_patches":len(patches),"n_surface_residues":surf,"top5_residues":len(top5),"input_sha256":row["sha256"]})
    for k in (1,3,5):
     selected=patches[:k]
     keys=normalized_member_list(patches,k)
     toprows.append({"pdb_id":pid,"chain":chain,"chemistry":chem,"top_k":k,"n_patches_used":len(selected),"n_union_residues":len(keys),"residue_keys":";".join(keys),"patch_centers":";".join(str(p["center_key"]) for p in selected),"pareto_fronts":";".join(str(p["pareto_front"]) for p in selected)})
    for p in patches:
     patchrows.append({"pdb_id":pid,"chain":chain,"chemistry":chem,"display_rank":p["display_rank"],"pareto_front":p["pareto_front"],"center_key":p["center_key"],"n_members":p["n_members"],"member_keys":";".join(str(k) for k in p["members"]),"seed_keys":";".join(str(k) for k in p["seed_members"]),"chemistry_support":p["chemistry_support"],"mean_accessibility":p["mean_accessibility"],"patch_coherence":p["patch_coherence"],"orientation_coherence":p["orientation_coherence"]})
   manifest_out.append({"pdb_id":pid,"chain":chain,"n_residues":row["n_residues"],"method":row["method"],"n_surface_residues":surf,"original_sha256":row["sha256"]})
   # JSON file per protein enables lazy-loaded 3D atlas; include raw PDB for offline viewer.
   payload={"pdb_id":pid,"chain":chain,"n_residues":int(row["n_residues"]),"pdb":text,"maps":{}}
   for chem in MAPS:
    ps=result["maps"][chem]["patches"]
    payload["maps"][chem]={"n_patches":len(ps),"top1":normalized_member_list(ps,1),"top3":normalized_member_list(ps,3),"top5":normalized_member_list(ps,5),"top5_centers":[p["center_key"] for p in ps[:5]],"top5_pareto_fronts":[p["pareto_front"] for p in ps[:5]]}
   dest=OUT/"atlas"/(pid+"_"+chain+".json")
   dest.parent.mkdir(parents=True,exist_ok=True)
   dest.write_text(json.dumps(payload,separators=(",",":"),ensure_ascii=False))
   print("ATLAS_PROTEIN_OK",idx,EXPECTED_N,pid,chain,flush=True)
  except Exception as exc:
   errors.append({"pdb_id":pid,"chain":chain,"error":repr(exc)})
   print("ATLAS_PROTEIN_ERROR",pid,chain,repr(exc),flush=True)
writecsv(OUT/"atlas_manifest.csv",manifest_out)
writecsv(OUT/"map_qc.csv",maprows)
writecsv(OUT/"top_k_residue_unions.csv",toprows)
writecsv(OUT/"all_ranked_patches.csv",patchrows)
if errors:writecsv(OUT/"errors.csv",errors)
summary={"n_input":300,"n_excluded":1,"excluded":[{"pdb_id":"3IYK","reason":"C-alpha-only representation; incompatible with side-chain SASA"}],"n_proteins":len(manifest_out),"n_maps":len(maprows),"n_topk_rows":len(toprows),"n_ranked_patches":len(patchrows),"n_errors":len(errors),"frozen_original_metric_agreement":len(errors)==0,"chemistry_maps":MAPS,"source":"GitHub Actions run 43869361584 artifact 14226897455","method_model":"2.4.0-selected-parameters","interpretation":"Candidate surface regions only, no experimental adsorption pose validation"}
(OUT/"atlas_summary.json").write_text(json.dumps(summary,indent=2))
print("ATLAS_COMPLETE",json.dumps(summary),flush=True)
if len(manifest_out)!=299 or len(maprows)!=3289 or len(toprows)!=9867 or errors:raise RuntimeError("Atlas completeness verification failed")
