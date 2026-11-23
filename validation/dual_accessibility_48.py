"""Experimental dual-accessibility model for the frozen 48-structure cohort.
All-residue RSA gates surface geometry; non-GLY side-chain RSA informs chemistry.
GLY is geometrically retained with zero side-chain chemistry. No production edits.
Tien et al. 2013 theoretical maximum ASA values, Table 1, DOI:
10.1371/journal.pone.0080635. 1.4 A water probe and 200 points/atom.
Surface threshold 0.05 is a sensitivity setting, NOT yet calibrated for tRSA.
"""
import sys,csv,json,hashlib,argparse,traceback
from io import StringIO
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"backend"),str(ROOT/"validation")]
import core,geometry,coarse_patch,strict_apo_compat_48 as prep,sensitivity_48 as sens
OUT=Path("dual_accessibility_48")
OUT.mkdir(exist_ok=True)
MAX_ASA_TIEN={
 "ALA":129.,"ARG":274.,"ASN":195.,"ASP":193.,"CYS":167.,
 "GLN":225.,"GLU":223.,"GLY":104.,"HIS":224.,"ILE":197.,
 "LEU":201.,"LYS":236.,"MET":224.,"PHE":240.,"PRO":159.,
 "SER":155.,"THR":172.,"TRP":285.,"TYR":263.,"VAL":174.
}
def dual_rows(pdb,baseline):
    structure=PDBParser(QUIET=True).get_structure("dual_accessibility",StringIO(pdb))
    ShrakeRupley(probe_radius=core.SASA_PROBE_A,n_points=core.SASA_POINTS).compute(structure,level="A")
    total={}
    for ch in next(structure.get_models()):
        for res in ch:
            name=res.resname.strip().upper()
            if name not in MAX_ASA_TIEN or "CA" not in res:continue
            key=f"{ch.id}:{int(res.id[1])}:{str(res.id[2]).strip()}"
            atoms=[a for a in res.get_atoms() if (getattr(a,"element","") or "").upper()!="H"]
            total[key]=sum(float(a.sasa) for a in atoms)
    if set(total)!={r["key"] for r in baseline}:raise ValueError("Residue-key mismatch")
    updated=[]; per_res=[]
    for row in baseline:
        rn=row["res_name"]; key=row["key"]; x=row.copy()
        tASA=total[key]; trsa=tASA/MAX_ASA_TIEN[rn]
        x["residue_sasa"]=round(tASA,3)
        x["trsa_raw"]=round(trsa,5)
        x["trsa"]=round(min(max(trsa,0),1),5)
        x["surface_exposed"]=bool(trsa>=.05)
        # GLY is present geometrically but never receives a pseudo-side-chain.
        if rn=="GLY":
            x["sidechain_sasa"]=0.0
            x["scrsa_raw"]=0.0
            x["scrsa"]=0.0
        updated.append(x)
        per_res.append({"key":key,"residue":rn,"tASA_A2":round(tASA,3),
            "tRSA":round(trsa,5),"sASA_A2":x["sidechain_sasa"],
            "scRSA":x["scrsa"],"surface_total_gate":int(x["surface_exposed"]),
            "surface_legacy_gate":int(row["surface_exposed"])})
    return updated,per_res

def run_dual(pdb,rows,nodes):
    # Unlike sensitivity_48.run_cached, never re-gate using side-chain RSA.
    # This is a strict two-channel calculation: total RSA gates, sidechain
    # RSA describes chemistry; glycine has zero chemistry but remains in distances.
    previous=core._surface_residues
    try:
        core._surface_residues=lambda pdb_text,pH:[r.copy() for r in rows]
        request=core.AnalyzeRequest(pdb_text=pdb,env=core.EnvParams(pH=7.4,ionic=150.,temp=298.))
        result=core.analyze(request)
    finally:
        core._surface_residues=previous
    maps={}
    for channel in sens.CHANNELS:
        patches=coarse_patch.build_coarse_patches(v1_result=result,chemistry=channel,ca_nodes=nodes)
        maps[channel]=sens.summarize_map(patches,result["chemistries"][channel]["residues"])
    return maps,{str(r["key"]) for r in result["surface_residues"]}

def compute(item):
    pid=item["pdb_id"].upper()
    prep.prepare(pid)
    cif=prep.OUT/f"{pid}-assembly1.cif"
    checksum=hashlib.sha256(cif.read_bytes()).hexdigest()
    if checksum!=item["assembly_sha256"]:raise ValueError("Frozen assembly mismatch")
    pdb=(prep.OUT/f"{pid}_prepared.pdb").read_text()
    nodes=geometry.extract_ca_nodes(pdb)
    baseline=core._surface_residues(pdb,7.4)
    dual,res=dual_rows(pdb,baseline)
    bm,bs=sens.run_cached(pdb,nodes,baseline,.05,(6.,9.),8.,7.4)
    dm,ds=run_dual(pdb,dual,nodes)
    assert {r["key"] for r in dual if r["surface_exposed"]}==ds
    gly={r["key"] for r in dual if r["res_name"]=="GLY" and r["surface_exposed"]}
    if not gly.issubset(ds):raise ValueError("Glycine omitted from surface geometry")
    for ch in sens.CHANNELS:
        # A glycine residue must not appear as favorable chemistry.
        records=dm[ch]["patches"]
        if any(p.get("center_key") not in ds for p in records):
            raise ValueError("Patch center outside geometric surface")
    comparisons=[]
    for ch in sens.CHANNELS:
        a,b=bm[ch],dm[ch]
        comparisons.append({"protein":pid,"channel":ch,"baseline_surface_count":len(bs),
            "dual_surface_count":len(ds),"gly_surface_geometry":len(gly),
            "surface_jaccard":sens.jacc(bs,ds),
            "top5_members_jaccard":sens.jacc(a["members"],b["members"]),
            "top5_centers_jaccard":sens.jacc(a["centers"],b["centers"]),
            "top1_members_jaccard":sens.jacc(a["top1_members"],b["top1_members"]),
            "front1_centers_jaccard":sens.jacc(a["front1_centers"],b["front1_centers"]),
            "baseline_n_patches":a["n_patches"],"dual_n_patches":b["n_patches"]})
    for rr in res:rr["protein"]=pid
    return comparisons,res,{"protein":pid,"source_sha256":checksum,"status":"PASS",
        "n_residues":len(res),"n_gly":sum(r["residue"]=="GLY" for r in res),
        "n_gly_geometrically_exposed":len(gly)}
def save(path,rows):
    if not rows:return
    with path.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
def shard(i,n):
    maps=[];residues=[];records=[];errors=[]
    for index,item in enumerate(sens.MANIFEST):
        if index%n!=i:continue
        try:
            a,b,c=compute(item)
            maps+=a;residues+=b;records.append(c)
            print("DUAL_PASS",json.dumps(c),flush=True)
        except Exception as exc:
            error={"protein":item["pdb_id"],"error":repr(exc),"trace":traceback.format_exc()[-2000:]}
            errors.append(error);print("DUAL_FAIL",json.dumps(error),flush=True)
    save(OUT/f"maps_{i}.csv",maps);save(OUT/f"residues_{i}.csv",residues)
    (OUT/f"provenance_{i}.json").write_text(json.dumps(records,indent=2))
    (OUT/f"errors_{i}.json").write_text(json.dumps(errors,indent=2))
    if errors:raise SystemExit(2)
def aggregate(path):
    import pandas as pd
    files=sorted(path.glob("maps_[0-9].csv"))
    if len(files)!=4:raise ValueError("Missing shards")
    df=pd.concat([pd.read_csv(p) for p in files],ignore_index=True)
    if len(df)!=528 or df.protein.nunique()!=48 or df.channel.nunique()!=11:
        raise ValueError("Incomplete 48 x 11 data")
    rr=pd.concat([pd.read_csv(p) for p in sorted(path.glob("residues_[0-9].csv"))],ignore_index=True)
    df.to_csv(path/"all_528_comparisons.csv",index=False)
    rr.to_csv(path/"all_residues.csv",index=False)
    metrics=["surface_jaccard","top5_members_jaccard","top5_centers_jaccard",
             "top1_members_jaccard","front1_centers_jaccard"]
    df.groupby("channel")[metrics].mean().to_csv(path/"channel_means.csv")
    rng=np.random.default_rng(20261123)
    byprot=df.groupby("protein")[metrics].mean()
    confidence={}
    for key in metrics:
        values=byprot[key].to_numpy()
        boot=np.mean(rng.choice(values,(3000,len(values)),replace=True),axis=1)
        confidence[key]={"mean":float(values.mean()),
           "ci95":[float(np.quantile(boot,.025)),float(np.quantile(boot,.975))]}
    report={"n_proteins":48,"n_maps":11,"n_map_comparisons":len(df),
        "n_residues":len(rr),"n_gly":int((rr.residue=="GLY").sum()),
        "n_gly_geometrically_exposed":int(((rr.residue=="GLY")&(rr.surface_total_gate==1)).sum()),
        "surface_gate_changes":int((rr.surface_total_gate!=rr.surface_legacy_gate).sum()),
        "protein_cluster_CI":confidence,
        "interpretation":"Exploratory total-RSA gate threshold 0.05; distinct from original scRSA threshold. No experimental adsorption outcomes evaluated."}
    (path/"SUMMARY.json").write_text(json.dumps(report,indent=2))
    print("DUAL_FINAL",json.dumps(report),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--shard",type=int)
    p.add_argument("--shards",type=int,default=4);p.add_argument("--aggregate",type=Path)
    a=p.parse_args()
    if a.aggregate:aggregate(a.aggregate)
    elif a.shard is not None:shard(a.shard,a.shards)
    else:p.error("Choose shard or aggregate")
