"""Compare three explicit GLY policies on frozen 48 assemblies and all 11 maps.
Does not mutate production algorithm or use adsorption labels.
Baseline: existing CA proxy / 32.28 A2. Alternative: full GLY residue SASA /
104 A2 (Tien et al. PLoS ONE 2013). Ablation: omit GLY.
"""
import sys, csv, json, hashlib, argparse, traceback
from pathlib import Path
from io import StringIO
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"backend"),str(ROOT/"validation")]
import core,geometry,strict_apo_compat_48 as prep,sensitivity_48 as sens
OUT=Path("glycine_three_policy_48");OUT.mkdir(exist_ok=True)
GLY_REF=104.0
POLICIES=("ca_proxy","whole_residue","exclude")
METRICS=("top5_members","centers","top1_members","front1_centers")
def run_case(item):
    pid=item["pdb_id"].upper()
    prep.prepare(pid)
    cif=prep.OUT/f"{pid}-assembly1.cif"
    sourcehash=hashlib.sha256(cif.read_bytes()).hexdigest()
    if sourcehash!=item["assembly_sha256"]:raise ValueError("Frozen assembly SHA mismatch")
    pdb=(prep.OUT/f"{pid}_prepared.pdb").read_text(encoding="utf-8")
    nodes=geometry.extract_ca_nodes(pdb)
    original=core._surface_residues(pdb,7.4)
    atoms=PDBParser(QUIET=True).get_structure("gly_residue",StringIO(pdb))
    ShrakeRupley(probe_radius=core.SASA_PROBE_A,n_points=core.SASA_POINTS).compute(atoms,level="A")
    full={}
    for chain in next(atoms.get_models()):
        for res in chain:
            if res.get_resname().strip().upper()=="GLY" and "CA" in res:
                key=f"{chain.id}:{int(res.id[1])}:{str(res.id[2]).strip()}"
                full[key]=sum(float(a.sasa) for a in res.get_atoms() if a.element!="H")
    base_gly={r["key"] for r in original if r["res_name"]=="GLY"}
    if set(full)!=base_gly:raise ValueError(f"Gly key mismatch {set(full)^base_gly}")
    alternative=[]
    gly_records=[]
    for r in original:
        if r["res_name"]!="GLY":
            alternative.append(r.copy());continue
        alt=r.copy()
        asa=full[r["key"]]
        raw=asa/GLY_REF
        alt["sidechain_sasa"]=round(asa,3)
        alt["scrsa_raw"]=round(raw,5)
        alt["scrsa"]=round(min(max(raw,0.),1.),5)
        alt["surface_exposed"]=bool(raw>=.05)
        alternative.append(alt)
        gly_records.append({"protein":pid,"key":r["key"],"ca_proxy_A2":r["sidechain_sasa"],
            "ca_proxy_scrsa":r["scrsa_raw"],"whole_residue_A2":round(asa,3),
            "whole_residue_RSA":round(raw,5),"ca_exposed":int(r["surface_exposed"]),
            "whole_exposed":int(alt["surface_exposed"])})
    policies={"ca_proxy":original,"whole_residue":alternative,
              "exclude":[r.copy() for r in original if r["res_name"]!="GLY"]}
    calculated={}
    for label,rr in policies.items():
        maps,surface=sens.run_cached(pdb,nodes,rr,.05,(6.,9.),8.,7.4)
        calculated[label]=(maps,surface)
    rows=[]
    for channel in sens.CHANNELS:
        for a,b in (("ca_proxy","whole_residue"),("ca_proxy","exclude"),("whole_residue","exclude")):
            left,leftsurface=calculated[a];right,rightsurface=calculated[b]
            al=left[channel];bl=right[channel]
            row={"protein":pid,"channel":channel,"policy_a":a,"policy_b":b,
                 "n_gly":len(base_gly),"n_gly_exposure_changed":sum(x["ca_exposed"]!=x["whole_exposed"] for x in gly_records),
                 "surface_jaccard":sens.jacc(leftsurface,rightsurface),
                 "patch_count_a":al["n_patches"],"patch_count_b":bl["n_patches"]}
            for key in METRICS:
                row[key+"_jaccard"]=sens.jacc(al[key],bl[key])
            rows.append(row)
    return rows,gly_records,{"protein":pid,"status":"PASS","source_sha256":sourcehash,
       "n_gly":len(base_gly),"n_gly_exposure_changed":sum(x["ca_exposed"]!=x["whole_exposed"] for x in gly_records)}
def write_csv(file,rows):
    if not rows:return
    with file.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
def shard(i,n):
    all_rows=[];gly=[];provenance=[];failed=[]
    for j,item in enumerate(sens.MANIFEST):
        if j%n!=i:continue
        try:
            rows,g,meta=run_case(item)
            all_rows+=rows;gly+=g;provenance.append(meta)
            print("GLY_THREE_OK",json.dumps(meta),flush=True)
        except Exception as e:
            fail={"protein":item["pdb_id"],"error":str(e),"trace":traceback.format_exc()[-1500:]}
            failed.append(fail);print("GLY_THREE_FAIL",json.dumps(fail),flush=True)
    write_csv(OUT/f"maps_{i}.csv",all_rows)
    write_csv(OUT/f"gly_{i}.csv",gly)
    (OUT/f"provenance_{i}.json").write_text(json.dumps(provenance,indent=2))
    (OUT/f"failures_{i}.json").write_text(json.dumps(failed,indent=2))
    print("THREE_SHARD",i,"pass",len(provenance),"fail",len(failed),flush=True)
    if failed:raise SystemExit(2)
def aggregate(path):
    import pandas as pd
    source=sorted(path.glob("maps_[0-9].csv"))
    if len(source)!=4:raise ValueError(f"Expected four shards, found {source}")
    df=pd.concat([pd.read_csv(x) for x in source],ignore_index=True)
    if len(df)!=48*11*3 or df.protein.nunique()!=48 or df.channel.nunique()!=11:
        raise ValueError("Incomplete three-way 48x11 comparisons")
    g=pd.concat([pd.read_csv(x) for x in sorted(path.glob("gly_[0-9].csv"))],ignore_index=True)
    df.to_csv(path/"all_1584_map_comparisons.csv",index=False)
    g.to_csv(path/"glycine_residue_comparison.csv",index=False)
    cols=["surface_jaccard"]+[k+"_jaccard" for k in METRICS]
    means=df.groupby(["policy_a","policy_b"])[cols].mean().reset_index()
    means.to_csv(path/"policy_means.csv",index=False)
    per_channel=df.groupby(["policy_a","policy_b","channel"])[cols].mean().reset_index()
    per_channel.to_csv(path/"channel_means.csv",index=False)
    rng=np.random.default_rng(20261123)
    ci=[]
    for (a,b),group in df.groupby(["policy_a","policy_b"]):
        agg=group.groupby("protein")[cols].mean()
        for m in cols:
            v=agg[m].to_numpy()
            draws=np.mean(rng.choice(v,(3000,len(v)),replace=True),axis=1)
            ci.append({"policy_a":a,"policy_b":b,"metric":m,"mean":float(v.mean()),
              "ci95_low":float(np.quantile(draws,.025)),"ci95_high":float(np.quantile(draws,.975))})
    pd.DataFrame(ci).to_csv(path/"protein_cluster_bootstrap.csv",index=False)
    result={"proteins":48,"maps":11,"comparisons":len(df),"gly_residue_records":len(g),
      "changed_gly_exposure_count":int((g.ca_exposed!=g.whole_exposed).sum()),
      "policy_means":means.to_dict("records"),
      "changed_top5_by_policy":df.groupby(["policy_a","policy_b"])["top5_members_jaccard"].apply(lambda v:int((v<1-1e-9).sum())).to_dict().__str__(),
      "interpretation_limit":"Whole-residue GLY RSA uses Tien 104 A2; other residues retain original sidechain RSA, therefore this is a GLY-only policy sensitivity, not a uniform whole-residue RSA framework."}
    (path/"SUMMARY.json").write_text(json.dumps(result,indent=2))
    print("GLY_THREE_FINAL",json.dumps(result),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--shard",type=int);p.add_argument("--shards",type=int,default=4)
    p.add_argument("--aggregate",type=Path);a=p.parse_args()
    if a.aggregate:aggregate(a.aggregate)
    elif a.shard is not None:shard(a.shard,a.shards)
    else:p.error("Specify --shard or --aggregate")
