"""Paired literature localization benchmark: legacy vs dual-accessibility.
Frozen 9 condition manifest; same structure, pH and all 11 channels. Output
per-case/channel paired scores; broad experimental ranges are not exact contacts.
This is exploratory cross-model evaluation, not an independent holdout after
model choice if manuscript labels were previously inspected.
"""
import sys,csv,json,hashlib,traceback,argparse,urllib.request
from pathlib import Path
from io import StringIO
import numpy as np
import gemmi
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"backend"),str(ROOT/"validation")]
import core,geometry,coarse_patch,sensitivity_48 as sens,dual_accessibility_48 as dual
MANIFEST=json.loads((ROOT/"validation"/"benchmark_manifest.json").read_text())
OUT=Path("dual_literature_benchmark");OUT.mkdir(exist_ok=True)
ASSEMBLY={"3ECA","1TGU"}
CANON=set(core.SIDECHAIN_REF_ASA)
def prepare(case):
    pid=case["pdb_id"];folder=OUT/pid;folder.mkdir(parents=True,exist_ok=True)
    url=f"https://files.rcsb.org/download/{pid}{'-assembly1' if pid in ASSEMBLY else ''}.cif"
    source=folder/"source.cif"
    if not source.exists():
        with urllib.request.urlopen(url,timeout=100) as r:source.write_bytes(r.read())
    st=gemmi.read_structure(str(source))
    st.remove_alternative_conformations();st.remove_hydrogens()
    want=[x.strip() for x in case["chain"].split(",")]
    chosen=list(st[0]) if pid in ASSEMBLY else [c for c in st[0] if c.name in want]
    if not chosen:raise ValueError(f"missing chain {want} among {[x.name for x in st[0]]}")
    clean=gemmi.Structure();clean.name=pid;clean.cell=st.cell;clean.spacegroup_hm=st.spacegroup_hm
    model=gemmi.Model(1);mapping=[]
    for i,ch in enumerate(chosen):
        new=gemmi.Chain(chr(65+i))
        for res in ch.get_polymer():
            if res.name in CANON:new.add_residue(res.clone())
        if len(new)==0:continue
        mapping.append({"original":ch.name,"prepared":new.name,"residues":len(new)})
        model.add_chain(new)
    clean.add_model(model)
    pdb=clean.make_pdb_string()
    if not pdb.count("\nATOM"):raise ValueError("empty protein-only output")
    (folder/"prepared.pdb").write_text(pdb)
    return pdb,{"source_url":url,"source_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
                "prepared_sha256":hashlib.sha256(pdb.encode()).hexdigest(),"mapping":mapping}
def compute_maps(pdb,pH,rows):
    previous=core._surface_residues
    try:
        core._surface_residues=lambda text,pH_input:[r.copy() for r in rows]
        response=core.analyze(core.AnalyzeRequest(pdb_text=pdb,env=core.EnvParams(pH=pH,ionic=150.,temp=298.)))
    finally:core._surface_residues=previous
    nodes=geometry.extract_ca_nodes(pdb)
    result={}
    for channel in sens.CHANNELS:
        patches=coarse_patch.build_coarse_patches(v1_result=response,chemistry=channel,ca_nodes=nodes)
        result[channel]=patches[:5]
    return result
def labels_for_case(case,rows):
    eligible=[r for r in rows if case["gt_kind"]=="equivalent_chain_ranges" or
              r["chain"] in {x.strip() for x in case["chain"].split(",")}]
    if case["gt_kind"]=="exact":
        positions=set(case["gt_exact"])
        labels={r["key"] for r in eligible if r["res_seq"] in positions}
    else:
        spans=case["gt_ranges"]
        labels={r["key"] for r in eligible if any(lo<=r["res_seq"]<=hi for lo,hi in spans)}
    if not labels:raise ValueError(f"no mapped benchmark labels {case['id']}")
    return labels
def comparison(case):
    pdb,prov=prepare(case)
    pH=float(case["pH"])
    baseline=core._surface_residues(pdb,pH)
    double,extras=dual.dual_rows(pdb,baseline)
    lab=labels_for_case(case,baseline)
    coords={r["key"]:np.array((r["x"],r["y"],r["z"])) for r in baseline}
    old=compute_maps(pdb,pH,baseline)
    new=compute_maps(pdb,pH,double)
    records=[]
    for ch in sens.CHANNELS:
        for name,patches in (("legacy",old[ch]),("dual",new[ch])):
            chosen=patches[:5]
            centers=[p.get("center_key") for p in chosen]
            centers=[k for k in centers if k in coords]
            members=set(str(k) for p in chosen for k in p.get("members",[]))
            near={k for k in lab if any(np.linalg.norm(coords[k]-coords[c])<=8. for c in centers)}
            first=next((rank for rank,p in enumerate(chosen,1) if
               p.get("center_key") in coords and any(np.linalg.norm(coords[k]-coords[p["center_key"]])<=8. for k in lab)),None)
            records.append({"case":case["id"],"analysis_set":case["analysis_set"],
                "pdb_id":case["pdb_id"],"protein":case["protein"],"channel":ch,
                "method":name,"pH":pH,"n_labels":len(lab),"n_top5":len(chosen),
                "top5_hit_8A":int(bool(near)),"top5_recall_8A":len(near)/len(lab),
                "top5_direct_overlap":int(bool(members & lab)),
                "top5_direct_recall":len(members & lab)/len(lab),
                "first_near_rank":first if first is not None else 999,
                "n_map_members":len(members),"n_near_label_residues":len(near),
                "label_kind":case["gt_kind"],"evidence_resolution":case["evidence_resolution"]})
    prov.update({"case":case["id"],"n_labeled_residues":len(lab),"pH":pH})
    return records,prov
def save(path,items):
    if items:
        with path.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(items[0]));w.writeheader();w.writerows(items)
def shard(i,n):
    rows=[];prov=[];failed=[]
    for j,item in enumerate(MANIFEST["cases"]):
        if j%n!=i:continue
        try:
            a,b=comparison(item);rows+=a;prov.append(b)
            print("BENCH_PASS",item["id"],len(a),flush=True)
        except Exception as e:
            fail={"case":item["id"],"error":str(e),"trace":traceback.format_exc()[-1800:]}
            failed.append(fail);print("BENCH_FAIL",json.dumps(fail),flush=True)
    save(OUT/f"rows_{i}.csv",rows)
    (OUT/f"provenance_{i}.json").write_text(json.dumps(prov,indent=2))
    (OUT/f"errors_{i}.json").write_text(json.dumps(failed,indent=2))
    if failed:raise SystemExit(2)
def aggregate(path):
    import pandas as pd
    files=sorted(path.glob("rows_[0-9].csv"))
    if len(files)!=3:raise ValueError("Expected three shards")
    d=pd.concat([pd.read_csv(f) for f in files],ignore_index=True)
    if len(d)!=9*11*2 or d.case.nunique()!=9 or d.channel.nunique()!=11:
        raise ValueError("Incomplete experimental comparison")
    d.to_csv(path/"all_198_condition_channel_scores.csv",index=False)
    metrics=["top5_hit_8A","top5_recall_8A","top5_direct_overlap","top5_direct_recall"]
    per_case=d.groupby(["analysis_set","case","method"])[metrics].mean().reset_index()
    per_case.to_csv(path/"per_case_all_channel_means.csv",index=False)
    delta=per_case.pivot(index=["analysis_set","case"],columns="method",values=metrics)
    deltas=[]
    for key,case in delta.iterrows():
        for metric in metrics:
            deltas.append({"analysis_set":key[0],"case":key[1],"metric":metric,
               "legacy":float(case[(metric,"legacy")]),"dual":float(case[(metric,"dual")]),
               "paired_delta":float(case[(metric,"dual")]-case[(metric,"legacy")])})
    pd.DataFrame(deltas).to_csv(path/"paired_case_deltas.csv",index=False)
    group=d.groupby(["analysis_set","channel","method"])[metrics].mean().reset_index()
    group.to_csv(path/"per_channel_outcomes.csv",index=False)
    report={"conditions":9,"primary_conditions":7,"secondary_conditions":2,
      "independent_proteins":len(set(case["pdb_id"] for case in MANIFEST["cases"])),
      "channels":11,"records":len(d),
      "per_case_method_means":per_case.to_dict("records"),
      "primary_summary":per_case[per_case.analysis_set=="primary"].groupby("method")[metrics].mean().to_dict("index"),
      "secondary_summary":per_case[per_case.analysis_set=="secondary"].groupby("method")[metrics].mean().to_dict("index"),
      "interpretation_limit":"No new holdout: this experimental manifest was previously used to report legacy performance. Broad-region labels and multiple pH cases are not independent atomic contacts. All-channel average is descriptive, not material-specific model selection."}
    (path/"SUMMARY.json").write_text(json.dumps(report,indent=2))
    print("BENCH_FINAL",json.dumps(report),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--shard",type=int);p.add_argument("--shards",type=int,default=3)
    p.add_argument("--aggregate",type=Path);a=p.parse_args()
    if a.aggregate:aggregate(a.aggregate)
    elif a.shard is not None:shard(a.shard,a.shards)
    else:p.error("shard or aggregate")
