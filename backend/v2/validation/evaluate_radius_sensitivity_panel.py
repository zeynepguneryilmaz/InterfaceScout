"""Label-independent radius-pair sensitivity across the 14-protein panel."""
from __future__ import annotations

import json
import urllib.request
from statistics import median

import numpy as np

import main as v1
from v2.chemistry_freeze import apply_publication_chemistry
from v2.prepare import prepare_pdb_text

STRUCTURES = [
    {"id":"albumin","pdb_id":"4F5S","chain":"A","group":"plasma_corona"},
    {"id":"fibrinogen","pdb_id":"3GHG","chain":"A,B,C","group":"plasma_corona"},
    {"id":"igg","pdb_id":"1IGT","chain":"A,B,C,D","group":"plasma_corona"},
    {"id":"transferrin","pdb_id":"3QYT","chain":"A","group":"plasma_corona"},
    {"id":"lysozyme","pdb_id":"1AKI","chain":"A","group":"high_resolution_adsorption_model"},
    {"id":"rnase_a","pdb_id":"7RSA","chain":"A","group":"high_resolution_adsorption_model"},
    {"id":"myoglobin","pdb_id":"1MBN","chain":"A","group":"high_resolution_adsorption_model"},
    {"id":"hcaii","pdb_id":"2CBA","chain":"A","group":"high_resolution_adsorption_model"},
    {"id":"ubiquitin","pdb_id":"1UBQ","chain":"A","group":"high_resolution_adsorption_model"},
    {"id":"acp","pdb_id":"1APS","chain":"A","group":"high_resolution_adsorption_model"},
    {"id":"catalase","pdb_id":"1TGU","chain":"A,B,C,D","group":"diverse_size_function"},
    {"id":"asparaginase","pdb_id":"3ECA","chain":"A,B,C,D","group":"diverse_size_function"},
    {"id":"beta2m","pdb_id":"1JNJ","chain":"A","group":"diverse_size_function"},
    {"id":"cytochrome_c","pdb_id":"1HRC","chain":"A","group":"diverse_size_function"},
]
PAIRS = [(5.0,8.0),(5.0,10.0),(8.0,10.0),(8.0,12.0),(10.0,12.0),(10.0,15.0)]
REF_PAIR = (5.0,8.0)


def fetch_pdb(pid):
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.pdb", timeout=45) as r:
        return r.read().decode("utf-8", errors="replace")


def analyze_base(raw, chain):
    prepared, _ = prepare_pdb_text(raw, chain=chain)
    req = v1.AnalyzeRequest(pdb_text=prepared, chain=None, env=v1.EnvParams(pH=7.4, ionic=150.0, temp=298.0))
    return v1.analyze(req)


def jaccard(a,b):
    u=a|b
    return len(a&b)/len(u) if u else 1.0


def top10_for_pair(surface, distances, chemistry_key, pair):
    meta=v1.CHEMISTRIES[chemistry_key]
    fav=meta["favorable"]
    local=np.zeros(len(surface), dtype=float)
    for i,r in enumerate(surface):
        rn=r["res_name"]
        if rn in fav:
            _,_,state=fav[rn]
            local[i]=float(r["scrsa"])*float(v1.state_availability(rn,state,7.4))
    if not len(surface) or np.max(local)<=0:
        return set()
    dens=[]
    norms=[]
    for R in pair:
        d=(distances<=R).astype(float) @ local
        dens.append(d)
        mx=float(np.max(d)) if d.size else 0.0
        norms.append(d/mx if mx>0 else np.zeros_like(d))
    persistence=100.0*np.minimum(norms[0],norms[1])
    geomean=100.0*np.sqrt(norms[0]*norms[1])
    rows=[]
    for i,r in enumerate(surface):
        if persistence[i]<=0: continue
        rows.append((float(persistence[i]),float(geomean[i]),str(r["key"])))
    rows.sort(key=lambda x:(-x[0],-x[1],x[2]))
    return {k for _,_,k in rows[:10]}


def main():
    apply_publication_chemistry(v1)
    rows=[]
    for spec in STRUCTURES:
        try:
            out=analyze_base(fetch_pdb(spec["pdb_id"]), spec["chain"])
            surface=out.get("surface_residues",[])
            coords=np.asarray([[r["x"],r["y"],r["z"]] for r in surface],dtype=float)
            dif=coords[:,None,:]-coords[None,:,:]
            distances=np.sqrt(np.sum(dif*dif,axis=2)) if len(coords) else np.empty((0,0))
            ref={c:top10_for_pair(surface,distances,c,REF_PAIR) for c in v1.CHEMISTRIES}
            by_pair={}
            for pair in PAIRS:
                cur={c:top10_for_pair(surface,distances,c,pair) for c in v1.CHEMISTRIES}
                vals=[jaccard(ref[c],cur[c]) for c in v1.CHEMISTRIES]
                by_pair[f"{int(pair[0])}/{int(pair[1])}"]=float(median(vals))
            row={**spec,"status":"ok","n_surface":len(surface),"median_top10_hotspot_jaccard_vs_5_8":by_pair}
        except Exception as exc:
            row={**spec,"status":"failed","error":f"{type(exc).__name__}: {exc}"}
        rows.append(row)
        print(json.dumps(row,indent=2,sort_keys=True))
    good=[r for r in rows if r["status"]=="ok"]
    summary={"n_requested":len(rows),"n_successful":len(good),"n_failed":len(rows)-len(good),"panel_median_top10_hotspot_jaccard_vs_5_8":{}}
    for pair in PAIRS:
        k=f"{int(pair[0])}/{int(pair[1])}"
        summary["panel_median_top10_hotspot_jaccard_vs_5_8"][k]=float(median(r["median_top10_hotspot_jaccard_vs_5_8"][k] for r in good)) if good else None
    payload={"purpose":"label-independent radius-pair sensitivity across the 14-protein panel","reference_pair_A":[5.0,8.0],"tested_pairs_A":[list(x) for x in PAIRS],"structures":rows,"summary":summary}
    with open("radius_sensitivity_panel.json","w",encoding="utf-8") as fh: json.dump(payload,fh,indent=2,sort_keys=True)
    print(json.dumps(summary,indent=2,sort_keys=True))

if __name__=="__main__": main()
