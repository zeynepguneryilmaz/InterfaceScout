from __future__ import annotations
import csv, json, math, sys, urllib.request
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
HERE=ROOT/"validation"/"level3_expansion_2026-10-10"
sys.path.insert(0,str(ROOT/"backend"))

from interface_engine import analyze_all_maps
from prepare import prepare_pdb_text
from geometry import extract_ca_nodes

CHANNEL_MAP={
    "anionic":"anionic",
    "H-bond acceptor":"hbond_acceptor",
    "cationic":"cationic",
    "H-bond donor":"hbond_donor",
    "hydrophobic":"hydrophobic",
    "pi/aromatic":"pi_carbon",
    "metal-oxide":"oxide",
    "calcium/phosphate charged sites":"calcium_phosphate",
    "transition-metal coordination":"metal_coord",
    "soft-metal sulfur affinity":"soft_metal_sulfur",
    "phosphate-rich":"phosphate",
}

def fetch_pdb(pid):
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.pdb",timeout=60) as r:
        return r.read().decode("utf-8",errors="replace")

def load_registry_channels():
    with (HERE/"candidate_registry.csv").open(encoding="utf-8-sig",newline="") as f:
        rows=list(csv.DictReader(f))
    out={}
    for r in rows:
        if r.get("status")!="ADMIT":
            continue
        channels=[x.strip() for x in str(r.get("chemistry_map","")).split(";") if x.strip()]
        unknown=[x for x in channels if x not in CHANNEL_MAP]
        if unknown:
            raise RuntimeError(f"{r['candidate_id']}: unknown frozen chemistry channel(s): {unknown}")
        out[r["candidate_id"]]=channels
    return out

def pool_patches(result,channels):
    pooled=[];seen=set()
    for ch in channels:
        pk=CHANNEL_MAP[ch]
        mm=result["maps"].get(pk,{})
        for p in mm.get("patches",[]):
            key=(str(p.get("center_key")),tuple(sorted(str(x) for x in p.get("members",[]))))
            if key in seen: continue
            seen.add(key)
            q=dict(p);q["source_channel"]=pk;pooled.append(q)
    pooled.sort(key=lambda p:(
        int(p.get("pareto_front",999)),
        -float(p.get("patch_coherence",0.0)),
        -float(p.get("chemistry_support",0.0)),
        str(p.get("source_channel","")),
        str(p.get("center_key",""))
    ))
    for i,p in enumerate(pooled,1):
        p["pooled_rank"]=i
    return pooled

def combined_propensity(result,channels):
    out={}
    for ch in channels:
        pk=CHANNEL_MAP[ch]
        for r in result["maps"].get(pk,{}).get("residues",[]):
            k=str(r.get("key"))
            if k:
                out[k]=max(out.get(k,0.0),float(r.get("propensity") or 0.0))
    return out

def coords(nodes):
    return {n["key"]:[float(x) for x in n["coord"]] for n in nodes}

def truth_keys(nodes,chain,exact):
    es={int(x) for x in exact}
    return {n["key"] for n in nodes if n["chain"]==chain and int(n["res_seq"]) in es}

def min_dist(members,truth,coord):
    best=math.inf
    for a in members:
        if a not in coord: continue
        for b in truth:
            if b in coord:
                best=min(best,math.dist(coord[a],coord[b]))
    return best

def first_rank(patches,truth,coord,cutoff=None):
    for i,p in enumerate(patches,1):
        m=set(p.get("members",[]))
        if cutoff is None:
            if m & truth:return i
        elif min_dist(m,truth,coord)<=cutoff:
            return i
    return None

def recall_near(patches,truth,coord,topk,cutoff):
    covered=set()
    for p in patches[:topk]:
        m=p.get("members",[])
        for t in truth:
            if t in coord and any(x in coord and math.dist(coord[x],coord[t])<=cutoff for x in m):
                covered.add(t)
    return len(covered)/len(truth) if truth else None

def recall_direct(patches,truth,topk):
    covered=set()
    for p in patches[:topk]:
        covered |= set(p.get("members",[])) & truth
    return len(covered)/len(truth) if truth else None

def best_front1(patches,truth,coord):
    vals=[]
    for p in patches:
        if int(p.get("pareto_front",999))==1:
            vals.append(recall_near([p],truth,coord,1,8.0))
    return max(vals) if vals else 0.0

def main():
    channel_by_id=load_registry_channels()
    with (HERE/"prospective_level3_ground_truth.csv").open(encoding="utf-8-sig",newline="") as f:
        rows=list(csv.DictReader(f))

    out=[]
    for r in rows:
        cid=r["candidate_id"]
        if cid not in channel_by_id:
            raise RuntimeError(f"{cid}: ground-truth row is not frozen ADMIT in candidate registry")
        channels=channel_by_id[cid]
        raw=fetch_pdb(r["pdb_id"])
        prepared,_=prepare_pdb_text(raw,chain=r["chain"])
        nodes=extract_ca_nodes(prepared)
        result=analyze_all_maps(
            pdb_text=raw,chain=r["chain"],pH=float(r["pH"]),
            ionic_mM=150.0,temp_K=298.0
        )
        patches=pool_patches(result,channels)
        prop=combined_propensity(result,channels)
        exact=[int(x) for x in str(r.get("exact_residues","")).split(";") if x.strip()]
        truth=truth_keys(nodes,r["chain"],exact)
        co=coords(nodes)

        requested=set(exact)
        mapped_res={int(n["res_seq"]) for n in nodes if n["chain"]==r["chain"] and int(n["res_seq"]) in requested}
        missing=sorted(requested-mapped_res)

        fd=first_rank(patches,truth,co,None)
        f5=first_rank(patches,truth,co,5.0)
        f8=first_rank(patches,truth,co,8.0)

        z=dict(r)
        z.update({
            "channels":";".join(channels),
            "n_truth_residues_requested":len(exact),
            "n_truth_residues_mapped":len(truth),
            "missing_truth_residue_numbers":";".join(map(str,missing)),
            "n_pooled_patches":len(patches),
            "first_direct_overlap_rank":fd or "",
            "first_near_5A_rank":f5 or "",
            "first_near_8A_rank":f8 or "",
            "top3_near_8A_hit":int((f8 or 10**9)<=3),
            "top5_near_8A_hit":int((f8 or 10**9)<=5),
            "top5_direct_overlap_hit":int((fd or 10**9)<=5),
            "top3_near_8A_recall":recall_near(patches,truth,co,3,8.0),
            "top5_near_8A_recall":recall_near(patches,truth,co,5,8.0),
            "top5_direct_overlap_recall":recall_direct(patches,truth,5),
            "best_front1_near_8A_recall":best_front1(patches,truth,co),
            "mean_truth_propensity":sum(prop.get(k,0.0) for k in truth)/len(truth) if truth else "",
            "max_truth_propensity":max([prop.get(k,0.0) for k in truth],default=""),
            "scoring_status":"scored" if truth and not missing else ("partial_truth_mapping" if truth else "non_evaluable")
        })
        out.append(z)

    path=HERE/"prospective_level3_results.csv"
    fields=[]
    for r in out:
        for k in r:
            if k not in fields:fields.append(k)
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader();w.writerows(out)

    scored=[x for x in out if x["scoring_status"]=="scored"]
    summary={
      "n_conditions":len(out),
      "n_scored":len(scored),
      "n_partial_truth_mapping":sum(x["scoring_status"]=="partial_truth_mapping" for x in out),
      "top3_near8_hits":sum(int(x["top3_near_8A_hit"]) for x in scored),
      "top5_near8_hits":sum(int(x["top5_near_8A_hit"]) for x in scored),
      "top5_direct_hits":sum(int(x["top5_direct_overlap_hit"]) for x in scored),
      "condition_results":{
          x["candidate_id"]:{
              "first_direct_overlap_rank":x["first_direct_overlap_rank"],
              "first_near_5A_rank":x["first_near_5A_rank"],
              "first_near_8A_rank":x["first_near_8A_rank"],
              "top3_near_8A_hit":x["top3_near_8A_hit"],
              "top5_near_8A_hit":x["top5_near_8A_hit"],
              "top5_direct_overlap_hit":x["top5_direct_overlap_hit"],
              "top3_near_8A_recall":x["top3_near_8A_recall"],
              "top5_near_8A_recall":x["top5_near_8A_recall"],
              "top5_direct_overlap_recall":x["top5_direct_overlap_recall"],
              "best_front1_near_8A_recall":x["best_front1_near_8A_recall"],
              "mean_truth_propensity":x["mean_truth_propensity"],
              "max_truth_propensity":x["max_truth_propensity"],
              "channels":x["channels"],
          } for x in scored
      }
    }
    (HERE/"prospective_level3_summary.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    print(json.dumps(out,indent=2))

if __name__=="__main__":
    main()
