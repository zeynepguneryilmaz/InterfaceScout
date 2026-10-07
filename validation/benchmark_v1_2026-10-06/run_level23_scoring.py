from __future__ import annotations
import csv, json, math, os, sys, urllib.request
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "validation" / "benchmark_v1_2026-10-06"
sys.path.insert(0, str(ROOT / "backend"))

from interface_engine import analyze_all_maps
from prepare import prepare_pdb_text
from geometry import extract_ca_nodes

CHANNEL_MAP = {
    "anionic":"anionic",
    "cationic":"cationic",
    "hydrophobic":"hydrophobic",
    "pi/aromatic":"pi_carbon",
    "H-bond donor":"hbond_donor",
    "H-bond acceptor":"hbond_acceptor",
    "metal-oxide":"oxide",
    "calcium/phosphate charged sites":"calcium_phosphate",
    "transition-metal coordination":"metal_coord",
    "soft-metal sulfur affinity":"soft_metal_sulfur",
    "phosphate-rich":"phosphate",
}
LEVEL3_CHANNELS = {
    "hsa_paa_fe3o4":["anionic","hbond_acceptor"],
    "ubiquitin_citrate_au":["anionic","hbond_acceptor"],
    "hcaii_silica_ph63":["oxide"],
    "hcaii_silica_ph83":["oxide"],
    "b2m_mpa_au":["anionic","hbond_acceptor"],
    "asparaginase_alumoh":["cationic","oxide"],
    "bsa_silica":["oxide"],
}
LEVEL2_CHANNELS = {
    ("polystyrene","none"):["hydrophobic","pi_carbon"],
    ("iron oxide","none"):["oxide"],
}
CASE_CHAIN = {
    "hsa_paa_fe3o4":"A",
    "ubiquitin_citrate_au":"A",
    "hcaii_silica_ph63":"A",
    "hcaii_silica_ph83":"A",
    "b2m_mpa_au":"A",
    "asparaginase_alumoh":"A,B,C,D",
    "bsa_silica":"A",
    "transferrin_PS":"C,E",
    "apo_transferrin_PS":"C,E",
    "catalase_PS":"A,B,C,D",
    "transferrin_IO":"C,E",
}

def fetch_pdb(pdb_id:str)->str:
    url=f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read().decode("utf-8", errors="replace")

def pool_patches(result:dict, channels:list[str])->list[dict]:
    pooled=[]
    seen=set()
    for ch in channels:
        public=CHANNEL_MAP.get(ch,ch)
        m=result["maps"].get(public)
        if not m: continue
        for p in m.get("patches",[]):
            key=(str(p.get("center_key")), tuple(sorted(str(x) for x in p.get("members",[]))))
            if key in seen: continue
            seen.add(key)
            q=dict(p)
            q["source_channel"]=public
            pooled.append(q)
    pooled.sort(key=lambda p:(
        int(p.get("pareto_front",999)),
        -float(p.get("patch_coherence",0.0)),
        -float(p.get("chemistry_support",0.0)),
        str(p.get("source_channel","")),
        str(p.get("center_key","")),
    ))
    for i,p in enumerate(pooled,1):
        p["pooled_rank"]=i
    return pooled

def combined_propensity(result:dict, channels:list[str])->dict[str,float]:
    out={}
    for ch in channels:
        public=CHANNEL_MAP.get(ch,ch)
        m=result["maps"].get(public)
        if not m: continue
        for r in m.get("residues",[]):
            k=str(r.get("key"))
            if not k: continue
            v=float(r.get("propensity") or 0.0)
            out[k]=max(out.get(k,0.0),v)
    return out

def truth_keys(nodes:list[dict], chain_sel:str, start, end, exact)->set[str]:
    chains={x.strip() for x in (chain_sel or "").split(",") if x.strip()}
    exact_set=set(exact or [])
    keys=set()
    for n in nodes:
        if chains and n["chain"] not in chains: continue
        rs=int(n["res_seq"])
        if exact_set and rs in exact_set:
            keys.add(n["key"])
        elif start is not None and end is not None and int(start)<=rs<=int(end):
            keys.add(n["key"])
    return keys

def min_dist_to_truth(member_keys, truth_keys_set, coord):
    best=math.inf
    for a in member_keys:
        if a not in coord: continue
        ax=coord[a]
        for b in truth_keys_set:
            if b not in coord: continue
            bx=coord[b]
            d=math.dist(ax,bx)
            if d<best: best=d
    return best

def recall_near(patches, truth, coord, topk, cutoff):
    if not truth: return None
    covered=set()
    for p in patches[:topk]:
        members=p.get("members",[])
        for t in truth:
            if t not in coord: continue
            tc=coord[t]
            if any(m in coord and math.dist(coord[m],tc)<=cutoff for m in members):
                covered.add(t); continue
    return len(covered)/len(truth) if truth else None

def recall_direct(patches, truth, topk):
    if not truth: return None
    covered=set()
    for p in patches[:topk]:
        covered |= set(p.get("members",[])) & truth
    return len(covered)/len(truth)

def first_rank(patches, truth, coord, cutoff=None):
    for i,p in enumerate(patches,1):
        members=set(p.get("members",[]))
        if cutoff is None:
            if members & truth: return i
        else:
            if min_dist_to_truth(members, truth, coord)<=cutoff: return i
    return None

def best_front1_recall(patches, truth, coord, cutoff=8.0):
    if not truth: return None
    vals=[]
    for p in patches:
        if int(p.get("pareto_front",999))!=1: continue
        vals.append(recall_near([p], truth, coord, 1, cutoff))
    return max(vals) if vals else 0.0

def coords_from_nodes(nodes):
    return {n["key"]:[float(x) for x in n["coord"]] for n in nodes}

def parse_exact(text):
    if not text: return []
    out=[]
    for x in str(text).split(";"):
        x=x.strip()
        if not x: continue
        try: out.append(int(float(x)))
        except: pass
    return out

def parse_int(v):
    if v is None or v=="": return None
    try: return int(float(v))
    except: return None

def score_record(base:dict, pdb_id:str, chain:str, pH:float, channels:list[str], start, end, exact):
    raw=fetch_pdb(pdb_id)
    prepared,_=prepare_pdb_text(raw,chain=chain)
    nodes=extract_ca_nodes(prepared)
    result=analyze_all_maps(pdb_text=raw,chain=chain,pH=pH,ionic_mM=150.0,temp_K=298.0)
    patches=pool_patches(result,channels)
    prop=combined_propensity(result,channels)
    truth=truth_keys(nodes,chain,start,end,exact)
    coord=coords_from_nodes(nodes)
    row=dict(base)
    row.update({
        "channels":";".join(channels),
        "n_truth_residues":len(truth),
        "n_pooled_patches":len(patches),
        "first_direct_overlap_rank":first_rank(patches,truth,coord,None) if truth else "",
        "first_near_5A_rank":first_rank(patches,truth,coord,5.0) if truth else "",
        "first_near_8A_rank":first_rank(patches,truth,coord,8.0) if truth else "",
        "top3_near_8A_hit":int((first_rank(patches,truth,coord,8.0) or 10**9)<=3) if truth else "",
        "top5_near_8A_hit":int((first_rank(patches,truth,coord,8.0) or 10**9)<=5) if truth else "",
        "top5_direct_overlap_hit":int((first_rank(patches,truth,coord,None) or 10**9)<=5) if truth else "",
        "top3_near_8A_recall":recall_near(patches,truth,coord,3,8.0) if truth else "",
        "top5_near_8A_recall":recall_near(patches,truth,coord,5,8.0) if truth else "",
        "top5_direct_overlap_recall":recall_direct(patches,truth,5) if truth else "",
        "best_front1_near_8A_recall":best_front1_recall(patches,truth,coord,8.0) if truth else "",
        "mean_truth_propensity":sum(prop.get(k,0.0) for k in truth)/len(truth) if truth else "",
        "max_truth_propensity":max([prop.get(k,0.0) for k in truth], default="") if truth else "",
        "structure_origin":"experimental PDB",
        "scoring_status":"scored" if truth else "non_evaluable_non_numeric_or_unmapped_truth",
    })
    return row

def write_csv(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    fields=[]
    for r in rows:
        for k in r:
            if k not in fields: fields.append(k)
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader(); w.writerows(rows)

def main():
    outdir=BENCH/"scoring_v1_2026-10-07"
    outdir.mkdir(exist_ok=True)
    l3=[]
    with (BENCH/"residue_level_ground_truth.csv").open(encoding="utf-8-sig",newline="") as f:
        for r in csv.DictReader(f):
            sid=r["study_id"]
            channels=LEVEL3_CHANNELS[sid]
            exact=parse_exact(r.get("exact_residues"))
            start=parse_int(r.get("region_start")); end=parse_int(r.get("region_end"))
            base={
                "benchmark_level":"LEVEL3","source_id":r["source_id"],"study_id":sid,
                "source_row_identifier":r["source_row_identifier"],"protein_name":r["protein_name"],
                "pdb_id":r["pdb_id"],"chain":r["chain"],"pH":r["pH"],
                "ground_truth_type":r["ground_truth_type"],"evidence_resolution":r["evidence_resolution"],
                "confidence_quality_flag":r["confidence_quality_flag"],"DOI":r["DOI"],
            }
            chain=(r["chain"] or "").replace("/",",")
            base["chain"]=chain
            try:
                l3.append(score_record(base,r["pdb_id"],chain,float(r["pH"]),channels,start,end,exact))
            except Exception as e:
                z=dict(base); z.update({"channels":";".join(channels),"scoring_status":"failed","error":repr(e)}); l3.append(z)
    write_csv(outdir/"level3_spatial_results.csv",l3)

    l2=[]
    with (BENCH/"structural_response_ground_truth.csv").open(encoding="utf-8-sig",newline="") as f:
        for r in csv.DictReader(f):
            channels=LEVEL2_CHANNELS[(r["material"].strip().lower(),r["surface_modification"].strip().lower())]
            exact=[]
            start=parse_int(r.get("region_start")); end=parse_int(r.get("region_end"))
            chain=CASE_CHAIN[r["study_id"]]
            base={
                "benchmark_level":"LEVEL2","source_id":r["source_id"],"study_id":r["study_id"],
                "source_row_identifier":r["source_row_identifier"],"protein_name":r["protein_name"],
                "pdb_id":r["pdb_id"],"chain":chain,"pH":r["pH"],
                "response_label":r["response_label"],"experimental_direction":r["experimental_direction"],
                "contact_inference_allowed":r["contact_inference_allowed"],"DOI":r["DOI"],
            }
            if start is None or end is None:
                z=dict(base); z.update({"channels":";".join(channels),"scoring_status":"non_evaluable_non_numeric_region"})
                l2.append(z); continue
            try:
                l2.append(score_record(base,r["pdb_id"],chain,float(r["pH"]),channels,start,end,exact))
            except Exception as e:
                z=dict(base); z.update({"channels":";".join(channels),"scoring_status":"failed","error":repr(e)}); l2.append(z)
    write_csv(outdir/"level2_spatial_results.csv",l2)

    summary={
        "level3_records":len(l3),
        "level3_scored":sum(r.get("scoring_status")=="scored" for r in l3),
        "level3_failed":sum(r.get("scoring_status")=="failed" for r in l3),
        "level2_records":len(l2),
        "level2_scored":sum(r.get("scoring_status")=="scored" for r in l2),
        "level2_non_evaluable":sum(str(r.get("scoring_status","")).startswith("non_evaluable") for r in l2),
        "level2_failed":sum(r.get("scoring_status")=="failed" for r in l2),
    }
    (outdir/"level23_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(json.dumps(summary,indent=2))

if __name__=="__main__":
    main()
