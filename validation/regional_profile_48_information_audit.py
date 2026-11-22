"""Prospective descriptive utility audit of Regional Chemistry Profile on frozen 48.
Tested question: does a computed descriptor distinguish candidate regions beyond
binary co-location and single-channel support? NOT adsorption prediction.
The 48-PDB panel is selected independently of material experimental labels.
Run --shard 0..7, --shards 8 on pinned GitHub Actions runners.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,random,statistics,sys
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"backend"))
sys.path.insert(0,str(ROOT/"validation"))
from interface_engine import analyze_all_maps
from research_workflows import combine_maps
from strict_apo_compat_48 import prepare
PANEL=ROOT/"validation"/"strict_apo_48_selected_candidates_v2_final.csv"
IDS=[r["pdb_id"] for r in csv.DictReader(PANEL.open())]
assert len(IDS)==48 and len(set(IDS))==48
# Fixed interpretable channel pairs chosen before examining protein outcomes.
PAIRS=[
 ("anionic","hbond_donor"),
 ("cationic","hbond_acceptor"),
 ("hydrophobic","pi_carbon"),
 ("hydrophobic","hbond_acceptor"),
 ("oxide","hbond_donor"),
 ("calcium_phosphate","phosphate"),
 ("metal_coord","soft_metal_sulfur"),
 ("anionic","phosphate"),
]
OUT=Path("regional48"); OUT.mkdir(exist_ok=True)
def asfloat(x):return float(x) if x is not None else None
def ranked_corr(xs,ys):
 if len(xs)<4 or len(set(xs))<2 or len(set(ys))<2:return None
 c=float(spearmanr(xs,ys).statistic)
 return c if math.isfinite(c) else None
def pairs_discordant(xs,ys):
 # Number of region pairs ordered in opposite directions by metrics. Ties ignored.
 n=opp=0
 for i in range(len(xs)):
  for j in range(i+1,len(xs)):
   dx=xs[i]-xs[j]; dy=ys[i]-ys[j]
   if abs(dx)<1e-9 or abs(dy)<1e-9:continue
   n+=1;opp+=(dx*dy<0)
 return (opp/n if n else None,n)
def summary(pid,comb,a):
 out=combine_maps(a,comb,max_regions=None)
 rs=out["regions"]
 prof=[r["regional_profile"] for r in rs if r["all_maps_colocated"]]
 near=[r["min_distinct_intermap_CA_distance_A"] for r in prof
       if r["min_distinct_intermap_CA_distance_A"] is not None]
 s={ch:[p["mean_local_support_by_map"][ch] for p in prof] for ch in comb}
 f={ch:[p["supporting_residue_fraction_by_map"][ch] for p in prof] for ch in comb}
 rho={ch:ranked_corr(s[ch],f[ch]) for ch in comb}
 flip={ch:pairs_discordant(s[ch],f[ch]) for ch in comb}
 # Contrast equals mean_local_support minus constant per channel.
 contrast_same={}
 for ch in comb:
  c=[p["relative_to_candidate_median_by_map"][ch] for p in prof]
  contrast_same[ch]= (max(abs((x-y)-(s[ch][i]-s[ch][0])) for i,(x,y) in enumerate(zip(c,[c[0]]*len(c)))) < 0.000004) if c else None
 # Nearby residue pairs: co-located by common patch but pair may be far apart.
 # C-alpha 8A is the existing neighborhood radius, not a newly fitted threshold.
 raw={
  "pdb_id":pid,"maps":"+".join(comb),"n_regions":len(rs),"n_all_map_regions":len(prof),
  "allmap_fraction":len(prof)/len(rs) if rs else None,
  "exact_overlap_count":out["n_exact_overlap"],
  "intermap_distance_available_fraction":len(near)/len(prof) if prof else None,
  "fraction_pairs_beyond_8A":sum(v>8 for v in near)/len(near) if near else None,
  "min_intermap_CA_A":min(near) if near else None,
  "median_intermap_CA_A":statistics.median(near) if near else None,
  "max_intermap_CA_A":max(near) if near else None,
  "n_distinct_intermap_distance_values":len(set(near)),
  "support_fraction_spearman_m1":rho[comb[0]],"support_fraction_spearman_m2":rho[comb[1]],
  "discordant_order_fraction_m1":flip[comb[0]][0],"discordant_order_fraction_m2":flip[comb[1]][0],
  "n_discordance_comparisons_m1":flip[comb[0]][1],"n_discordance_comparisons_m2":flip[comb[1]][1],
  "region_variability_support_m1":len(set(s[comb[0]])),
  "region_variability_support_m2":len(set(s[comb[1]])),
  "region_variability_fraction_m1":len(set(f[comb[0]])),
  "region_variability_fraction_m2":len(set(f[comb[1]])),
  "contrast_affine_check_m1":contrast_same[comb[0]],
  "contrast_affine_check_m2":contrast_same[comb[1]],
 }
 per=[]
 for r in rs:
  p=r["regional_profile"]
  per.append({"pdb_id":pid,"maps":"+".join(comb),
   "center_key":r["source_center_key"],"source_map":r["source_map"],
   "member_count":len(r["member_keys"]),
   "all_maps_colocated":r["all_maps_colocated"],"n_represented_maps":r["n_represented_maps"],
   "mean_support_1":p["mean_local_support_by_map"][comb[0]],
   "mean_support_2":p["mean_local_support_by_map"][comb[1]],
   "fraction_1":p["supporting_residue_fraction_by_map"][comb[0]],
   "fraction_2":p["supporting_residue_fraction_by_map"][comb[1]],
   "nearest_distinct_intermap_CA_A":p["min_distinct_intermap_CA_distance_A"],
   "contrast_1":p["relative_to_candidate_median_by_map"][comb[0]],
   "contrast_2":p["relative_to_candidate_median_by_map"][comb[1]]})
 return raw,per
def write_csv(path,rows):
 if not rows:return
 keys=list(rows[0])
 with path.open("w",newline="",encoding="utf8") as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--shard",type=int,default=0)
 ap.add_argument("--shards",type=int,default=8);ap.add_argument("--aggregate",action="store_true")
 args=ap.parse_args()
 if args.aggregate:
  import glob
  files=sorted(glob.glob("regional48/regional48-shard-*/summary.csv"))
  assert len(files)==args.shards,("Incomplete shards",files)
  rows=[r for p in files for r in csv.DictReader(open(p))]
  ids={r["pdb_id"] for r in rows}
  assert ids==set(IDS),(len(ids),len(IDS),set(IDS)-ids)
  assert len(rows)==len(IDS)*len(PAIRS)
  write_csv(OUT/"all_48_summary.csv",rows)
  alld={}
  for col in ["allmap_fraction","intermap_distance_available_fraction",
    "fraction_pairs_beyond_8A","median_intermap_CA_A",
    "support_fraction_spearman_m1","support_fraction_spearman_m2",
    "discordant_order_fraction_m1","discordant_order_fraction_m2"]:
   x=[float(r[col]) for r in rows if r[col] not in ("",None)]
   alld[col]={"n":len(x),"mean":float(np.mean(x)) if x else None,"median":float(np.median(x)) if x else None,
      "q25":float(np.percentile(x,25)) if x else None,"q75":float(np.percentile(x,75)) if x else None}
  tests={
   "n_proteins":len(ids),"n_prespecified_map_pairs":len(PAIRS),
   "n_pair_evaluations":len(rows),"frozen_input":PANEL.name,
   "all_map_pairs":PAIRS,"stats":alld,
   "regions_with_more_than_one_distinct_proximity_value":
     sum(int(r["n_distinct_intermap_distance_values"])>1 for r in rows),
   "n_contrast_affine_checks_pass":
     sum(r["contrast_affine_check_m1"]=="True" and r["contrast_affine_check_m2"]=="True" for r in rows),
   "n_intermap_distance_unavailable":sum(not r["min_intermap_CA_A"] for r in rows),
   "no_experimental_labels":True,
   "interpretation":"All fields are deterministic descriptions of existing maps and candidate patches; no independent-information or adsorption prediction claim.",
   "caveats":["48 single structures do not establish across-conformation structural robustness",
     "Spearman between mean support and supporting fraction evaluates summary redundancy, not statistical independence",
     "8A limit is current patch radius; distinct pair distance can exceed 8A despite co-location within center neighborhood",
     "All feature comparisons use same-protein overlapping regions (not independent observations)"]
  }
  (OUT/"aggregate_summary.json").write_text(json.dumps(tests,indent=2))
  print("AGGREGATE",json.dumps(tests),flush=True);return
 shard=OUT/f"shard_{args.shard}";shard.mkdir(parents=True,exist_ok=True)
 rows=[];details=[];errors=[]
 for pid in IDS[args.shard::args.shards]:
  try:
   meta=prepare(pid)
   source=Path("strict_apo_panel_48")/f"{pid}_prepared.pdb"
   source=source.read_text()
   a=analyze_all_maps(pdb_text=source,chain=None,pH=7.4,ionic_mM=150,temp_K=298)
   assert a["n_maps"]==11
   for channels in PAIRS:
    r,d=summary(pid,channels,a);rows.append(r);details.extend(d)
   print("OK",pid,"regions",len(details),"maps",len(a["maps"]),flush=True)
  except Exception as exc:
   errors.append({"pdb":pid,"type":type(exc).__name__,"error":str(exc)})
   print("FAIL",pid,type(exc).__name__,str(exc),flush=True)
 write_csv(shard/"summary.csv",rows)
 write_csv(shard/"region_level.csv",details)
 (shard/"errors.json").write_text(json.dumps(errors,indent=2))
 (shard/"provenance.json").write_text(json.dumps({
   "shard":args.shard,"shards":args.shards,"selected_proteins":IDS[args.shard::args.shards],
   "pairs":PAIRS,"source_repo":"InterfaceScout","pH":7.4,"ionic_mM":150,
   "selection":"strict_apo_48_selected_candidates_v2_final.csv independent of adsorption labels"},indent=2))
 if errors or len(rows)!=len(IDS[args.shard::args.shards])*len(PAIRS):raise RuntimeError("Partial cohort; analysis refuses to claim 48 complete")
if __name__=="__main__":main()
