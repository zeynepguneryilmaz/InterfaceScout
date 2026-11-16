import csv,io,json,os,platform,random,resource,statistics,subprocess,sys,time,zipfile
from pathlib import Path
R=Path(__file__).resolve().parents[1];O=R/"validation/performance_results";O.mkdir(parents=True,exist_ok=True)
Z=R/"validation/input_300_proteins.zip";SEED=20261117;BINS=[(60,119),(120,199),(200,349),(350,499),(500,800)]
def save(path,rows):
 with open(path,"w",newline="") as f:
  w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
def worker(inp,chain,out):
 sys.path.insert(0,str(R/"backend"))
 from interface_engine import analyze_all_maps
 pdb=Path(inp).read_text()
 start=time.perf_counter()
 res=analyze_all_maps(pdb_text=pdb,chain=chain,pH=7.4,ionic_mM=150,temp_K=298)
 d={"wall_s":time.perf_counter()-start,"peak_rss_mib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
    "n_maps":res["n_maps"],"n_patches":sum(m["n_patches"] for m in res["maps"].values())}
 Path(out).write_text(json.dumps(d))
def main():
 with zipfile.ZipFile(Z) as z:
  rows=list(csv.DictReader(io.StringIO(z.read("manifest.csv").decode("utf-8-sig"))))
  rng=random.Random(SEED);selected=[]
  for lo,hi in BINS:
   eligible=sorted([r for r in rows if r["pdb_id"]!="3IYK" and lo<=int(r["n_residues"])<=hi],key=lambda x:x["pdb_id"])
   selected.extend(rng.sample(eligible,8))
  rng.shuffle(selected);(O/"input").mkdir(exist_ok=True)
  for r in selected:(O/"input"/(r["pdb_id"]+"_"+r["chain"]+".pdb")).write_bytes(z.read("pdb/"+r["pdb_id"]+"_"+r["chain"]+".pdb"))
 manifest=[{"pdb_id":r["pdb_id"],"chain":r["chain"],"n_residues":int(r["n_residues"]),"original_sha256":r["sha256"]} for r in selected]
 save(O/"sample_manifest.csv",manifest)
 env={"platform":platform.platform(),"python":sys.version,"cpu_count_visible":os.cpu_count(),"protein_count":40,"replicates":3,
 "runner_cpu_model":next((s.split(":",1)[1].strip() for s in Path("/proc/cpuinfo").read_text().splitlines() if s.startswith("model name")),"not detected"),
 "time_definition":"wall time around complete analyze_all_maps; excludes Python startup and input file read; includes preparation,SASA,11 maps,patches and Pareto",
 "memory_definition":"peak resident memory (MiB) for isolated Python process including interpreter and imports, from Linux ru_maxrss",
 "execution":"sequential single process per protein/replicate; no parallel benchmark tasks","pH":7.4,"ionic_mM":150,"temperature_K":298.0,"seed":SEED}
 (O/"hardware_protocol.json").write_text(json.dumps(env,indent=2))
 raw=[]
 for i,r in enumerate(selected,1):
  inp=O/"input"/(r["pdb_id"]+"_"+r["chain"]+".pdb")
  for rep in range(1,4):
   out=O/"tmp.json"
   subprocess.run([sys.executable,__file__,"--worker",str(inp),r["chain"],str(out)],check=True,timeout=180,capture_output=True)
   d=json.loads(out.read_text());assert d["n_maps"]==11
   raw.append({"pdb_id":r["pdb_id"],"chain":r["chain"],"n_residues":int(r["n_residues"]),"replicate":rep,**d})
  save(O/"raw_repeats.csv",raw)
  print("PERF_OK",i,"/40",r["pdb_id"],flush=True)
 summary=[]
 for r in selected:
  rr=[q for q in raw if q["pdb_id"]==r["pdb_id"] and q["chain"]==r["chain"]]
  summary.append({"pdb_id":r["pdb_id"],"chain":r["chain"],"n_residues":int(r["n_residues"]),
   "median_wall_s":statistics.median(q["wall_s"] for q in rr),
   "median_peak_rss_mib":statistics.median(q["peak_rss_mib"] for q in rr),
   "min_wall_s":min(q["wall_s"] for q in rr),"max_wall_s":max(q["wall_s"] for q in rr)})
 save(O/"protein_performance.csv",summary)
 import numpy as np
 from scipy.stats import spearmanr,linregress
 lengths=np.array([r["n_residues"] for r in summary],float)
 results={}
 for name,vals in [("wall_time",np.array([r["median_wall_s"] for r in summary])),("peak_rss",np.array([r["median_peak_rss_mib"] for r in summary]))]:
  rho=spearmanr(lengths,vals);reg=linregress(np.log(lengths),np.log(vals))
  rng=np.random.default_rng(SEED);boots=[];brho=[]
  for _ in range(2500):
   ids=rng.integers(0,40,40)
   if len(set(lengths[ids]))<3:continue
   boots.append(linregress(np.log(lengths[ids]),np.log(vals[ids])).slope)
   brho.append(spearmanr(lengths[ids],vals[ids]).statistic)
  results[name]={"median":float(np.median(vals)),"q25":float(np.quantile(vals,.25)),"q75":float(np.quantile(vals,.75)),
   "min":float(vals.min()),"max":float(vals.max()),"spearman_rho":float(rho.statistic),"spearman_p":float(rho.pvalue),
   "loglog_slope":float(reg.slope),"loglog_r_squared":float(reg.rvalue**2),
   "slope_bootstrap_CI95":np.quantile(boots,[.025,.975]).tolist(),"rho_bootstrap_CI95":np.quantile(brho,[.025,.975]).tolist()}
 results["n_proteins"]=40;results["n_runs"]=len(raw)
 (O/"statistics.json").write_text(json.dumps(results,indent=2))
 import matplotlib
 matplotlib.use("Agg")
 import matplotlib.pyplot as plt
 fig,axs=plt.subplots(1,2,figsize=(10,4.2),layout="constrained")
 for ax,vals,label in [(axs[0],np.array([r["median_wall_s"] for r in summary]),"Analysis time (s)"),(axs[1],np.array([r["median_peak_rss_mib"] for r in summary]),"Peak RSS (MiB)")]:
  ax.scatter(lengths,vals,alpha=.8,s=32);ax.set_xlabel("Protein length (amino acids)");ax.set_ylabel(label);ax.grid(alpha=.2)
 fig.savefig(O/"Figure_Performance_Benchmark.png",dpi=350);fig.savefig(O/"Figure_Performance_Benchmark.svg")
 print("PERFORMANCE_COMPLETE",json.dumps(results),flush=True)
if __name__=="__main__":
 if len(sys.argv)>1 and sys.argv[1]=="--worker":worker(*sys.argv[2:5])
 else:main()
