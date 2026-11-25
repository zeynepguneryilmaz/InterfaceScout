"""Paired solver wall time and peak process RSS on original PQR/grid; not whole application timing."""
import csv,json,hashlib,os,sys,subprocess,shutil,time
from pathlib import Path
ref,out=map(Path,sys.argv[1:]);out.mkdir(exist_ok=True)
out=out.resolve()
rows=[]
for s in sorted(ref.rglob("summary.json")):
 d=json.loads(s.read_text());wd=out/d["pdb"];wd.mkdir(exist_ok=True)
 src=s.parent;pqr=src/"structure.pqr"
 assert hashlib.sha256(pqr.read_bytes()).hexdigest()==d["pqr_sha256"]
 shutil.copyfile(pqr,wd/"structure.pqr")
 a=(src/"lpbe.in").read_text();b=(src/"npbe.in").read_text()
 assert a.replace("    lpbe\n","    npbe\n")==b
 # Reverse first solver across structures, sequential within each structure.
 modes=["lpbe","npbe"] if sum(map(ord,d["pdb"]))%2==0 else ["npbe","lpbe"]
 for order,mode in enumerate(modes):
  for f in wd.glob("potential*.dx"):f.unlink()
  shutil.copyfile(src/f"{mode}.in",wd/"potential.in")
  start=time.monotonic()
  with (wd/f"{mode}.log").open("w") as log:
   p=subprocess.run(["/usr/bin/time","-f","%e,%M,%x","-o",f"{mode}_time.csv","apbs","potential.in"],cwd=wd,stdout=log,stderr=subprocess.STDOUT,timeout=2400,env={**os.environ,"OMP_NUM_THREADS":"1","OPENBLAS_NUM_THREADS":"1"})
  secs,rss,code=(wd/f"{mode}_time.csv").read_text().strip().splitlines()[-1].split(",")
  dx=list(wd.glob("potential*.dx"))
  r=dict(pdb=d["pdb"],solver=mode,order=order,wall_s=float(secs),peak_rss_KiB=int(rss),exit_code=p.returncode,dx_count=len(dx),grid=str(d["grid"]),pqr_sha256=d["pqr_sha256"],python_wall_s=time.monotonic()-start)
  rows.append(r);print(json.dumps(r),flush=True)
  for f in dx:f.unlink()
  with (out/"performance.csv").open("w",newline="") as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
assert rows and all(r["exit_code"]==0 and r["dx_count"]==1 for r in rows)
