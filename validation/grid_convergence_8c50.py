"""8C50 grid study on exact frozen PQR and archived shell-supported residue set."""
import csv,json,hashlib,math,re,sys,shutil,time
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'backend'))
from electrostatics import _pqr_atoms,_run,_read_dx
from surface_potential import sample_residue_sas_potential
mode,solver,reference,out=sys.argv[1:];ref=Path(reference).resolve();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
src=next(ref.rglob('8C50/structure.pqr'));src=src.parent
summary=json.loads((src/'summary.json').read_text());pqr=src/'structure.pqr'
assert hashlib.sha256(pqr.read_bytes()).hexdigest()==summary['pqr_sha256']
shutil.copyfile(pqr,out/'structure.pqr')
with (src/'paired.csv').open() as f: old={r['key']:r for r in csv.DictReader(f)}
assert len(old)==683
atoms=_pqr_atoms(pqr);xyz=np.asarray([a['xyz'] for a in atoms]);lo=xyz.min(0);hi=xyz.max(0);span=np.maximum(hi-lo,1)
center=(lo+hi)/2
# Vary grid spacing, box size and sub-grid translation separately. No PQR change.
configs={'bbox_h08':(.8,10,20,0.),'bbox_h06':(.6,10,20,0.),
 'bbox_h05':(.5,10,20,0.),'expanded_h06':(.6,20,40,0.),
 'shifted_h06':(.6,10,20,.3)}
inp=(src/f'{solver}.in').read_text()
if mode!='archived':
 h,pad,outer,shift=configs[mode];fine=span+2*pad;coarse=span+2*outer
 dims=[max(65,32*math.ceil(v/h/32)+1) for v in fine]
 assert np.prod(dims)<16_000_000, dims
 for key,val in [('dime',' '.join(map(str,dims))),('fglen',' '.join(f'{v:.4f}' for v in fine)),('cglen',' '.join(f'{v:.4f}' for v in coarse)),('fgcent',' '.join(f'{v:.4f}' for v in center+shift)),('cgcent',' '.join(f'{v:.4f}' for v in center+shift))]:
  inp,n=re.subn(r'(?m)^    '+key+r' .+$','    '+key+' '+val,inp);assert n==1,key
(out/'potential.in').write_text(inp)
t=time.monotonic()
try:
 log=_run([shutil.which('apbs'),'potential.in'],out,timeout=2400)
 (out/'apbs.log').write_text(log)
 dx=list(out.glob('potential*.dx'));assert len(dx)==1,dx
 g,o,s=_read_dx(dx[0])
 data,meta=sample_residue_sas_potential(atoms,g,o,s,set(old))
 assert set(data)==set(old)
 for k in old:assert int(old[k]['shell_samples'])==data[k]['surface_potential_samples'],k
 keys=sorted(old);v=np.array([data[k]['surface_potential_kT_per_e'] for k in keys])
 a=np.array([float(old[k][solver.upper()+'_kT_e']) for k in keys]);d=v-a
 warnings=[l.strip() for l in log.splitlines() if re.search(r'trapped exp overflows|exponential overflow|convergence not guaranteed|focusFillBound.*warning|failed to converge|did not converge|non.?finite',l,re.I)]
 with (out/'potentials.csv').open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['key','potential_kT_e','archived_kT_e','delta','samples','area_A2'])
  for k,x,y,z in zip(keys,v,a,d):w.writerow([k,x,y,z,data[k]['surface_potential_samples'],data[k]['surface_sampling_area_A2']])
 q={'status':'computed','config':mode,'solver':solver,'pqr_sha256':summary['pqr_sha256'],'grid_shape':list(g.shape),'origin_A':o.tolist(),'spacing_A':s.tolist(),'max_spacing_A':float(max(s)),'n_paired':len(v),'elapsed_seconds':time.monotonic()-t,'median_potential':float(np.median(v)),'min_potential':float(min(v)),'max_potential':float(max(v)),'n_abs_gt10':int(sum(abs(v)>10)),'median_abs_delta_vs_archived':float(np.median(abs(d))),'p95_abs_delta_vs_archived':float(np.quantile(abs(d),.95)),'max_abs_delta_vs_archived':float(max(abs(d))),'spearman_vs_archived':float(spearmanr(a,v).statistic),'sign_disagreement_vs_archived':float(np.mean(np.sign(a)!=np.sign(v))),'diagnostic_warnings':warnings,'surface_sampling':meta}
 (out/'summary.json').write_text(json.dumps(q,indent=2,allow_nan=False))
 print(json.dumps(q),flush=True)
 dx[0].unlink()
except Exception as e:
 (out/'failure.json').write_text(json.dumps({'config':mode,'solver':solver,'error':str(e)},indent=2));raise
