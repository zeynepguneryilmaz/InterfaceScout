from pathlib import Path
import sys, json, urllib.request
import numpy as np, pandas as pd
from Bio.PDB import PDBParser,PDBIO,Select
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.SASA import ShrakeRupley
sys.path.insert(0,str(Path(__file__).resolve().parents[1])); import backend.main as ism

OUT=Path('validation/prospective_holdout/results_surface_is_v0'); OUT.mkdir(parents=True,exist_ok=True)
PDBDIR=OUT/'pdb_cache'; PDBDIR.mkdir(exist_ok=True)
SYSTEMS=[
 dict(system='GB3_CIT_AUNP_MUT',protein='GB3',pdb='2OED',chain='A',pH=6.4,truth=[4,13,50],truth_type='residue'),
 dict(system='B2M_CIT_AUNP_WT',protein='beta-2-microglobulin',pdb='1JNJ',chain='A',pH=7.7,truth=[2,3,26,28,29,30,33,55,56,58,59],truth_type='residue'),
 dict(system='HSA_PAA_FE3O4_XLMS',protein='HSA',pdb='2VUF',chain='A',pH=7.0,truth=list(range(373,390))+list(range(403,411))+list(range(414,429)),truth_type='peptide')]

class ChainAA(Select):
 def __init__(self,c,strip_h=False): self.c=c; self.strip_h=strip_h
 def accept_model(self,m): return int(m.id==0)
 def accept_chain(self,c): return int(c.id==self.c)
 def accept_residue(self,r): return int(is_aa(r,standard=True))
 def accept_atom(self,a):
  if self.strip_h and str(a.element).upper() in {'H','D'}: return 0
  return 1

def ensure(pid,chain):
 raw=PDBDIR/f'{pid}_raw.pdb'; cur=PDBDIR/f'{pid}_curated.pdb'; hs=PDBDIR/f'{pid}_hstrip.pdb'
 if not raw.exists(): urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pid}.pdb',raw)
 s=PDBParser(QUIET=True).get_structure(pid,str(raw)); io=PDBIO(); io.set_structure(s); io.save(str(cur),ChainAA(chain,False)); io.save(str(hs),ChainAA(chain,True)); return cur,hs

def auc(y,s):
 if y.sum()==0 or y.sum()==len(y): return np.nan
 p=s[y==1]; n=s[y==0]; return float(np.mean([(a>b)+.5*(a==b) for a in p for b in n]))
def ap(y,s):
 if y.sum()==0:return np.nan
 o=np.argsort(-s,kind='stable'); yy=y[o]; pr=np.cumsum(yy)/(np.arange(len(yy))+1); return float((pr*yy).sum()/yy.sum())
def rec(y,s,k):
 if y.sum()==0:return np.nan
 return float(y[np.argsort(-s,kind='stable')[:min(k,len(s))]].sum()/y.sum())
def spatial(seq,coords,y,score):
 top=np.argsort(-score,kind='stable')[:min(10,len(score))]; ti=np.where(y==1)[0]
 if not len(ti): return (np.nan,)*4
 nearest=[float(np.min(np.linalg.norm(coords[top]-coords[i],axis=1))) for i in ti]
 a=np.array(nearest); return float(np.mean(a<=5)),float(np.mean(a<=8)),float(np.mean(a<=10)),float(np.median(a))

def residue_core(pdb,pH):
 old=ism.SASA_POINTS; ism.SASA_POINTS=200
 try: _,res,_,_=ism.build_surface_residues(pdb,pH); surf=[r for r in res if r['surface_exposed']]; D=ism.build_distances(surf); mp=ism.chemistry_map(surf,D,'anionic',pH)
 finally: ism.SASA_POINTS=old
 core={x['center_key']:x['multiscale_persistence']/100 for x in mp['patch_centers']}
 allrows=[]
 for r in res:
  allrows.append((r['res_seq'],r['res_name'],r['x'],r['y'],r['z'],core.get(r['key'],0.0),r['surface_exposed'],r['scrsa'],r['key']))
 return allrows

def atomic_surface(pdb,pH):
 s=PDBParser(QUIET=True).get_structure('p',str(pdb)); sr=ShrakeRupley(probe_radius=1.4,n_points=200); sr.compute(s,level='A')
 atoms=[]; residues={}
 for model in s:
  for ch in model:
   for r in ch:
    if not is_aa(r,standard=True) or 'CA' not in r: continue
    rn=r.get_resname().strip().upper(); seq=int(r.id[1]); key=f'{ch.id}:{seq}:{str(r.id[2]).strip()}'
    ca=np.array(r['CA'].coord,float); residues[key]={'seq':seq,'rn':rn,'ca':ca,'atoms':[]}
    for a in r.get_atoms():
     if str(a.element).upper() in {'H','D'}: continue
     sasa=float(getattr(a,'sasa',0.0));
     if sasa<=0: continue
     name=a.get_name().strip().upper(); q=0.0
     if rn=='LYS' and name=='NZ': q=sasa*ism.state_availability(rn,'protonated',pH)
     elif rn=='ARG' and name in {'NE','NH1','NH2'}: q=sasa*ism.state_availability(rn,'protonated',pH)
     elif rn=='HIS' and name in {'ND1','NE2'}: q=sasa*ism.state_availability(rn,'protonated',pH)
     elif rn=='SER' and name=='OG': q=sasa
     elif rn=='THR' and name=='OG1': q=sasa
     rec={'key':key,'seq':seq,'rn':rn,'atom':name,'coord':np.array(a.coord,float),'sasa':sasa,'q':q}; atoms.append(rec); residues[key]['atoms'].append(rec)
  break
 X=np.vstack([a['coord'] for a in atoms]); q=np.array([a['q'] for a in atoms]); D=np.linalg.norm(X[:,None,:]-X[None,:,:],axis=2)
 vals=[]
 for R in (5.0,8.0): vals.append((D<=R).astype(float)@q)
 norms=[]
 for v in vals: norms.append(v/v.max() if v.max()>0 else np.zeros_like(v))
 M=100*np.minimum(norms[0],norms[1])
 score_by={k:0.0 for k in residues}
 for i,a in enumerate(atoms): score_by[a['key']]=max(score_by[a['key']],float(M[i]/100.0))
 rows=[]; atomdiag=[]
 for k,r in residues.items():
  rows.append((r['seq'],r['rn'],*r['ca'],score_by[k],True,np.nan,k))
  for a in r['atoms']:
   if a['q']>0: atomdiag.append({'key':k,'res_seq':r['seq'],'res_name':r['rn'],'atom':a['atom'],'atom_sasa':a['sasa'],'chem_weighted_sasa':a['q']})
 return rows,atomdiag

metrics=[]; residue_rows=[]; atom_rows=[]
for sysd in SYSTEMS:
 cur,hs=ensure(sysd['pdb'],sysd['chain'])
 arms=[('RAW_RESIDUE_CORE',residue_core(cur,sysd['pH'])),('HSTRIP_RESIDUE_CORE',residue_core(hs,sysd['pH']))]
 ar,adiag=atomic_surface(hs,sysd['pH']); arms.append(('ATOMIC_SURFACE_CORE',ar))
 for x in adiag: x['system']=sysd['system']; atom_rows.append(x)
 for arm,rows in arms:
  seq=np.array([int(r[0]) for r in rows]); rn=[r[1] for r in rows]; coords=np.array([[r[2],r[3],r[4]] for r in rows],float); score=np.array([float(r[5]) for r in rows]); y=np.isin(seq,sysd['truth']).astype(int)
  sp5,sp8,sp10,med=spatial(seq,coords,y,score)
  order=np.argsort(-score,kind='stable'); top10=';'.join(f'A:{seq[i]}:' for i in order[:10])
  metrics.append({'system':sysd['system'],'protein':sysd['protein'],'arm':arm,'n_residues':len(seq),'truth_present':int(y.sum()),'AUROC':auc(y,score),'AP':ap(y,score),'Recall@5':rec(y,score,5),'Recall@10':rec(y,score,10),'Recall@20':rec(y,score,20),'Spatial@5A_top10':sp5,'Spatial@8A_top10':sp8,'Spatial@10A_top10':sp10,'median_nearest_top10_A':med,'top10':top10})
  ranks=np.empty(len(order),int); ranks[order]=np.arange(1,len(order)+1)
  for i in range(len(seq)): residue_rows.append({'system':sysd['system'],'arm':arm,'res_seq':int(seq[i]),'res_name':rn[i],'score':float(score[i]),'rank':int(ranks[i]),'ground_truth':int(y[i])})
pd.DataFrame(metrics).to_csv(OUT/'metrics.csv',index=False); pd.DataFrame(residue_rows).to_csv(OUT/'residue_scores.csv',index=False); pd.DataFrame(atom_rows).to_csv(OUT/'functional_group_surface.csv',index=False)
(OUT/'METADATA.json').write_text(json.dumps({'protocol':'SURFACE_IS_PROTOCOL_V0','model_status':'exploratory; frozen IS-v1 unchanged','sasa_points':200,'probe_A':1.4,'surface_arm':'H-stripped exposed heavy atoms; functional-group SASA; 5/8A persistence'},indent=2))
print(pd.DataFrame(metrics).to_string(index=False))
print('\nGB3 functional groups:')
print(pd.DataFrame(atom_rows).query("system=='GB3_CIT_AUNP_MUT' and res_seq in [4,13,50]").to_string(index=False))
