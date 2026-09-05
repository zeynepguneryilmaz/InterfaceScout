from pathlib import Path
import sys, urllib.request, json
import numpy as np, pandas as pd
from Bio.PDB import PDBParser,PDBIO,Select
from Bio.PDB.Polypeptide import is_aa
sys.path.insert(0,str(Path(__file__).resolve().parents[1])); import backend.main as ism
OUT=Path('validation/prospective_holdout/results_gb3_experimental_v1'); OUT.mkdir(parents=True,exist_ok=True)
PDBDIR=OUT/'pdb_cache'; PDBDIR.mkdir(exist_ok=True)
PDB='2OED'; CHAIN='A'; PH=7.0; CHEM='anionic'; TRUTH=[4,13,50]
class S(Select):
 def accept_model(self,m): return int(m.id==0)
 def accept_chain(self,c): return int(c.id==CHAIN)
 def accept_residue(self,r): return int(is_aa(r,standard=True))
raw=PDBDIR/f'{PDB}_raw.pdb'; clean=PDBDIR/f'{PDB}_curated.pdb'
if not raw.exists(): urllib.request.urlretrieve(f'https://files.rcsb.org/download/{PDB}.pdb',raw)
s=PDBParser(QUIET=True).get_structure(PDB,str(raw)); io=PDBIO(); io.set_structure(s); io.save(str(clean),S())
def auc(y,s):
 p=s[y==1]; n=s[y==0]; return float(np.mean([(a>b)+.5*(a==b) for a in p for b in n])) if len(p) and len(n) else np.nan
def ap(y,s):
 if y.sum()==0:return np.nan
 o=np.argsort(-s,kind='stable'); yy=y[o]; pr=np.cumsum(yy)/(np.arange(len(yy))+1); return float((pr*yy).sum()/yy.sum())
def rec(y,s,k): return float(y[np.argsort(-s,kind='stable')[:min(k,len(s))]].sum()/y.sum()) if y.sum() else np.nan
old=ism.SASA_POINTS; ism.SASA_POINTS=200
try:
 _,res,atoms,_=ism.build_surface_residues(clean,PH)
finally: ism.SASA_POINTS=old
truth_diag=pd.DataFrame([{k:r.get(k) for k in ['key','res_name','res_seq','sidechain_sasa','scrsa_raw','scrsa','surface_exposed']} for r in res if int(r['res_seq']) in TRUTH])
truth_diag.to_csv(OUT/'ground_truth_exposure_diagnostics.csv',index=False)
surf=[r for r in res if r['surface_exposed']]; D=ism.build_distances(surf); mp=ism.chemistry_map(surf,D,CHEM,PH)
core={x['center_key']:x['multiscale_persistence']/100 for x in mp['patch_centers']}
keys=[r['key'] for r in surf]; seq=np.array([int(r['res_seq']) for r in surf]); score=np.array([core.get(k,0.) for k in keys]); y=np.isin(seq,TRUTH).astype(int)
top=np.argsort(-score,kind='stable')[:10]; truth_idx=np.where(y==1)[0]; nearest=[float(np.min(D[ti,top])) for ti in truth_idx] if len(truth_idx) else []
row={'candidate_id':'GB3_CIT_AUNP_MUT','protein':'GB3','pdb':PDB,'pH':PH,'chemistry':CHEM,'n_total_residues':len(res),'n_surface':len(surf),'n_truth_requested':len(TRUTH),'n_truth_present_structure':len(truth_diag),'n_truth_mapped_surface':int(y.sum()),'AUROC':auc(y,score),'AP':ap(y,score),'Recall@5':rec(y,score,5),'Recall@10':rec(y,score,10),'Recall@20':rec(y,score,20),'Spatial@5A_top10':float(np.mean(np.array(nearest)<=5)) if nearest else np.nan,'Spatial@8A_top10':float(np.mean(np.array(nearest)<=8)) if nearest else np.nan,'Spatial@10A_top10':float(np.mean(np.array(nearest)<=10)) if nearest else np.nan,'median_nearest_top10_A':float(np.median(nearest)) if nearest else np.nan,'top10':';'.join(keys[i] for i in top)}
pd.DataFrame([row]).to_csv(OUT/'metrics.csv',index=False)
order=np.argsort(-score,kind='stable'); ranks=np.empty(len(order),int); ranks[order]=np.arange(1,len(order)+1)
pd.DataFrame([{'key':k,'res_seq':int(r),'score':float(sc),'rank':int(ranks[i]),'ground_truth':int(yy)} for i,(k,r,sc,yy) in enumerate(zip(keys,seq,score,y))]).to_csv(OUT/'residue_scores.csv',index=False)
(OUT/'METADATA.json').write_text(json.dumps({'model':'IS-v1-Core-primary','ground_truth_frozen_before_scoring':True,'SASA_POINTS':200,'source':'J Phys Chem C 2016 doi:10.1021/acs.jpcc.6b08469; K4/K13/K50 Ala variants significantly reduce GB3 adsorption to citrate-AuNP','note':'No model or GT changes after scoring. Diagnostics report whether experimental residues are lost at the frozen scRSA exposure gate.'},indent=2))
print(pd.DataFrame([row]).to_string(index=False)); print(truth_diag.to_string(index=False))
