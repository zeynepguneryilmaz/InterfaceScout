from __future__ import annotations
import json, math, sys
from pathlib import Path
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, PDBIO

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import backend.main as ism

BASE = Path('validation/packages/all15/results_ready_batch')
PDBDIR = BASE / 'pdb_cache'
OUT = Path('validation/packages/all15/results_anm_exploratory')
OUT.mkdir(parents=True, exist_ok=True)
CONF_DIR = OUT / 'conformers'
CONF_DIR.mkdir(exist_ok=True)

# Same predeclared benchmark set as first batch, excluding negative control from primary summary.
BENCH = [
    {'id':'FN_COO','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'anionic','pH':7.0,'truth':[1469]},
    {'id':'FN_NH3','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'cationic','pH':7.0,'truth':[1312,1509]},
    {'id':'FN_CH3_HEAD','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'hydrophobic','pH':7.0,'truth':[1454,1455,1456,1457,1478,1479,1480,1481,1509]},
    {'id':'FN_CH3_SIDE','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'hydrophobic','pH':7.0,'truth':[1275,1276,1355,1376,1378,1454,1455,1456,1457]},
    {'id':'FN_CH3_BETA','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'hydrophobic','pH':7.0,'truth':[1431,1446,1457,1458,1459,1461,1463,1464,1466,1475,1476]},
    {'id':'CYTC_CH3','protein':'Cytochrome c','pdb':'3NWV','chem':'hydrophobic','pH':7.0,'truth':[1,2,3,4,96,99,100]},
    {'id':'CYTC_COOH','protein':'Cytochrome c','pdb':'3NWV','chem':'anionic','pH':7.0,'truth':[1,2,4,61,99,100,103]},
    {'id':'LYZ_SWCNT','protein':'Lysozyme','pdb':'1LYZ','chem':'pi_carbon','pH':7.4,'truth':[61,62,63,68,69,70,73,103,106,107,111,112,113,116]},
    {'id':'CHT_CNT_HYDRO','protein':'Alpha-chymotrypsin','pdb':'4CHA','chem':'hydrophobic','pH':7.0,'truth':[3,5,6,7,8,10,77,114,115,116]},
    {'id':'MB_CIT_AUNP','protein':'Myoglobin','pdb':'1MBN','chem':'anionic','pH':7.0,'truth':[43,45,96,97,98,99,146,147,148,149,150,151,152,153]},
]

N_MODES = 5
TARGET_CA_RMSD_A = 1.0
SASA_POINTS = 120  # exploratory ensemble pass; static comparison uses same setting for fairness
TOP_FRAC = 0.10


def residue_key(chain_id, res):
    het, seq, icode = res.id
    return f"{chain_id}:{int(seq)}:{str(icode).strip()}"


def load_structure(pid):
    path=PDBDIR/f'{pid}.pdb'
    if not path.exists(): raise FileNotFoundError(path)
    return PDBParser(QUIET=True).get_structure(pid,str(path))


def ca_nodes(structure):
    model=next(structure.get_models())
    nodes=[]
    for chain in model:
        for res in chain:
            if 'CA' in res:
                nodes.append((residue_key(chain.id,res),res,np.array(res['CA'].coord,float)))
    return nodes


def build_anm_modes(coords, cutoff=15.0, gamma=1.0, n_modes=5):
    n=len(coords); H=np.zeros((3*n,3*n),float)
    for i in range(n-1):
        for j in range(i+1,n):
            d=coords[j]-coords[i]; r=np.linalg.norm(d)
            if r<=0 or r>cutoff: continue
            K=-gamma*np.outer(d,d)/(r*r)
            si=slice(3*i,3*i+3); sj=slice(3*j,3*j+3)
            H[si,sj]+=K; H[sj,si]+=K
            H[si,si]-=K; H[sj,sj]-=K
    vals,vecs=np.linalg.eigh(H)
    pos=np.where(vals>1e-8)[0]
    if len(pos)<n_modes: raise RuntimeError('Too few nonzero ANM modes')
    idx=pos[:n_modes]
    return vals[idx],vecs[:,idx]


def write_conformer(pid, base_structure, nodes, mode_vec, sign, idx):
    # clone by reparsing clean PDB so every conformer starts from identical coordinates
    s=load_structure(pid); model=next(s.get_models())
    node_map={k:i for i,(k,_,_) in enumerate(nodes)}
    disp=mode_vec.reshape(len(nodes),3).copy()
    rms=np.sqrt(np.mean(np.sum(disp*disp,axis=1)))
    if rms==0: raise RuntimeError('zero mode RMS')
    disp*= (TARGET_CA_RMSD_A/rms) * sign
    for chain in model:
        for res in chain:
            k=residue_key(chain.id,res)
            if k not in node_map: continue
            dv=disp[node_map[k]]
            for atom in res:
                atom.coord=np.asarray(atom.coord,float)+dv
    path=CONF_DIR/f'{pid}_m{idx+1}_{"plus" if sign>0 else "minus"}.pdb'
    io=PDBIO(); io.set_structure(s); io.save(str(path))
    return path


def get_scores(pdb_path, chem, ph):
    old=ism.SASA_POINTS; ism.SASA_POINTS=SASA_POINTS
    try:
        _,allres,_,_=ism.build_surface_residues(pdb_path,ph)
    finally:
        ism.SASA_POINTS=old
    surf=[r for r in allres if r['surface_exposed']]
    D=ism.build_distances(surf)
    mp=ism.chemistry_map(surf,D,chem,ph)
    sm={x['center_key']:x['multiscale_persistence']/100.0 for x in mp['patch_centers']}
    return {r['key']:{'score':sm.get(r['key'],0.0),'scrsa':r['scrsa'],'coord':np.asarray([r['x'],r['y'],r['z']],float),'seq':int(r['res_seq'])} for r in surf}


def auc(y,s):
    p=s[y==1]; n=s[y==0]
    if not len(p) or not len(n): return np.nan
    return float(np.mean([(a>b)+0.5*(a==b) for a in p for b in n]))

def ap(y,s):
    if not y.sum(): return np.nan
    o=np.argsort(-s,kind='stable'); yy=y[o]; prec=np.cumsum(yy)/(np.arange(len(yy))+1)
    return float((prec*yy).sum()/yy.sum())

def recall(y,s,k):
    if not y.sum(): return np.nan
    o=np.argsort(-s,kind='stable')[:min(k,len(s))]
    return float(y[o].sum()/y.sum())

def eval_metric(ref, truth_seq, scores):
    keys=list(ref.keys()); s=np.asarray([scores.get(k,0.0) for k in keys],float)
    y=np.asarray([int(ref[k]['seq'] in truth_seq) for k in keys],int)
    return {'AUROC':auc(y,s),'AP':ap(y,s),'Recall@5':recall(y,s,5),'Recall@10':recall(y,s,10),'Recall@20':recall(y,s,20)}

# Build conformers once per protein.
protein_confs={}
for pid in sorted(set(b['pdb'] for b in BENCH)):
    s=load_structure(pid); nodes=ca_nodes(s); coords=np.vstack([x[2] for x in nodes])
    vals,vecs=build_anm_modes(coords,n_modes=N_MODES)
    paths=[]
    for m in range(N_MODES):
        for sign in (-1,1): paths.append(write_conformer(pid,s,nodes,vecs[:,m],sign,m))
    protein_confs[pid]={'paths':paths,'eigenvalues':vals.tolist(),'n_ca':len(nodes)}

rows=[]; residue_rows=[]
for b in BENCH:
    static=get_scores(PDBDIR/f"{b['pdb']}.pdb",b['chem'],b['pH'])
    keys=list(static.keys())
    conf_scores=[]
    for p in protein_confs[b['pdb']]['paths']:
        sc=get_scores(p,b['chem'],b['pH'])
        conf_scores.append(np.asarray([sc.get(k,{}).get('score',0.0) for k in keys],float))
    A=np.vstack(conf_scores)
    mean=A.mean(axis=0)
    # per-conformer top 10% persistence
    persistence=np.zeros(len(keys),float)
    for row in A:
        k=max(1,int(math.ceil(TOP_FRAC*len(row))))
        idx=np.argsort(-row,kind='stable')[:k]; persistence[idx]+=1
    persistence/=A.shape[0]
    consensus=mean*persistence
    static_arr=np.asarray([static[k]['score'] for k in keys],float)
    for method,arr in [('static_sameSASA',static_arr),('anm_mean',mean),('anm_persistence',persistence),('anm_consensus',consensus)]:
        met=eval_metric(static,b['truth'],dict(zip(keys,arr)))
        rows.append({'condition_id':b['id'],'protein':b['protein'],'pdb':b['pdb'],'method':method,**met})
    for i,k in enumerate(keys):
        residue_rows.append({'condition_id':b['id'],'protein':b['protein'],'key':k,'res_seq':static[k]['seq'],'truth':int(static[k]['seq'] in b['truth']),'static':static_arr[i],'anm_mean':mean[i],'anm_persistence':persistence[i],'anm_consensus':consensus[i]})

res=pd.DataFrame(rows)
res.to_csv(OUT/'anm_method_comparison.csv',index=False)
pd.DataFrame(residue_rows).to_csv(OUT/'anm_residue_scores.csv',index=False)

# Paired deltas versus same-SASA static baseline.
wide=res.pivot(index=['condition_id','protein','pdb'],columns='method',values=['AUROC','AP','Recall@5','Recall@10','Recall@20'])
flat=[]
for idx,row in wide.iterrows():
    d={'condition_id':idx[0],'protein':idx[1],'pdb':idx[2]}
    for metric in ['AUROC','AP','Recall@5','Recall@10','Recall@20']:
        base=row[(metric,'static_sameSASA')]
        for method in ['anm_mean','anm_persistence','anm_consensus']:
            d[f'{metric}_delta_{method}']=row[(metric,method)]-base
    flat.append(d)
pd.DataFrame(flat).to_csv(OUT/'anm_paired_deltas.csv',index=False)

summary={
    'design':{'n_modes':N_MODES,'n_conformers_per_protein':2*N_MODES,'target_ca_rmsd_A':TARGET_CA_RMSD_A,'sasa_points':SASA_POINTS,'top_fraction_for_persistence':TOP_FRAC},
    'proteins':protein_confs,
    'median_metrics_by_method':res.groupby('method')[['AUROC','AP','Recall@5','Recall@10','Recall@20']].median().to_dict(orient='index'),
    'mean_metrics_by_method':res.groupby('method')[['AUROC','AP','Recall@5','Recall@10','Recall@20']].mean().to_dict(orient='index')
}
(OUT/'anm_summary.json').write_text(json.dumps(summary,indent=2))
print(res.to_string(index=False))
print(json.dumps(summary,indent=2))
