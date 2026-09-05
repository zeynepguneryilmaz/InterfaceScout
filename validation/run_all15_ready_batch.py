from __future__ import annotations
import json, math, urllib.request, sys
from pathlib import Path
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import backend.main as ism

OUT = Path('validation/packages/all15/results_ready_batch')
OUT.mkdir(parents=True, exist_ok=True)
PDBDIR = OUT / 'pdb_cache'
PDBDIR.mkdir(exist_ok=True)

# Pre-declared literature ground truth. No benchmark-specific retuning.
# Only conditions with explicit residue-level GT and already-supported chemistry channels are run here.
BENCH = [
    {'id':'FN_COO','protein':'Fibronectin III8-10','pdb':'1FNF','chains':['A'],'chem':'anionic','pH':7.0,
     'truth':[1469], 'source':'Materials 2018 doi:10.3390/ma11122570'},
    {'id':'FN_NH3','protein':'Fibronectin III8-10','pdb':'1FNF','chains':['A'],'chem':'cationic','pH':7.0,
     'truth':[1312,1509], 'source':'Materials 2018 doi:10.3390/ma11122570'},
    {'id':'FN_CH3_HEAD','protein':'Fibronectin III8-10','pdb':'1FNF','chains':['A'],'chem':'hydrophobic','pH':7.0,
     'truth':[1454,1455,1456,1457,1478,1479,1480,1481,1509], 'source':'Materials 2018 doi:10.3390/ma11122570'},
    {'id':'FN_CH3_SIDE','protein':'Fibronectin III8-10','pdb':'1FNF','chains':['A'],'chem':'hydrophobic','pH':7.0,
     'truth':[1275,1276,1355,1376,1378,1454,1455,1456,1457], 'source':'Materials 2018 doi:10.3390/ma11122570'},
    {'id':'FN_CH3_BETA','protein':'Fibronectin III8-10','pdb':'1FNF','chains':['A'],'chem':'hydrophobic','pH':7.0,
     'truth':[1431,1446,1457,1458,1459,1461,1463,1464,1466,1475,1476], 'source':'Materials 2018 doi:10.3390/ma11122570'},
    {'id':'CYTC_CH3','protein':'Cytochrome c','pdb':'3NWV','chains':['A'],'chem':'hydrophobic','pH':7.0,
     'truth':[1,2,3,4,96,99,100], 'source':'PLOS ONE 2014 doi:10.1371/journal.pone.0107696'},
    {'id':'CYTC_COOH','protein':'Cytochrome c','pdb':'3NWV','chains':['A'],'chem':'anionic','pH':7.0,
     'truth':[1,2,4,61,99,100,103], 'source':'PLOS ONE 2014 doi:10.1371/journal.pone.0107696'},
    {'id':'LYZ_SWCNT','protein':'Lysozyme','pdb':'1LYZ','chains':['A'],'chem':'pi_carbon','pH':7.4,
     'truth':[61,62,63,68,69,70,73,103,106,107,111,112,113,116], 'source':'Scientific Reports 2025 doi:10.1038/s41598-025-96435-3'},
    {'id':'CHT_CNT_HYDRO','protein':'Alpha-chymotrypsin','pdb':'4CHA','chains':['A','B','C'],'chem':'hydrophobic','pH':7.0,
     'truth':[3,5,6,7,8,10,77,114,115,116], 'source':'Scientific Reports 2015 doi:10.1038/srep09297'},
    {'id':'CHT_CNT_PI_NEG','protein':'Alpha-chymotrypsin','pdb':'4CHA','chains':['A','B','C'],'chem':'pi_carbon','pH':7.0,
     'truth':[3,5,6,7,8,10,77,114,115,116], 'source':'Scientific Reports 2015 doi:10.1038/srep09297'},
    {'id':'MB_CIT_AUNP','protein':'Myoglobin','pdb':'1MBN','chains':['A'],'chem':'anionic','pH':7.0,
     'truth':[43,45,96,97,98,99,146,147,148,149,150,151,152,153], 'source':'Int J Mol Sci 2019 doi:10.3390/ijms20143539'},
]

CHAINSETS={b['pdb']:set(b['chains']) for b in BENCH}

class ProteinFirstModelSelect(Select):
    def __init__(self,pid,first_model_id):
        self.pid=pid; self.first_model_id=first_model_id
    def accept_model(self,m): return 1 if m.id==self.first_model_id else 0
    def accept_chain(self,c): return 1 if str(c.id) in CHAINSETS[self.pid] else 0
    def accept_residue(self,r): return 1 if is_aa(r,standard=True) else 0

def getpdb(pid):
    raw=PDBDIR/f'{pid}_raw.pdb'; clean=PDBDIR/f'{pid}.pdb'
    if not raw.exists(): urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pid}.pdb', raw)
    s=PDBParser(QUIET=True).get_structure(pid,str(raw)); first_model=next(s.get_models())
    io=PDBIO(); io.set_structure(s); io.save(str(clean),ProteinFirstModelSelect(pid,first_model.id))
    return clean

def auc(y,s):
    p=s[y==1]; n=s[y==0]
    if not len(p) or not len(n): return np.nan
    return float(np.mean([(a>b)+0.5*(a==b) for a in p for b in n]))

def ap(y,s):
    if not y.sum(): return np.nan
    o=np.argsort(-s,kind='stable'); yy=y[o]
    prec=np.cumsum(yy)/(np.arange(len(yy))+1)
    return float((prec*yy).sum()/yy.sum())

def recall(y,s,k):
    if not y.sum(): return np.nan
    o=np.argsort(-s,kind='stable')[:min(k,len(s))]
    return float(y[o].sum()/y.sum())

def spatial(truth,top,coords,R):
    if not truth: return np.nan
    vals=[]
    for t in truth:
        if t not in coords: continue
        ds=[np.linalg.norm(coords[t]-coords[p]) for p in top if p in coords]
        vals.append(bool(ds and min(ds)<=R))
    return float(np.mean(vals)) if vals else np.nan

def nearest_patch_distance(truth,top,coords):
    ds=[]
    for t in truth:
        if t not in coords: continue
        q=[np.linalg.norm(coords[t]-coords[p]) for p in top if p in coords]
        if q: ds.append(min(q))
    return float(np.median(ds)) if ds else np.nan

def perm(y,s,exp,n=10000,seed=20260905):
    rng=np.random.default_rng(seed)
    obs=float(s[y==1].mean())
    q=pd.qcut(exp,4,labels=False,duplicates='drop')
    bins=np.unique(q)
    need={int(b):int(((q==b)&(y==1)).sum()) for b in bins}
    pools={b:np.where(q==b)[0] for b in need}
    null=np.empty(n)
    for z in range(n):
        idx=[]
        for b,c in need.items():
            if c: idx.extend(rng.choice(pools[b],c,replace=False).tolist())
        null[z]=s[np.asarray(idx,dtype=int)].mean() if idx else 0.0
    return float((1+(null>=obs).sum())/(n+1)),obs,float(null.mean())

rows=[]; detail=[]; cache={}
for b in BENCH:
    ck=(b['pdb'],tuple(b['chains']),b['pH'])
    if ck not in cache:
        old=ism.SASA_POINTS; ism.SASA_POINTS=200
        try:
            _,allres,_,_=ism.build_surface_residues(getpdb(b['pdb']),b['pH'])
        finally: ism.SASA_POINTS=old
        surf=[r for r in allres if r['surface_exposed']]
        cache[ck]=(surf,ism.build_distances(surf))
    surf,D=cache[ck]
    mp=ism.chemistry_map(surf,D,b['chem'],b['pH'])
    sm={x['center_key']:x['multiscale_persistence']/100.0 for x in mp['patch_centers']}
    keys=[r['key'] for r in surf]
    scores=np.asarray([sm.get(k,0.0) for k in keys],float)
    exp=np.asarray([r['scrsa'] for r in surf],float)
    coords={r['key']:np.asarray([r['x'],r['y'],r['z']],float) for r in surf}
    truthkeys=[r['key'] for r in surf if int(r['res_seq']) in b['truth']]
    y=np.asarray([int(k in truthkeys) for k in keys],int)
    order=np.argsort(-scores,kind='stable'); top=[keys[i] for i in order]
    ranks={k:top.index(k)+1 for k in truthkeys}
    pp,obs,null=perm(y,scores,exp) if y.sum() else (np.nan,np.nan,np.nan)
    row={**b,'n_surface':len(keys),'n_truth_requested':len(b['truth']),'n_truth_mapped':int(y.sum()),
         'AUROC':auc(y,scores),'AP':ap(y,scores),'Recall@5':recall(y,scores,5),'Recall@10':recall(y,scores,10),'Recall@20':recall(y,scores,20),
         'Spatial@5A_top10':spatial(truthkeys,top[:10],coords,5.0),'Spatial@8A_top10':spatial(truthkeys,top[:10],coords,8.0),'Spatial@10A_top10':spatial(truthkeys,top[:10],coords,10.0),
         'median_nearest_top10_patch_A':nearest_patch_distance(truthkeys,top[:10],coords),
         'permutation_p':pp,'truth_mean_M':obs,'null_mean_M':null,
         'best_truth_rank':min(ranks.values()) if ranks else None,'median_truth_rank':float(np.median(list(ranks.values()))) if ranks else None,
         'truth_ranks_json':json.dumps(ranks,sort_keys=True)}
    rows.append(row)
    for r,k,yy,ss in zip(surf,keys,y,scores):
        detail.append([b['id'],b['protein'],k,r['res_name'],r['res_seq'],r['scrsa'],yy,ss])

df=pd.DataFrame(rows)
df.to_csv(OUT/'ready_tier1_results.csv',index=False)
pd.DataFrame(detail,columns=['condition_id','protein','key','res_name','res_seq','scRSA','truth','M']).to_csv(OUT/'ready_tier1_residue_scores.csv',index=False)
(OUT/'ready_tier1_results.json').write_text(df.to_json(orient='records',indent=2))

# Macro summaries by condition and by unique protein (negative control excluded from primary aggregate)
primary=df[df['id']!='CHT_CNT_PI_NEG'].copy()
metric_cols=['AUROC','AP','Recall@5','Recall@10','Recall@20','Spatial@5A_top10','Spatial@8A_top10','Spatial@10A_top10','median_nearest_top10_patch_A','permutation_p']
condition_macro={m:float(primary[m].median()) for m in metric_cols}
protein_macro_df=primary.groupby('protein')[metric_cols].median(numeric_only=True)
protein_macro={m:float(protein_macro_df[m].median()) for m in metric_cols}
summary={'n_conditions_run':int(len(df)),'n_primary_conditions':int(len(primary)),'n_unique_proteins':int(df['protein'].nunique()),'condition_macro_medians':condition_macro,'protein_macro_medians':protein_macro}
(OUT/'ready_tier1_summary.json').write_text(json.dumps(summary,indent=2))
print(df[['id','protein','n_truth_mapped','AUROC','AP','Recall@10','Spatial@8A_top10','Spatial@10A_top10','permutation_p','best_truth_rank','median_truth_rank']].to_string(index=False))
print(json.dumps(summary,indent=2))
