from pathlib import Path
import tempfile, json
import numpy as np
import pandas as pd

# Reuse the source-corrected benchmark definitions, structure filters, metrics, GNM/RIN implementation.
ns={'__file__':'validation/run_curated_validation_v2.py','__name__':'curated_v2_definitions'}
src=Path('validation/run_curated_validation_v2.py').read_text()
exec(src.split('# Structure QC before metrics.')[0],ns)

ROOT=ns['ROOT']; OUT=ROOT/'results_curated_validation_v3_fast'; OUT.mkdir(parents=True,exist_ok=True)
PDBDIR=OUT/'pdb_cache'; ns['PDBDIR']=PDBDIR; PDBDIR.mkdir(exist_ok=True)
ism=ns['ism']; pct=ns['pct']; PRIMARY=ns['PRIMARY']; STRESS=ns['STRESS']; NEG=ns['NEG']; structure=ns['structure']; getpdb=ns['getpdb']
auc=ns['auc']; ap=ns['ap']; recall=ns['recall']; spatial=ns['spatial']; nearest=ns['nearest']; SC={}; BASIC={}; CHARGED={}
SASA_POINTS=200

def basic_env(pid,ph):
 key=(pid,ph)
 if key in BASIC:return BASIC[key]
 old=ism.SASA_POINTS;ism.SASA_POINTS=SASA_POINTS
 try: _,res,atoms,_=ism.build_surface_residues(getpdb(pid),ph)
 finally: ism.SASA_POINTS=old
 surf=[r for r in res if r['surface_exposed']]; BASIC[key]=(surf,ism.build_distances(surf),atoms);return BASIC[key]

def charged_env(pid,ph,ionic,temp):
 key=(pid,ph,ionic,temp)
 if key in CHARGED:return CHARGED[key]
 old=ism.SASA_POINTS;ism.SASA_POINTS=SASA_POINTS
 try:
  _,res,atoms,_=ism.build_surface_residues(getpdb(pid),ph);e=ism.EnvParams(pH=ph,ionic=ionic,temp=temp)
  with tempfile.TemporaryDirectory(prefix='is_v3_') as td:status=ism.attach_apbs_auxiliary(getpdb(pid),res,atoms,e,Path(td))
 finally: ism.SASA_POINTS=old
 surf=[r for r in res if r['surface_exposed']];CHARGED[key]=(surf,ism.build_distances(surf),status);return CHARGED[key]

def score(case):
 pid=case['pdb'];ph=case['pH'];chem=case['chem'];ionic=case['ionic'];temp=case['temp']
 if pid not in SC:SC[pid]=structure(pid)[0]
 expected=ism.CHEMISTRIES[chem].get('expected_phi_sign')
 if expected in ('positive','negative'):
  surf,D,status=charged_env(pid,ph,ionic,temp);apbs_used=True
 else:
  surf,D,_=basic_env(pid,ph);status='not_run_non_electrostatic_channel';apbs_used=False
 mp=ism.chemistry_map(surf,D,chem,ph);core={x['center_key']:x['multiscale_persistence']/100 for x in mp['patch_centers']};rows=[]
 for r in surf:
  phi=r.get('phi');comp=None
  if expected=='positive' and phi is not None:comp=max(float(phi),0.)
  elif expected=='negative' and phi is not None:comp=max(-float(phi),0.)
  z={'key':r['key'],'res_name':r['res_name'],'res_seq':int(r['res_seq']),'chain':r['chain'],'scrsa':float(r['scrsa']),'core_raw':float(core.get(r['key'],0.)),'phi':phi,'apbs_compat':comp,'x':r['x'],'y':r['y'],'z':r['z']};z.update(SC[pid].get(r['key'],{}));rows.append(z)
 df=pd.DataFrame(rows);cr=pct(df.core_raw);gr=pct(df.gnm);rr=pct(df.radial)
 ar=pct(df.apbs_compat.fillna(0)) if apbs_used and df.apbs_compat.notna().any() else np.ones(len(df))
 dr=1-pct(df.rin_degree);clr=1-pct(df.rin_closeness);br=1-pct(df.rin_betweenness)
 vals={'Core':cr,'Core_GNM':cr*gr,'Core_APBS':cr*ar,'Core_GNM_APBS':cr*gr*ar,'IS_candidate':cr*gr*ar*rr,'IS_RIN_degree':cr*gr*ar*rr*dr,'IS_RIN_closeness':cr*gr*ar*rr*clr,'IS_RIN_betweenness':cr*gr*ar*rr*br}
 for k,v in vals.items():df[k]=100*v/v.max() if len(v) and np.max(v)>0 else np.zeros(len(v))
 return df,status,apbs_used

def evaluate(case,role):
 df,status,apbs_used=score(case);truth=set(case['truth']);truthkeys=df[df.res_seq.isin(truth)].key.tolist();y=df.key.isin(truthkeys).astype(int).to_numpy();coords={r.key:np.array([r.x,r.y,r.z],float) for _,r in df.iterrows()};rows=[]
 for method in ['Core','Core_GNM','Core_APBS','Core_GNM_APBS','IS_candidate','IS_RIN_degree','IS_RIN_closeness','IS_RIN_betweenness']:
  s=df[method].to_numpy(float);o=np.argsort(-s,kind='stable');top=df.key.iloc[o[:10]].tolist();rows.append({'condition_id':case['id'],'protein':case['protein'],'pdb':case['pdb'],'role':role,'method':method,'ionic_mM':case['ionic'],'pH':case['pH'],'n_surface':len(df),'n_truth_requested':len(truth),'n_truth_mapped':int(y.sum()),'AUROC':auc(y,s),'AP':ap(y,s),'Recall@5':recall(y,s,5),'Recall@10':recall(y,s,10),'Recall@20':recall(y,s,20),'Spatial@5A_top10':spatial(truthkeys,top,coords,5),'Spatial@8A_top10':spatial(truthkeys,top,coords,8),'Spatial@10A_top10':spatial(truthkeys,top,coords,10),'median_nearest_top10_A':nearest(truthkeys,top,coords),'apbs_status':status,'apbs_used':apbs_used,'source':case['source'],'note':case.get('condition_note',case.get('note',''))})
 df['condition_id']=case['id'];df['protein']=case['protein'];df['role']=role;df['truth']=y;return rows,df

qc=[]
for pid in ns['CHAIN_RULES']:
 st,n=structure(pid);p=getpdb(pid);s=ns['PDBParser'](QUIET=True).get_structure(pid,str(p));m=next(s.get_models());chains=sorted({c.id for c in m});seqs=[int(r.id[1]) for c in m for r in c if ns['is_aa'](r,standard=True)]
 qc.append({'pdb':pid,'chains':';'.join(chains),'n_CA':n,'min_resseq':min(seqs),'max_resseq':max(seqs),'range_rule':str(ns['RANGE_RULES'].get(pid,'all'))})
pd.DataFrame(qc).to_csv(OUT/'structure_qc.csv',index=False)

allm=[];alld=[]
for c in PRIMARY:
 r,d=evaluate(c,'PRIMARY');allm+=r;alld.append(d)
for c in STRESS:
 r,d=evaluate(c,'STRESS');allm+=r;alld.append(d)
r,d=evaluate(NEG,'NEGATIVE_CONTROL');allm+=r;alld.append(d)
MET=pd.DataFrame(allm);MET.to_csv(OUT/'all_metrics.csv',index=False);pd.concat(alld,ignore_index=True).to_csv(OUT/'residue_scores.csv',index=False)
PRI=MET[MET.role=='PRIMARY'].copy();metrics=['AUROC','AP','Recall@5','Recall@10','Recall@20','Spatial@5A_top10','Spatial@8A_top10','Spatial@10A_top10','median_nearest_top10_A']
agg=[]
for method,g in PRI.groupby('method'):
 pm=g.groupby('protein')[metrics].median();row={'method':method,'n_conditions':g.condition_id.nunique(),'n_proteins':g.protein.nunique()}
 for x in metrics:
  row[f'condition_median_{x}']=g[x].median();row[f'condition_mean_{x}']=g[x].mean();row[f'protein_macro_median_{x}']=pm[x].median();row[f'protein_macro_mean_{x}']=pm[x].mean()
 agg.append(row)
AGG=pd.DataFrame(agg);AGG.to_csv(OUT/'primary_aggregate.csv',index=False)

rng=np.random.default_rng(20260905);proteins=sorted(PRI.protein.unique());B=10000;boots=[];deltas=[]
for method,g in PRI.groupby('method'):
 pm=g.groupby('protein')[metrics].median()
 for x in metrics:
  arr=np.empty(B)
  for b in range(B):
   samp=rng.choice(proteins,size=len(proteins),replace=True);arr[b]=np.median([pm.loc[p,x] for p in samp])
  boots.append({'method':method,'metric':x,'protein_macro_median':float(pm[x].median()),'ci2.5':float(np.quantile(arr,.025)),'ci97.5':float(np.quantile(arr,.975))})
pd.DataFrame(boots).to_csv(OUT/'protein_clustered_bootstrap.csv',index=False)
base=PRI[PRI.method=='Core'].groupby('protein')[metrics].median()
for method in [m for m in PRI.method.unique() if m!='Core']:
 cur=PRI[PRI.method==method].groupby('protein')[metrics].median()
 for x in metrics:
  d=(cur[x]-base[x]).dropna();arr=np.empty(B)
  for b in range(B):arr[b]=np.median(rng.choice(d.to_numpy(),size=len(d),replace=True))
  deltas.append({'method':method,'metric':x,'median_delta_vs_Core':float(np.median(d)),'proteins_improved':int((d>0).sum()),'proteins_worsened':int((d<0).sum()),'proteins_tied':int((d==0).sum()),'ci2.5':float(np.quantile(arr,.025)),'ci97.5':float(np.quantile(arr,.975))})
pd.DataFrame(deltas).to_csv(OUT/'component_deltas_by_protein.csv',index=False)
loo=[]
for left in proteins:
 for method,g in PRI.groupby('method'):
  hold=g[g.protein==left];rest=g[g.protein!=left]
  for x in metrics:loo.append({'left_out_protein':left,'method':method,'metric':x,'holdout_condition_median':hold[x].median(),'remaining_protein_macro_median':rest.groupby('protein')[x].median().median()})
pd.DataFrame(loo).to_csv(OUT/'leave_one_protein_out.csv',index=False)
meta={'version':'v3-fast-scientifically-equivalent-primary','primary_conditions':len(PRIMARY),'primary_proteins':len(proteins),'stress_conditions':len(STRESS),'negative_controls':1,'SASA_POINTS':200,'GNM_cutoff_A':10.0,'RIN_cutoff_A':8.0,'APBS_policy':'run only for chemistry channels with expected electrostatic sign; neutral/not run for hydrophobic, pi-carbon, oxide','FN':'1FNF A:1236-1509; pH7; COO 1M; NH3 0.8M; CH3 0.05M; 300K','WNV':'2HG0 A:8-109; frequent-contact GT 55-61,106-109','PLOS':'Cyt-c/RNase pH7.4 source-aligned','statistical_unit':'protein; repeated conditions within protein summarized before macro inference','status':'development benchmark, not prospective validation'}
(OUT/'METADATA.json').write_text(json.dumps(meta,indent=2))
print(json.dumps(meta,indent=2));print(AGG[['method','n_conditions','n_proteins','protein_macro_median_AUROC','protein_macro_median_AP','protein_macro_median_Recall@10','protein_macro_median_Recall@20']].to_string(index=False))
