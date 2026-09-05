from pathlib import Path
import pandas as pd, json
R=Path('validation/packages/all15/results_curated_validation_v2')
O=R/'decision_summary'; O.mkdir(exist_ok=True)
m=pd.read_csv(R/'all_metrics.csv')
pri=m[m.role=='PRIMARY'].copy()
methods=['Core','Core_GNM','Core_APBS','Core_GNM_APBS','IS_candidate','IS_RIN_degree','IS_RIN_closeness','IS_RIN_betweenness']
metrics=['AUROC','AP','Recall@5','Recall@10','Recall@20','Spatial@5A_top10','Spatial@8A_top10','Spatial@10A_top10','median_nearest_top10_A']
# Condition wide table
w=pri.pivot(index=['protein','condition_id'],columns='method',values=metrics)
w.columns=[f'{x}__{y}' for x,y in w.columns]; w.reset_index().to_csv(O/'per_condition_wide.csv',index=False)
# Protein macro table (median across conditions in same protein)
pm=pri.groupby(['protein','method'])[metrics].median().reset_index()
pm.to_csv(O/'protein_macro_by_method.csv',index=False)
# Deltas vs core by protein and method
base=pm[pm.method=='Core'].set_index('protein')
rows=[]
for method in methods[1:]:
 cur=pm[pm.method==method].set_index('protein')
 for p in sorted(base.index):
  row={'protein':p,'method':method}
  for x in metrics: row[f'delta_{x}']=cur.loc[p,x]-base.loc[p,x]
  rows.append(row)
pd.DataFrame(rows).to_csv(O/'protein_deltas_vs_core.csv',index=False)
# Stress and negative control
m[m.role!='PRIMARY'].to_csv(O/'stress_and_negative_metrics.csv',index=False)
# Compact decision matrix
agg=pd.read_csv(R/'primary_aggregate.csv').set_index('method')
d=pd.read_csv(R/'component_deltas_by_protein.csv')
dec=[]
for method in methods:
 a=agg.loc[method]
 r={'method':method,'protein_macro_AUROC':a['protein_macro_median_AUROC'],'protein_macro_AP':a['protein_macro_median_AP'],'protein_macro_R5':a['protein_macro_median_Recall@5'],'protein_macro_R10':a['protein_macro_median_Recall@10'],'protein_macro_R20':a['protein_macro_median_Recall@20'],'protein_macro_Spatial10':a['protein_macro_median_Spatial@10A_top10'],'protein_macro_nearestA':a['protein_macro_median_median_nearest_top10_A']}
 if method!='Core':
  for x in ['AUROC','AP','Recall@10','Recall@20','Spatial@10A_top10']:
   q=d[(d.method==method)&(d.metric==x)].iloc[0];r[f'{x}_proteins_improved']=int(q.proteins_improved);r[f'{x}_proteins_worsened']=int(q.proteins_worsened);r[f'{x}_median_delta']=q.observed_median_delta_vs_Core;r[f'{x}_delta_CI']=f"[{q['bootstrap_ci2.5']:.4f},{q['bootstrap_ci97.5']:.4f}]"
 dec.append(r)
pd.DataFrame(dec).to_csv(O/'model_decision_matrix.csv',index=False)
print(pd.DataFrame(dec).to_string(index=False))
