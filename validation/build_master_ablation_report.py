from pathlib import Path
import pandas as pd

ROOT=Path('validation/packages/all15')
base=pd.read_csv(ROOT/'results_is_apbs_gnm_exploratory'/'is_apbs_gnm_method_comparison.csv')
geom=pd.read_csv(ROOT/'results_light_geometry_exploratory'/'light_geometry_method_comparison.csv')
rin=pd.read_csv(ROOT/'results_rin_light_exploratory'/'rin_method_comparison.csv')
reg=pd.read_csv(ROOT/'condition_registry.csv')
prot=pd.read_csv(ROOT/'protein_manifest.csv')

selected=[]

def take(df, method, label, stage):
    x=df[df['method']==method].copy()
    if x.empty: return
    x['model_label']=label; x['stage']=stage
    selected.append(x)

take(base,'IS_static','Core IS','core')
take(base,'IS_x_GNM','Core + GNM','gnm')
take(base,'IS_x_APBS','Core + conditional APBS','apbs')
take(base,'IS_x_GNM_x_APBS','Core + GNM + conditional APBS','gnm_apbs')
take(geom,'base_x_radial','Core + GNM + conditional APBS + radial prominence','radial')
take(rin,'IS_x_degree_peripheral','Current IS + RIN degree peripheralness','rin_degree')
take(rin,'IS_x_closeness_peripheral','Current IS + RIN closeness peripheralness','rin_closeness')
take(rin,'IS_x_betweenness_peripheral','Current IS + RIN betweenness peripheralness','rin_betweenness')

M=pd.concat(selected,ignore_index=True,sort=False)
keep=['condition_id','protein','pdb','model_label','stage','AUROC','AP','Recall@5','Recall@10','Recall@20']
M=M[[c for c in keep if c in M.columns]].copy()
core=M[M.stage=='core'][['condition_id','AUROC','AP','Recall@5','Recall@10','Recall@20']].rename(columns={c:f'core_{c}' for c in ['AUROC','AP','Recall@5','Recall@10','Recall@20']})
M=M.merge(core,on='condition_id',how='left')
for c in ['AUROC','AP','Recall@5','Recall@10','Recall@20']:
    M[f'delta_{c}_vs_core']=M[c]-M[f'core_{c}']

# add registry metadata
rcols=['condition_id','surface_or_condition','chemistry_channel','tier','ground_truth_type','ground_truth_status','analysis_role','package_status']
M=M.merge(reg[rcols],on='condition_id',how='left')
OUT=ROOT/'master_ablation'; OUT.mkdir(exist_ok=True)
M.to_csv(OUT/'master_ablation_long.csv',index=False)

# wide AUROC/AP delta table
wide=M.pivot_table(index=['condition_id','protein','pdb','surface_or_condition','chemistry_channel'],columns='stage',values=['AUROC','AP','Recall@10','Recall@20'],aggfunc='first')
wide.columns=[f'{a}__{b}' for a,b in wide.columns]
wide=wide.reset_index()
wide.to_csv(OUT/'master_ablation_wide.csv',index=False)

# package status of all 15 proteins + all conditions
reg.to_csv(OUT/'all_condition_status.csv',index=False)
prot.to_csv(OUT/'all15_protein_status.csv',index=False)

# consistency counts vs core
rows=[]
for stage,g in M[M.stage!='core'].groupby('stage'):
    rows.append({
        'stage':stage,
        'n_conditions':len(g),
        'AUROC_improved':int((g['delta_AUROC_vs_core']>0).sum()),
        'AUROC_worsened':int((g['delta_AUROC_vs_core']<0).sum()),
        'AP_improved':int((g['delta_AP_vs_core']>0).sum()),
        'AP_worsened':int((g['delta_AP_vs_core']<0).sum()),
        'R10_improved':int((g['delta_Recall@10_vs_core']>0).sum()),
        'R10_worsened':int((g['delta_Recall@10_vs_core']<0).sum()),
        'R20_improved':int((g['delta_Recall@20_vs_core']>0).sum()),
        'R20_worsened':int((g['delta_Recall@20_vs_core']<0).sum()),
        'median_AUROC':g['AUROC'].median(),'median_AP':g['AP'].median(),'median_R10':g['Recall@10'].median(),'median_R20':g['Recall@20'].median()
    })
cons=pd.DataFrame(rows); cons.to_csv(OUT/'component_consistency_summary.csv',index=False)

# markdown report
lines=['# InterfaceScout master ablation report','',
'## Interpretation rule','',
'All currently scored conditions are developmental/prefreeze evidence. Do not call them blind prospective validation. The 15-protein package contains Tier-1 residue-level, Tier-2 orientation/spatial, Tier-3 stress-test, and pending-audit conditions; only residue-level conditions with sufficiently explicit ground truth enter AUROC/AP tables.','',
'## Component consistency versus Core IS','',cons.to_markdown(index=False,floatfmt='.3f'),'','## Per-condition AUROC / AP','']
show=M[['condition_id','protein','model_label','AUROC','AP','Recall@10','Recall@20','delta_AUROC_vs_core','delta_AP_vs_core']]
lines.append(show.to_markdown(index=False,floatfmt='.3f'))
lines+=['','## All 15 proteins','',prot.to_markdown(index=False),'','## Condition registry','',reg[['condition_id','protein','surface_or_condition','tier','ground_truth_status','package_status']].to_markdown(index=False)]
(OUT/'MASTER_ABLATION_REPORT.md').write_text('\n'.join(lines))
print(cons.to_string(index=False))
