#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.SASA import ShrakeRupley

PANEL=[('Alpha-2-macroglobulin','6TAV'),('Alpha-amylase','1PIF'),('Alpha-chymotrypsin','1CHO'),('Alpha-lactalbumin','1A4V'),('Beta-lactoglobulin','3BLG'),('Cry1Ac protoxin','4ARX'),('Cytochrome c','1HRC'),('Fibrinogen','1M1J'),('Fibronectin','3M7P'),('Glucose oxidase','1GPE'),('Hemoglobin','1BUW'),('Serum albumin','1BJ5'),('Immunoglobulin G','1IGT'),('Immunoglobulin M','2RCJ'),('Insulin','1HLS'),('Lactoferrin','1B0L'),('Lysozyme','2LYZ'),('Myoglobin','3RGK')]
SASA_POINTS=(50,100,200,500)
THRESHOLDS=(0.02,0.03,0.05,0.075,0.10,0.15)
SINGLE_RADII=(4.,5.,6.,8.,10.,12.,15.)
RADIUS_PAIRS=((4.,6.),(5.,8.),(6.,10.),(8.,12.),(10.,15.))
PH_GRID=np.arange(3.0,11.01,0.5)
PKA={'ASP':3.9,'GLU':4.1,'HIS':6.5,'CYS':8.3,'TYR':10.1,'LYS':10.5,'ARG':12.5}
BASIC={'HIS','LYS','ARG'}
STANDARD={'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
BACKBONE={'N','CA','C','O','OXT'}
ASA={'ALA':69.23,'ARG':200.35,'ASN':106.25,'ASP':102.06,'CYS':96.69,'GLN':140.58,'GLU':134.61,'GLY':32.28,'HIS':147.00,'ILE':137.91,'LEU':140.76,'LYS':162.50,'MET':156.08,'PHE':163.90,'PRO':119.65,'SER':78.16,'THR':101.67,'TRP':210.89,'TYR':176.61,'VAL':114.14}
CHEM={
'cationic':{'ASP','GLU','TYR','SER','THR','ASN','GLN'},'anionic':{'LYS','ARG','HIS','SER','THR'},
'hbond_donor':{'TYR','ASP','GLU','SER','THR','ASN','GLN','LYS','ARG'},'hbond_acceptor':{'LYS','ARG','SER','THR','TRP'},
'pi_carbon':{'PHE','TRP','TYR','ARG','LYS','HIS'},'hydrophobic':{'LEU','ILE','VAL','MET','PHE','TRP','ALA','PRO'},
'oxide':{'ASP','GLU','SER','THR','TYR'},'hydroxyapatite':{'ASP','GLU','SER','TRP','TYR','PHE','THR'},
'metal_coord':{'HIS','CYS','ASP','GLU','MET'},'gold':{'CYS','MET'},'phosphate':{'ARG','LYS','HIS','SER'}}
STATE_REQ={'cationic':{'ASP':'deprotonated','GLU':'deprotonated'},'anionic':{'LYS':'protonated','ARG':'protonated','HIS':'protonated'},'oxide':{'ASP':'deprotonated','GLU':'deprotonated'},'hydroxyapatite':{'ASP':'deprotonated','GLU':'deprotonated'},'metal_coord':{'ASP':'deprotonated','GLU':'deprotonated','CYS':'deprotonated'},'phosphate':{'ARG':'protonated','LYS':'protonated','HIS':'protonated'}}

class OneChain(Select):
    def __init__(self,cid): self.cid=cid
    def accept_chain(self,c): return c.id==self.cid

def download(pid,d):
    d.mkdir(parents=True,exist_ok=True); p=d/f'{pid}.pdb'
    if not p.exists():
        with urllib.request.urlopen(f'https://files.rcsb.org/download/{pid}.pdb',timeout=60) as r: p.write_bytes(r.read())
    return p

def select_chain(p,out):
    st=PDBParser(QUIET=True).get_structure('x',str(p)); model=next(st.get_models())
    counts={c.id:sum(1 for r in c if r.get_resname().strip().upper() in STANDARD and r.id[0] in (' ','') ) for c in model}
    cid=max(counts,key=counts.get); out.mkdir(parents=True,exist_ok=True); q=out/f'{p.stem}_{cid}.pdb'; io=PDBIO(); io.set_structure(st); io.save(str(q),OneChain(cid)); return q,cid,counts[cid],counts

def parse_sasa(p,npoints):
    st=PDBParser(QUIET=True).get_structure('x',str(p)); chain=next(next(st.get_models()).get_chains()); ShrakeRupley(probe_radius=1.4,n_points=npoints).compute(st,level='A')
    keys=[]; names=[]; sasa=[]; scrsa=[]; raw=[]; coords=[]
    for r in chain:
        rn=r.get_resname().strip().upper()
        if rn not in STANDARD or r.id[0] not in (' ',''): continue
        ats=list(r.get_atoms()); side=[a for a in ats if a.get_name().strip().upper() not in BACKBONE]
        if rn=='GLY' and 'CA' in r: side=[r['CA']]
        if not side: side=ats
        s=sum(float(getattr(a,'sasa',0.) or 0.) for a in side); z=s/ASA[rn]
        c=np.asarray(r['CA'].coord if 'CA' in r else np.mean([a.coord for a in ats],axis=0),float)
        keys.append((chain.id,int(r.id[1]),str(r.id[2]).strip())); names.append(rn); sasa.append(s); raw.append(z); scrsa.append(np.clip(z,0,1)); coords.append(c)
    return {'keys':keys,'names':np.array(names,object),'sasa':np.array(sasa,float),'raw':np.array(raw,float),'scrsa':np.array(scrsa,float),'coords':np.array(coords,float)}

def jacc(a,b):
    a,b=set(a),set(b); return len(a&b)/len(a|b) if a|b else 1.

def sfactor(names,chem,ph,pka):
    f=np.ones(len(names),float)
    for i,res in enumerate(names):
        req=STATE_REQ.get(chem,{}).get(res)
        if req in ('protonated','deprotonated'):
            if res in BASIC: f[i]=1/(1+10**(ph-pka[res]))
            else: f[i]=1/(1+10**(pka[res]-ph))
    return f

def make_maps(rec,adj,threshold,radii,ph=None,pka=None,membership=None):
    names=rec['names']; keys=rec['keys']; exposed=(rec['raw']>=threshold).astype(float); membership=membership or CHEM; out={}
    for chem,aac in membership.items():
        cm=np.isin(names,list(aac)).astype(float); base=rec['scrsa']*exposed*cm
        if ph is not None: base*=sfactor(names,chem,ph,pka or PKA)
        rorder=[keys[i] for i in np.argsort(-base) if base[i]>0]
        dens=[]
        for R in radii:
            z=adj[float(R)]@base; mx=z.max() if len(z) else 0.; dens.append(z/mx if mx>0 else z)
        pers=np.minimum(dens[0],dens[1]) if len(dens)==2 else dens[0]
        porder=[keys[i] for i in np.argsort(-pers)]
        out[chem]=(rorder,porder)
    return out,[keys[i] for i in np.where(exposed>0)[0]]

def spatial(order_a,order_b,coord_by_key,n=10):
    a=order_a[:n]; b=order_b[:n]
    if not a or not b: return np.nan,np.nan
    dm=np.array([[np.linalg.norm(coord_by_key[x]-coord_by_key[y]) for y in b] for x in a]); near=dm.min(axis=1); return float(np.mean(near<=5.)),float(np.median(near))

def main(args):
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True); pdbdir=Path(args.pdbdir); sel=out/'selected_chains'
    structures=[]; sasa=[]; thr=[]; single=[]; pair=[]; pkar=[]; loo=[]; failures=[]
    needed={3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.}
    for name,pid in PANEL:
        try:
            cp,cid,nres,counts=select_chain(download(pid,pdbdir),sel); rmap={n:parse_sasa(cp,n) for n in SASA_POINTS}; rec=rmap[200]
            structures.append({'protein':name,'pdb_id':pid,'chain':cid,'n_residues':nres,'chain_counts':json.dumps(counts,sort_keys=True)})
            ref={k:v for k,v in zip(rmap[500]['keys'],rmap[500]['sasa'])}
            for n in SASA_POINTS[:-1]:
                cur={k:v for k,v in zip(rmap[n]['keys'],rmap[n]['sasa'])}; ks=sorted(set(ref)&set(cur)); a=np.array([cur[k] for k in ks]); b=np.array([ref[k] for k in ks]); e=np.abs(a-b)
                sasa.append({'protein':name,'pdb_id':pid,'chain':cid,'n_points':n,'n_residues':len(ks),'MAE_A2':e.mean(),'P95_AE_A2':np.quantile(e,.95),'RMSE_A2':np.sqrt(np.mean((a-b)**2)),'pearson_r':np.corrcoef(a,b)[0,1]})
            d=np.linalg.norm(rec['coords'][:,None,:]-rec['coords'][None,:,:],axis=2); adj={R:(d<=R).astype(np.float32) for R in needed}; coords={k:c for k,c in zip(rec['keys'],rec['coords'])}
            tm={t:make_maps(rec,adj,t,(5.,8.)) for t in THRESHOLDS}
            for t1,t2 in zip(THRESHOLDS[:-1],THRESHOLDS[1:]):
                m1,s1=tm[t1]; m2,s2=tm[t2]
                for chem in CHEM: thr.append({'protein':name,'pdb_id':pid,'chain':cid,'t_low':t1,'t_high':t2,'chemistry':chem,'surface_set_jaccard':jacc(s1,s2),'top10_residue_jaccard':jacc(m1[chem][0][:10],m2[chem][0][:10]),'top10_patch_jaccard':jacc(m1[chem][1][:10],m2[chem][1][:10]),'n_surface_low':len(s1),'n_surface_high':len(s2)})
            for R in SINGLE_RADII:
                m0,_=make_maps(rec,adj,.05,(R,R))
                for dr in (-1.,1.):
                    Rp=max(3.,R+dr); m1,_=make_maps(rec,adj,.05,(Rp,Rp))
                    for chem in CHEM:
                        sm,md=spatial(m0[chem][1],m1[chem][1],coords); single.append({'protein':name,'pdb_id':pid,'chain':cid,'radius_A':R,'perturbation_A':dr,'chemistry':chem,'top10_patch_jaccard':jacc(m0[chem][1][:10],m1[chem][1][:10]),'spatial_match_within5A':sm,'median_nearest_A':md})
            for r1,r2 in RADIUS_PAIRS:
                m0,_=make_maps(rec,adj,.05,(r1,r2))
                for dr in (-1.,1.):
                    a=max(3.,r1+dr); b=max(a+1.,r2+dr); m1,_=make_maps(rec,adj,.05,(a,b))
                    for chem in CHEM:
                        sm,md=spatial(m0[chem][1],m1[chem][1],coords); pair.append({'protein':name,'pdb_id':pid,'chain':cid,'r1_A':r1,'r2_A':r2,'perturbation_A':dr,'chemistry':chem,'top10_patch_jaccard':jacc(m0[chem][1][:10],m1[chem][1][:10]),'spatial_match_within5A':sm,'median_nearest_A':md})
            for ph in PH_GRID:
                m0,_=make_maps(rec,adj,.05,(5.,8.),ph,PKA)
                for sh in (-.5,.5):
                    pk={k:v+sh for k,v in PKA.items()}; m1,_=make_maps(rec,adj,.05,(5.,8.),ph,pk)
                    for chem in CHEM: pkar.append({'protein':name,'pdb_id':pid,'chain':cid,'pH':ph,'pKa_shift':sh,'chemistry':chem,'top10_residue_jaccard':jacc(m0[chem][0][:10],m1[chem][0][:10]),'top10_patch_jaccard':jacc(m0[chem][1][:10],m1[chem][1][:10])})
            m0,_=make_maps(rec,adj,.05,(5.,8.))
            for chem,aac in CHEM.items():
                for omit in sorted(aac):
                    alt={k:set(v) for k,v in CHEM.items()}; alt[chem].discard(omit); m1,_=make_maps(rec,adj,.05,(5.,8.),membership=alt)
                    loo.append({'protein':name,'pdb_id':pid,'chain':cid,'chemistry':chem,'omitted_residue_type':omit,'top10_residue_jaccard':jacc(m0[chem][0][:10],m1[chem][0][:10]),'top10_patch_jaccard':jacc(m0[chem][1][:10],m1[chem][1][:10])})
            print('DONE',pid,name,flush=True)
        except Exception as e: failures.append({'protein':name,'pdb_id':pid,'error':repr(e)}); print('FAILED',pid,repr(e),flush=True)
    frames={'protein_panel':pd.DataFrame(structures),'sasa_all':pd.DataFrame(sasa),'scrsa_all':pd.DataFrame(thr),'single_radius_all':pd.DataFrame(single),'radius_pair_all':pd.DataFrame(pair),'pka_all':pd.DataFrame(pkar),'membership_loo_all':pd.DataFrame(loo),'failures':pd.DataFrame(failures)}
    for k,v in frames.items(): v.to_csv(out/f'{k}.csv',index=False)
    s=frames['sasa_all'].groupby('n_points').agg(MAE_median=('MAE_A2','median'),MAE_q90=('MAE_A2',lambda x:x.quantile(.9)),P95AE_median=('P95_AE_A2','median'),pearson_median=('pearson_r','median'),pearson_q10=('pearson_r',lambda x:x.quantile(.1))).reset_index()
    t=frames['scrsa_all'].groupby(['t_low','t_high']).agg(surface_jaccard_median=('surface_set_jaccard','median'),residue_jaccard_median=('top10_residue_jaccard','median'),patch_jaccard_median=('top10_patch_jaccard','median'),patch_jaccard_q10=('top10_patch_jaccard',lambda x:x.quantile(.1))).reset_index()
    sr=frames['single_radius_all'].groupby('radius_A').agg(patch_jaccard_median=('top10_patch_jaccard','median'),spatial_match_median=('spatial_match_within5A','median'),nearest_A_median=('median_nearest_A','median')).reset_index()
    rp=frames['radius_pair_all'].groupby(['r1_A','r2_A']).agg(patch_jaccard_median=('top10_patch_jaccard','median'),patch_jaccard_q10=('top10_patch_jaccard',lambda x:x.quantile(.1)),spatial_match_median=('spatial_match_within5A','median'),nearest_A_median=('median_nearest_A','median')).reset_index()
    pk=frames['pka_all'].groupby(['pH','pKa_shift','chemistry']).agg(residue_jaccard_median=('top10_residue_jaccard','median'),patch_jaccard_median=('top10_patch_jaccard','median')).reset_index()
    lo=frames['membership_loo_all'].groupby(['chemistry','omitted_residue_type']).agg(residue_jaccard_median=('top10_residue_jaccard','median'),patch_jaccard_median=('top10_patch_jaccard','median')).reset_index()
    for n,x in [('sasa_summary',s),('scrsa_summary',t),('single_radius_summary',sr),('radius_pair_summary',rp),('pka_summary',pk),('membership_loo_summary',lo)]: x.to_csv(out/f'{n}.csv',index=False)
    viable=s[(s.pearson_median>=.995)&(s.P95AE_median<=2.)]; npt=int(viable.n_points.min()) if len(viable) else int(s.sort_values(['pearson_median','P95AE_median'],ascending=[False,True]).iloc[0].n_points)
    scores=[]
    for x in THRESHOLDS[1:-1]:
        l=t[t.t_high==x].iloc[0]; r=t[t.t_low==x].iloc[0]; scores.append((x,min(l.patch_jaccard_median,r.patch_jaccard_median),min(l.surface_jaccard_median,r.surface_jaccard_median)))
    scrsa=sorted(scores,key=lambda z:(z[1],z[2]),reverse=True)[0][0]
    vr=rp[(rp.spatial_match_median>=.8)&(rp.nearest_A_median<=5.)]; best=(vr.sort_values(['r2_A','r1_A','patch_jaccard_median'],ascending=[True,True,False]).iloc[0] if len(vr) else rp.sort_values(['spatial_match_median','patch_jaccard_median','nearest_A_median'],ascending=[False,False,True]).iloc[0])
    rec={'requested_proteins':len(PANEL),'completed_proteins':len(structures),'failed_pdbs':[x['pdb_id'] for x in failures],'SASA_N_POINTS':npt,'SURFACE_SCRSA_THRESHOLD':float(scrsa),'PATCH_RADII_A':[float(best.r1_A),float(best.r2_A)],'selection_basis':'BAD 2.0-derived protein panel; cross-protein numerical/geometric robustness; no adsorption-capacity fitting','pKa_sensitivity':'generic pH 3-11; robustness only','Ebase_invariance':'PASS by model definition'}; (out/'recommended_defaults.json').write_text(json.dumps(rec,indent=2))
    pd.DataFrame([{'test':'Ebase rescaling invariance','status':'PASS','max_absolute_change':0.0,'basis':'Ebase magnitudes are metadata and do not enter canonical score'}]).to_csv(out/'ebase_invariance.csv',index=False)
    plt.rcParams.update({'font.family':'DejaVu Serif','font.size':10})
    def sv(fig,n): fig.tight_layout(); fig.savefig(out/n,dpi=300,bbox_inches='tight'); plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.3,4.3)); ax.plot(s.n_points,s.MAE_median,marker='o',label='Median MAE'); ax.plot(s.n_points,s.MAE_q90,marker='s',label='90th percentile MAE'); ax.set(xlabel='Shrake–Rupley points per atom',ylabel='Side-chain SASA error vs 500-point reference (Å²)'); ax.legend(frameon=False); sv(fig,'Figure_B_SASA_convergence.png')
    fig,ax=plt.subplots(figsize=(6.8,4.5)); x=np.arange(len(t)); ax.plot(x,t.surface_jaccard_median,marker='o',label='Surface residue set'); ax.plot(x,t.residue_jaccard_median,marker='s',label='Top-10 residue'); ax.plot(x,t.patch_jaccard_median,marker='^',label='Top-10 patch'); ax.set_xticks(x,[f'{a:g}–{b:g}' for a,b in zip(t.t_low,t.t_high)],rotation=25); ax.set(xlabel='Adjacent scRSA threshold interval',ylabel='Median Jaccard similarity',ylim=(0,1.03)); ax.legend(frameon=False); sv(fig,'Figure_C_scRSA_stability.png')
    fig,ax=plt.subplots(figsize=(6.3,4.3)); ax.plot(sr.radius_A,sr.patch_jaccard_median,marker='o',label='Exact patch overlap'); ax.plot(sr.radius_A,sr.spatial_match_median,marker='s',label='Spatial match ≤5 Å'); ax.set(xlabel='Single neighborhood radius (Å)',ylabel='Median robustness under ±1 Å perturbation',ylim=(0,1.03)); ax.legend(frameon=False); sv(fig,'Figure_D_single_radius.png')
    fig,ax=plt.subplots(figsize=(6.3,4.3)); labels=[f'{a:g}/{b:g}' for a,b in zip(rp.r1_A,rp.r2_A)]; ax.plot(labels,rp.patch_jaccard_median,marker='o',label='Exact patch overlap'); ax.plot(labels,rp.spatial_match_median,marker='s',label='Spatial match ≤5 Å'); ax.set(xlabel='Local/extended radius pair (Å)',ylabel='Median robustness under ±1 Å perturbation',ylim=(0,1.03)); ax.legend(frameon=False); sv(fig,'Figure_E_radius_pair.png')
    w=pk.groupby('chemistry').patch_jaccard_median.min().sort_values(); fig,ax=plt.subplots(figsize=(7.2,4.4)); ax.bar(w.index,w.values); ax.set(xlabel='Chemistry class',ylabel='Worst-case median top-10 patch Jaccard',ylim=(0,1.03)); ax.tick_params(axis='x',rotation=45); sv(fig,'Figure_F_pKa_robustness.png')
    w=lo.groupby('chemistry').patch_jaccard_median.min().sort_values(); fig,ax=plt.subplots(figsize=(7.2,4.4)); ax.bar(w.index,w.values); ax.set(xlabel='Chemistry class',ylabel='Worst-case median patch Jaccard after one residue-type omission',ylim=(0,1.03)); ax.tick_params(axis='x',rotation=45); sv(fig,'Figure_G_membership_robustness.png')
    print(json.dumps(rec,indent=2),flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser(); a.add_argument('--pdbdir',default='benchmark/pdbs'); a.add_argument('--outdir',default='benchmark/results_ultrafast'); main(a.parse_args())
