"""Independent production-model sensitivity and paired PB validation."""
import os,sys,csv,json,hashlib,traceback,shutil,math
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'backend'),str(ROOT/'validation')]
import core,coarse_patch,interface_engine as ie,model_settings as ms
assert ms.MODEL_VERSION=='2.5.0-dual-accessibility'
mode,shard,shards,reference,out=sys.argv[1:];shard,shards=int(shard),int(shards)
ref=Path(reference).resolve();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
def unique_file(name):
    p=list(ref.rglob(name));assert len(p)==1,(name,p);return p[0]
def sha(b):return hashlib.sha256(b).hexdigest()
def csvread(path):
    with path.open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))
def savecsv(path,rows):
    if rows:
        with path.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def dump(path,obj):
    path.write_text(json.dumps(obj,indent=2,allow_nan=False))
def jac(x,y):return len(x&y)/len(x|y) if x|y else 1.
def members(m,k):return {x for p in m['patches'][:k] for x in p['members']}
def centers(m,k):return {p['center_key'] for p in m['patches'][:k]}
BASE=(200,(6.,9.),8.,.05,7.4)
CONFIG=[('baseline',BASE)]
CONFIG += [(f'SASA_{n}',(n,*BASE[1:])) for n in (100,500,1000)]
CONFIG += [(f'radii_{r[0]}_{r[1]}',(200,r,*BASE[2:])) for r in ms.RADIUS_SELECTION_CANDIDATES_A if r!=BASE[1]]
CONFIG += [(f'patch_{r}',(*BASE[:2],r,*BASE[3:])) for r in (6.,10.,12.)]
CONFIG += [(f'tRSA_{v}',(*BASE[:3],v,7.4)) for v in (.01,.03,.075,.10)]
CONFIG += [(f'pH_{v}',(*BASE[:4],v)) for v in (5.5,9.)]
assert len(CONFIG)==22
def configure(c):
    n,r,p,t,ph=c
    core.SASA_POINTS=ie.SASA_POINTS=n
    core.PATCH_RADII_A=ie.MULTISCALE_RADII_A=r
    core.COARSE_PATCH_RADIUS_A=coarse_patch.PATCH_SCALE_A=ie.PATCH_SCALE_A=p
    core.SC_RSA_THRESHOLD=ie.SC_RSA_THRESHOLD=t
def run_sens(row):
    pid=row['pdb_id'].upper()
    prep.prepare(pid)
    cif=prep.OUT/f'{pid}-assembly1.cif'
    assert sha(cif.read_bytes())==row['assembly_sha256']==provenance[pid]['source_sha256'],'assembly hash mismatch'
    text=(prep.OUT/f'{pid}_prepared.pdb').read_text()
    assert sha(text.encode())==provenance[pid]['prepared_sha256'],'prepared hash mismatch'
    result=[];base=None
    for label,c in CONFIG:
        configure(c)
        q=ie.analyze_all_maps(pdb_text=text,pH=c[-1],chain=None)
        assert q['version']==ms.MODEL_VERSION and len(q['maps'])==11
        surface={x['key'] for x in q['residue_profiles'] if x['surface_exposed']}
        if base is None:
            base=q;surface0=surface
            for ch,m in q['maps'].items():
                old=oldmaps[(pid,ch)]
                assert members(m,5)==set(filter(None,old['top5_members'].split(';'))),'archived baseline members mismatch'
                assert centers(m,5)==set(filter(None,old['top5_centers'].split(';'))),'archived baseline centers mismatch'
                assert len(m['patches'])==int(old['n_patches']),'archived baseline patch count mismatch'
        for ch,m in q['maps'].items():
            orig=base['maps'][ch]
            d=dict(pdb=pid,channel=ch,configuration=label,prepared_sha256=sha(text.encode()),sasa_points=c[0],inner=c[1][0],outer=c[1][1],patch_radius=c[2],trsa_gate=c[3],pH=c[4],n_patches=len(m['patches']),baseline_n_patches=len(orig['patches']),n_exposed=len(surface),exposed_jaccard=jac(surface,surface0))
            for k in (1,3,5):
                x,y=members(orig,k),members(m,k);xc,yc=centers(orig,k),centers(m,k)
                d.update({f'jaccard_top{k}':jac(x,y),f'center_jaccard_top{k}':jac(xc,yc),f'baseline_n_top{k}':len(x),f'variant_n_top{k}':len(y),f'empty_union_top{k}':int(not(x|y))})
            result.append(d)
    assert len(result)==242
    return result
def run_pb(row):
    from electrostatics import _run,_pqr_atoms,_apbs_input,_read_dx,_pdb_orphan_n_only_residues,_RESIDUE_ALIASES,ElectrostaticsError
    from surface_potential import sample_residue_sas_potential
    pid=row['proposed_pdb_id'];wd=out/pid;wd.mkdir(exist_ok=True)
    meta,path=prepare(row);text=path.read_text();old=provenance[pid]
    assert meta['assembly_cif_sha256']==old['prepared_source_sha256'] and sha(text.encode())==old['prepared_PDB_sha256'],'assembly or prepared hash mismatch'
    context=ie._prepare_shared_context(pH=7.4,ionic_mM=150.,temp_K=298.,pdb_id=None,pdb_text=text,chain=None)
    profiles=ie.profile_rows(context['core_result']);text=context['prepared']
    valid={r['key']:r['res_name'] for r in profiles}
    (wd/'structure.pdb').write_text(text);omitted=[]
    cmd=[shutil.which('pdb2pqr'),'--ff=AMBER','--keep-chain','--titration-state-method=propka','--whitespace','--with-ph=7.4000','structure.pdb','structure.pqr']
    try:log=_run(cmd,wd)
    except ElectrostaticsError as e:
        (wd/'pqr_initial_failure.txt').write_text(str(e));cleaned,omitted=_pdb_orphan_n_only_residues(text,valid)
        if not omitted:raise
        (wd/'structure.pdb').write_text(cleaned);(wd/'structure.pqr').unlink(missing_ok=True);log=_run(cmd,wd)
    (wd/'pqr.log').write_text(log);atoms=_pqr_atoms(wd/'structure.pqr')
    mapped={f'{x["chain"]}:{x["seq"]}:{x["icode"]}' for x in atoms if valid.get(f'{x["chain"]}:{x["seq"]}:{x["icode"]}')==_RESIDUE_ALIASES.get(x['res_name'],x['res_name'])}
    assert mapped==set(valid)-{x['key'] for x in omitted},'mapping mismatch'
    exposed={r['key'] for r in profiles if r['surface_exposed'] and r['key'] in mapped}
    inp=_apbs_input(atoms,150.,298.);assert inp.count('    lpbe\n')==1
    vals={};grids={};samples={}
    for typ in ('lpbe','npbe'):
        for f in wd.glob('potential*.dx'):f.unlink()
        (wd/f'{typ}.in').write_text(inp.replace('    lpbe\n',f'    {typ}\n'))
        log=_run([shutil.which('apbs'),f'{typ}.in'],wd,timeout=2400);(wd/f'{typ}.log').write_text(log)
        files=list(wd.glob('potential*.dx'));assert len(files)==1,'missing or ambiguous DX'
        grid,origin,spacing=_read_dx(files[0]);grids[typ]=(grid.shape,origin.tolist(),spacing.tolist())
        assert np.isfinite(grid).all(),'nonfinite grid'
        d,sm=sample_residue_sas_potential(atoms,grid,origin,spacing,exposed)
        vals[typ]={k:v['surface_potential_kT_per_e'] for k,v in d.items() if v.get('surface_potential_kT_per_e') is not None}
        samples[typ]={k:v.get('surface_potential_samples') for k,v in d.items()}
        dump(wd/f'{typ}_sampling.json',sm)
        files[0].unlink()
    assert grids['lpbe']==grids['npbe'],'grid mismatch'
    assert samples['lpbe']==samples['npbe'],'sampling mismatch'
    keys=sorted(vals['lpbe']);assert keys and set(keys)==set(vals['npbe'])
    x=np.array([vals['lpbe'][k] for k in keys]);y=np.array([vals['npbe'][k] for k in keys])
    assert np.isfinite(x).all() and np.isfinite(y).all()
    oldpath=list(ref.rglob(f'proteins/{pid}/real_APBS_SAS_per_residue.csv'))
    archived={}
    if len(oldpath)==1:
        for z in csvread(oldpath[0]):
            v=z.get('surface_potential_kT_per_e','')
            if v not in ('','None','nan'):archived[z['key']]=float(v)
    common=sorted(set(keys)&set(archived))
    olddelta=[abs(vals['lpbe'][k]-archived[k]) for k in common]
    savecsv(wd/'paired.csv',[dict(key=k,LPBE_kT_e=float(u),NPBE_kT_e=float(v),delta_NPBE_minus_LPBE=float(v-u),archived_LPBE_kT_e=archived.get(k,''),shell_samples=samples['lpbe'][k]) for k,u,v in zip(keys,x,y)])
    rho=float(spearmanr(x,y).statistic)
    q=dict(pdb=pid,status='paired_solver_outputs_review_required',pqr_sha256=sha((wd/'structure.pqr').read_bytes()),source_sha256=meta['assembly_cif_sha256'],prepared_sha256=old['prepared_PDB_sha256'],grid=grids['lpbe'],omitted=omitted,n_exposed=len(exposed),n_paired=len(keys),median_lpbe=float(np.median(x)),median_npbe=float(np.median(y)),median_abs_delta=float(np.median(abs(x-y))),p95_abs_delta=float(np.quantile(abs(x-y),.95)),sign_disagreement=float(np.mean(np.sign(x)!=np.sign(y))),spearman=rho if math.isfinite(rho) else None,n_lpbe_abs_gt10=int(sum(abs(x)>10)),n_npbe_abs_gt10=int(sum(abs(y)>10)),n_archived_LPBE_overlap=len(common),archived_LPBE_max_abs_difference=max(olddelta) if olddelta else None)
    dump(wd/'summary.json',q);return q
if mode=='sens':
    import sensitivity_48 as original,strict_apo_compat_48 as prep
    panel=original.MANIFEST;assert len(panel)==48
    provenance={x['protein']:x for x in json.loads(unique_file('provenance.json').read_text())}
    oldmaps={(x['protein'],x['channel']):x for x in csvread(unique_file('all_528_maps.csv'))}
    assert len(provenance)==48 and len(oldmaps)==528
else:
    from isv_full_80_qa import prepare
    panel=csvread(ROOT/'validation/ISV_active_candidates_76_v3.csv');assert len(panel)==76
    provenance={x['pdb_id']:x for x in csvread(unique_file('ISV76_ALL_STRUCTURE_STATUS.csv'))};assert len(provenance)==76
rows=[];failures=[]
for i,row in enumerate(panel):
    if i%shards!=shard:continue
    pid=row.get('pdb_id',row.get('proposed_pdb_id'))
    try:
        if mode=='sens':rows.extend(run_sens(row))
        else:rows.append(run_pb(row))
        print('COMPUTED',pid,flush=True)
    except Exception as e:
        failures.append(dict(pdb=pid,error=str(e),trace=traceback.format_exc()))
        print('FAILED',pid,str(e),flush=True)
    if mode=='sens':savecsv(out/f'sensitivity_{shard}.csv',rows)
    else:dump(out/f'paired_{shard}.json',rows)
    dump(out/f'failures_{shard}.json',failures)
sys.exit(bool(failures))
