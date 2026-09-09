"""Independent, label-free selection of only the two numerical/model parameters that require empirical choice.

Selected here:
1) Shrake-Rupley sampling density (points/atom) by numerical convergence to 1000 points/atom.
2) Multiscale spatial radius pair by local-perturbation robustness on an adsorption-label-free structural panel.

Not selected here: probe radius, scRSA threshold, pKa values, same-face gate, pH, ionic strength, temperature.
Those are fixed by literature convention, chemical reference values, or parameter-free method definition.
"""
from __future__ import annotations

import json
import math
import urllib.request
from pathlib import Path
from statistics import mean, median

import numpy as np

import main as v1
from v2.chemistry_freeze import apply_publication_chemistry
from v2.prepare import prepare_pdb_text

# Independent development panel: none of these proteins occurs in the planned
# experimental interface-localization benchmark.
STRUCTURES = [
    {"id":"crambin","label":"Crambin","pdb_id":"1CRN","chain":"A","class":"very small compact"},
    {"id":"repressor434","label":"Phage 434 repressor N-terminal domain","pdb_id":"1R69","chain":"A","class":"small all-alpha"},
    {"id":"sh3","label":"Alpha-spectrin SH3 domain","pdb_id":"1SHG","chain":"A","class":"small all-beta"},
    {"id":"fkbp12","label":"FKBP12","pdb_id":"1FKK","chain":"A","class":"small alpha-beta"},
    {"id":"adenylate_kinase","label":"Adenylate kinase","pdb_id":"1AKE","chain":"A","class":"medium alpha-beta enzyme"},
    {"id":"maltose_binding_protein","label":"Maltodextrin-binding protein","pdb_id":"1OMP","chain":"A","class":"large multidomain"},
    {"id":"tim_dimer","label":"Triosephosphate isomerase dimer","pdb_id":"1TIM","chain":"A,B","class":"dimeric enzyme"},
    {"id":"citrate_synthase","label":"Citrate synthase","pdb_id":"2CTS","chain":"A","class":"large enzyme"},
]

SASA_POINTS = [100, 200, 500, 1000]
SASA_REFERENCE = 1000

# Pre-specified numerical convergence tolerances. These are not adsorption-performance
# thresholds and do not use experimental interface labels.
SASA_CRITERIA = {
    "median_surface_jaccard_min": 0.98,
    "worst_surface_jaccard_min": 0.95,
    "median_scrsa_mae_max": 0.01,
    "median_top10_hotspot_jaccard_min": 0.90,
}

# Physically local C-alpha neighborhood range only. Candidate pairs are judged without
# designating any one pair as a reference.
RADIUS_PAIRS = [
    (5.0,7.0),(5.0,8.0),(5.0,9.0),(5.0,10.0),
    (6.0,8.0),(6.0,9.0),(6.0,10.0),
    (7.0,9.0),(7.0,10.0),
    (8.0,10.0),
]


def fetch_pdb(pid: str) -> str:
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.pdb", timeout=60) as r:
        return r.read().decode("utf-8", errors="replace")


def run_case(raw: str, chain: str, points: int) -> tuple[dict,str,dict]:
    prepared, prep = prepare_pdb_text(raw, chain=chain)
    old_points = v1.SASA_POINTS
    try:
        v1.SASA_POINTS = int(points)
        req = v1.AnalyzeRequest(pdb_text=prepared, chain=None, env=v1.EnvParams(pH=7.4, ionic=150.0, temp=298.0))
        out = v1.analyze(req)
    finally:
        v1.SASA_POINTS = old_points
    return out, prepared, prep


def jaccard(a: set[str], b: set[str]) -> float:
    u = a | b
    return len(a & b) / len(u) if u else 1.0


def surface_set(out: dict) -> set[str]:
    return {str(r['key']) for r in out.get('surface_residues', []) if r.get('key')}


def scrsa_map(out: dict) -> dict[str,float]:
    d = {}
    for r in out.get('all_residues', []):
        if not r.get('key'): continue
        val = r.get('scrsa_raw', r.get('scrsa'))
        if val is not None: d[str(r['key'])] = float(val)
    return d


def top_sets(out: dict, k: int=10) -> dict[str,set[str]]:
    d = {}
    for c,payload in out.get('chemistries', {}).items():
        d[c] = {str(x['center_key']) for x in payload.get('top_patches', [])[:k] if x.get('center_key')}
    return d


def median_informative_jaccard(a: dict[str,set[str]], b: dict[str,set[str]]) -> float | None:
    vals = []
    for c in sorted(set(a)|set(b)):
        aa,bb = a.get(c,set()),b.get(c,set())
        if not (aa|bb):
            continue
        vals.append(jaccard(aa,bb))
    return float(median(vals)) if vals else None


def percentile(vals, q):
    if not vals: return None
    return float(np.percentile(np.asarray(vals,dtype=float), q))


def compute_pair_topsets(surface: list[dict], distances: np.ndarray, pair: tuple[float,float], k: int) -> dict[str,set[str]]:
    result = {}
    for c,meta in v1.CHEMISTRIES.items():
        fav = meta['favorable']
        local = np.zeros(len(surface), dtype=float)
        for i,r in enumerate(surface):
            rn = r['res_name']
            if rn in fav:
                _,_,state = fav[rn]
                local[i] = float(r['scrsa']) * float(v1.state_availability(rn,state,7.4))
        if len(surface)==0 or float(np.max(local))<=0:
            result[c]=set(); continue
        norms=[]
        for R in pair:
            d = (distances <= R).astype(float) @ local
            mx = float(np.max(d)) if d.size else 0.0
            norms.append(d/mx if mx>0 else np.zeros_like(d))
        persistence = 100.0*np.minimum(norms[0], norms[1])
        geomean = 100.0*np.sqrt(norms[0]*norms[1])
        rows=[]
        for i,r in enumerate(surface):
            if persistence[i] <= 0: continue
            rows.append((float(persistence[i]),float(geomean[i]),str(r['key'])))
        rows.sort(key=lambda x:(-x[0],-x[1],x[2]))
        result[c] = {x[2] for x in rows[:k]}
    return result


def adjacent_pairs(pair):
    r1,r2 = pair
    out=[]
    for dr1,dr2 in [(-1,0),(1,0),(0,-1),(0,1)]:
        p=(r1+dr1,r2+dr2)
        if p in RADIUS_PAIRS:
            out.append(p)
    return out


def main():
    apply_publication_chemistry(v1)
    Path('selection_pdbs/full').mkdir(parents=True, exist_ok=True)
    Path('selection_pdbs/selected').mkdir(parents=True, exist_ok=True)

    raw_by_id={}
    prep_by_id={}
    outputs={}
    pdb_manifest=[]

    # Download once and retain exact structures used.
    for s in STRUCTURES:
        raw=fetch_pdb(s['pdb_id']); raw_by_id[s['id']]=raw
        Path(f"selection_pdbs/full/{s['pdb_id']}.pdb").write_text(raw,encoding='utf-8')
        # 100-point run supplies prepared PDB; preparation does not depend on sampling density.
        out,prepared,prep=run_case(raw,s['chain'],100)
        prep_by_id[s['id']]=prepared
        Path(f"selection_pdbs/selected/{s['id']}_{s['pdb_id']}_chains_{s['chain'].replace(',','-')}.pdb").write_text(prepared,encoding='utf-8')
        outputs[(s['id'],100)] = out
        pdb_manifest.append({**s, **prep})

    # SASA convergence against 1000-point numerical reference.
    for s in STRUCTURES:
        for pts in SASA_POINTS:
            if (s['id'],pts) not in outputs:
                outputs[(s['id'],pts)] = run_case(raw_by_id[s['id']], s['chain'], pts)[0]

    sasa_rows=[]
    for s in STRUCTURES:
        ref=outputs[(s['id'],SASA_REFERENCE)]
        ref_surf=surface_set(ref); ref_scrsa=scrsa_map(ref); ref_top=top_sets(ref,10)
        for pts in SASA_POINTS:
            cur=outputs[(s['id'],pts)]
            cur_surf=surface_set(cur); cur_scrsa=scrsa_map(cur); cur_top=top_sets(cur,10)
            keys=sorted(set(ref_scrsa)&set(cur_scrsa))
            dif=np.asarray([cur_scrsa[k]-ref_scrsa[k] for k in keys],dtype=float) if keys else np.asarray([],dtype=float)
            sasa_rows.append({
                **s,'points_per_atom':pts,'reference_points':SASA_REFERENCE,
                'n_residues_compared':len(keys),'n_surface':len(cur_surf),'n_surface_reference':len(ref_surf),
                'surface_jaccard_vs_1000':jaccard(cur_surf,ref_surf),
                'surface_crossing_count_vs_1000':len(cur_surf ^ ref_surf),
                'scrsa_mae_vs_1000':float(np.mean(np.abs(dif))) if dif.size else None,
                'scrsa_rmse_vs_1000':float(math.sqrt(np.mean(dif*dif))) if dif.size else None,
                'top10_hotspot_jaccard_vs_1000':median_informative_jaccard(cur_top,ref_top),
            })

    sasa_summary=[]
    selected_points=None
    for pts in SASA_POINTS:
        rows=[r for r in sasa_rows if r['points_per_atom']==pts]
        sj=[r['surface_jaccard_vs_1000'] for r in rows]
        ma=[r['scrsa_mae_vs_1000'] for r in rows]
        hj=[r['top10_hotspot_jaccard_vs_1000'] for r in rows if r['top10_hotspot_jaccard_vs_1000'] is not None]
        summary={
            'points_per_atom':pts,
            'median_surface_jaccard':float(median(sj)),
            'mean_surface_jaccard':float(mean(sj)),
            'std_surface_jaccard':float(np.std(sj,ddof=1)),
            'worst_surface_jaccard':float(min(sj)),
            'median_scrsa_mae':float(median(ma)),
            'mean_scrsa_mae':float(mean(ma)),
            'std_scrsa_mae':float(np.std(ma,ddof=1)),
            'median_top10_hotspot_jaccard':float(median(hj)),
            'mean_top10_hotspot_jaccard':float(mean(hj)),
            'std_top10_hotspot_jaccard':float(np.std(hj,ddof=1)),
        }
        summary['passes_predefined_convergence'] = (
            summary['median_surface_jaccard'] >= SASA_CRITERIA['median_surface_jaccard_min'] and
            summary['worst_surface_jaccard'] >= SASA_CRITERIA['worst_surface_jaccard_min'] and
            summary['median_scrsa_mae'] <= SASA_CRITERIA['median_scrsa_mae_max'] and
            summary['median_top10_hotspot_jaccard'] >= SASA_CRITERIA['median_top10_hotspot_jaccard_min']
        )
        sasa_summary.append(summary)
        if selected_points is None and summary['passes_predefined_convergence']:
            selected_points=pts
    if selected_points is None:
        selected_points=SASA_REFERENCE

    # Spatial scale selection at selected SASA density. No pair is treated as reference.
    radius_case_data={}
    for s in STRUCTURES:
        out=outputs.get((s['id'],selected_points)) or run_case(raw_by_id[s['id']],s['chain'],selected_points)[0]
        surface=out.get('surface_residues',[])
        coords=np.asarray([[r['x'],r['y'],r['z']] for r in surface],dtype=float)
        dif=coords[:,None,:]-coords[None,:,:] if len(coords) else np.empty((0,0,3))
        dist=np.sqrt(np.sum(dif*dif,axis=2)) if len(coords) else np.empty((0,0))
        radius_case_data[s['id']]={'surface':surface,'dist':dist,'top5':{},'top10':{}}
        for pair in RADIUS_PAIRS:
            radius_case_data[s['id']]['top5'][pair]=compute_pair_topsets(surface,dist,pair,5)
            radius_case_data[s['id']]['top10'][pair]=compute_pair_topsets(surface,dist,pair,10)

    radius_rows=[]
    for pair in RADIUS_PAIRS:
        neighbors=adjacent_pairs(pair)
        top5_vals=[]; top10_vals=[]; locality_vals=[]; neighbor_counts=[]
        per_protein=[]
        for s in STRUCTURES:
            data=radius_case_data[s['id']]; surface=data['surface']; dist=data['dist']
            p5=[]; p10=[]
            for q in neighbors:
                for c in sorted(v1.CHEMISTRIES):
                    a5,b5=data['top5'][pair].get(c,set()),data['top5'][q].get(c,set())
                    if a5|b5: p5.append(jaccard(a5,b5))
                    a10,b10=data['top10'][pair].get(c,set()),data['top10'][q].get(c,set())
                    if a10|b10: p10.append(jaccard(a10,b10))
            top5_vals.extend(p5); top10_vals.extend(p10)
            if len(surface):
                counts=np.sum(dist <= pair[1],axis=1)
                fracs=counts/len(surface)
                locality_vals.extend([float(x) for x in fracs])
                neighbor_counts.extend([float(x) for x in counts])
            per_protein.append({
                **s,'radius_inner_A':pair[0],'radius_outer_A':pair[1],
                'n_adjacent_pairs':len(neighbors),
                'median_local_top5_jaccard':float(median(p5)) if p5 else None,
                'median_local_top10_jaccard':float(median(p10)) if p10 else None,
                'median_outer_neighborhood_fraction':float(median(fracs)) if len(surface) else None,
                'median_outer_neighbor_count':float(median(counts)) if len(surface) else None,
            })
        row={
            'radius_inner_A':pair[0],'radius_outer_A':pair[1],
            'adjacent_pairs':[list(x) for x in neighbors],
            'n_informative_top5_comparisons':len(top5_vals),
            'top5_local_robustness_q10':percentile(top5_vals,10),
            'top5_local_robustness_median':float(median(top5_vals)) if top5_vals else None,
            'top5_local_robustness_mean':float(mean(top5_vals)) if top5_vals else None,
            'top5_local_robustness_std':float(np.std(top5_vals,ddof=1)) if len(top5_vals)>1 else 0.0,
            'top10_local_robustness_q10':percentile(top10_vals,10),
            'top10_local_robustness_median':float(median(top10_vals)) if top10_vals else None,
            'top10_local_robustness_mean':float(mean(top10_vals)) if top10_vals else None,
            'top10_local_robustness_std':float(np.std(top10_vals,ddof=1)) if len(top10_vals)>1 else 0.0,
            'median_outer_neighborhood_fraction':float(median(locality_vals)) if locality_vals else None,
            'median_outer_neighbor_count':float(median(neighbor_counts)) if neighbor_counts else None,
            'per_protein':per_protein,
        }
        radius_rows.append(row)

    # Objective, label-free lexicographic selection: prioritize lower-tail stability,
    # then median stability, then locality. No candidate is privileged as a reference.
    ranked=sorted(radius_rows,key=lambda r:(
        -float(r['top5_local_robustness_q10']),
        -float(r['top5_local_robustness_median']),
        float(r['median_outer_neighborhood_fraction']),
        float(r['radius_outer_A']),
        float(r['radius_inner_A']),
    ))
    for i,r in enumerate(ranked,1): r['selection_rank']=i
    selected_pair=[ranked[0]['radius_inner_A'],ranked[0]['radius_outer_A']]

    payload={
        'purpose':'independent label-free core parameter selection',
        'development_panel':STRUCTURES,
        'fixed_environment_for_selection':{'pH':7.4,'ionic_mM':150.0,'temp_K':298.0,'scrsa_threshold':v1.SC_RSA_THRESHOLD,'probe_A':v1.SASA_PROBE_A},
        'sasa_selection':{
            'tested_points_per_atom':SASA_POINTS,'reference_points_per_atom':SASA_REFERENCE,
            'predefined_convergence_criteria':SASA_CRITERIA,
            'per_protein':sasa_rows,'summary':sasa_summary,'selected_points_per_atom':selected_points,
        },
        'radius_selection':{
            'tested_pairs_A':[list(x) for x in RADIUS_PAIRS],
            'selection_rule':'lexicographic: maximize 10th-percentile top-5 local-perturbation Jaccard; maximize median top-5 Jaccard; minimize median outer-neighborhood fraction; minimize outer then inner radius',
            'note':'empty-empty chemistry comparisons are excluded; no experimental adsorption labels and no designated reference radius pair are used',
            'candidates':ranked,'selected_pair_A':selected_pair,
        },
        'pdb_manifest':pdb_manifest,
    }
    with open('core_parameter_selection.json','w',encoding='utf-8') as f:
        json.dump(payload,f,indent=2,sort_keys=True)
    print(json.dumps({'selected_points_per_atom':selected_points,'selected_pair_A':selected_pair,'sasa_summary':sasa_summary,'radius_ranking':[{k:r[k] for k in ['selection_rank','radius_inner_A','radius_outer_A','top5_local_robustness_q10','top5_local_robustness_median','median_outer_neighborhood_fraction']} for r in ranked]},indent=2))

if __name__=='__main__':
    main()
