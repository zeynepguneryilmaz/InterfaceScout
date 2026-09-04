#!/usr/bin/env python3
from __future__ import annotations
import json, shutil, sys, tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BACKEND=ROOT/'backend'
sys.path.insert(0,str(BACKEND))

import main as core
from structural_context import AnalyzeRequest, EnvParams, prepare_context
from structure_repair import repair_for_forcefield
from nonpolar_explicit import scan_explicit

core.PDB2PQR=None; core.APBS=None; core.MKDSSP=None
PDB='1PGA'; CHAIN='A'
MODE_A=set(range(13,20)); MODE_B=set(range(42,47)); HELIX=set(range(23,37))


def ids(o):
    return {int(r['res_seq']) for r in o.get('contact_residues',[]) if r.get('chain')==CHAIN}

def rec(s,t): return len(s&t)/len(t)

def summarize(o):
    s=ids(o)
    return {
        'rank':o['rank'],'orientation_index':o['orientation_index'],
        'energy_kj_mol':o['total_energy_change_kj_mol'],
        'vdw_kj_mol':o['vdw_energy_kj_mol'],
        'solvation_kj_mol':o['solvation_energy_change_kj_mol'],
        'minimum_separation_A':o['minimum_separation_A'],
        'mode_A_recall':rec(s,MODE_A),'mode_B_recall':rec(s,MODE_B),'helix_recall':rec(s,HELIX),
        'contact_residues':sorted(s),'n_contact_residues':len(s),
        'protein_buried_sasa_A2':o['protein_buried_sasa_A2'],'graphene_buried_sasa_A2':o['graphene_buried_sasa_A2'],
    }

def run(n):
    work=Path(tempfile.mkdtemp(prefix='is_harrison_explicit_'))
    try:
        req=AnalyzeRequest(pdb_id=PDB,chain=CHAIN,structure_context='selected_chain_legacy',protrusion=False,
                           env=EnvParams(pH=7.0,ionic=150.0,temp=298.0))
        pdb,_,context,_=prepare_context(req,work)
        fixed=work/'forcefield_ready.pdb'
        repair=repair_for_forcefield(pdb,fixed)
        if repair.get('status')!='ok': raise RuntimeError(repair)
        struct,_,_,_=core.build_surface_residues(fixed,7.0)
        result=scan_explicit(fixed,struct,pH=7.0,n_orientations=n,separations_A=(2.8,3.2,3.6,4.0,4.4),sasa_points=100)
        if result.get('status')!='ok': raise RuntimeError(result)
        rows=[summarize(o) for o in result['top_orientations']]
        bestA=max(rows,key=lambda r:(r['mode_A_recall'],r['helix_recall'],-r['rank']))
        bestB=max(rows,key=lambda r:(r['mode_B_recall'],r['helix_recall'],-r['rank']))
        firstA=next((r for r in rows if r['mode_A_recall']>r['mode_B_recall']),None)
        firstB=next((r for r in rows if r['mode_B_recall']>r['mode_A_recall']),None)
        union=set()
        for r in rows[:10]: union.update(r['contact_residues'])
        return {
            'n_orientations':n,'structure_context':context,'model':result['method'],
            'best_energy_kj_mol':result['best_energy_change_kj_mol'],'graphene':result['graphene'],'sasa':result['sasa'],
            'top1':rows[0],'best_mode_A_top20':bestA,'best_mode_B_top20':bestB,
            'first_A_dominant_top20':firstA,'first_B_dominant_top20':firstB,
            'two_published_faces_present_top20':bool(firstA and firstB),
            'union_top10_A_recall':rec(union,MODE_A),'union_top10_B_recall':rec(union,MODE_B),'union_top10_helix_recall':rec(union,HELIX),
            'top20':rows,'structure_repair':repair,
        }
    finally:
        shutil.rmtree(work,ignore_errors=True)

def main():
    r128=run(128)
    r256=run(256)
    report={
        'status':'COMPLETED',
        'classification':'literature mechanistic reproduction / consistency validation; not independent held-out validation',
        'source':{'citation':'Harrison ET, Weidner T, Castner DG, Interlandi G. Biointerphases 2017;12(2):02D401.',
                  'doi':'10.1116/1.4971381','pdb':'1PGA','mode_A':list(range(13,20)),'mode_B':list(range(42,47)),'helix':list(range(23,37)),'contact_cutoff_A':6.0},
        'frozen_before_run':True,
        'no_parameter_tuning_after_literature_targets':True,
        'grid_128':r128,'grid_256':r256,
    }
    out=ROOT/'validation'/'literature_harrison_1pga_explicit_report.json'
    out.write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
