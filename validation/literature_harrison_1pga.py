#!/usr/bin/env python3
"""Literature-consistency validation of InterfaceScout 2.0 nonpolar physics.

This validation is intentionally NOT called held-out validation because the
Harrison et al. model motivated the v2 nonpolar physics architecture. It asks a
more limited question: does the lightweight InterfaceScout adaptation recover
the two published protein-G B1 contact modes when run from the published 1PGA
structure without tuning to those modes?

Predeclared literature targets (Harrison et al., Biointerphases 2017,
DOI 10.1116/1.4971381):
  mode A hallmark beta strand: residues 13-19
  mode B hallmark beta strand: residues 42-46
  shared alpha helix: residues 23-36
  contact definition in source study: any atom within 6 A of graphene

Primary outputs are continuous/descriptive. The script does not alter model
parameters or declare a scientific PASS based on an arbitrary recovery cutoff.
"""
from __future__ import annotations

import json, shutil, sys, tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BACKEND=ROOT/'backend'
sys.path.insert(0,str(BACKEND))

import main as core
from structural_context import AnalyzeRequest, EnvParams, prepare_context
from structure_repair import repair_for_forcefield
from nonpolar_energy import scan

core.PDB2PQR=None; core.APBS=None; core.MKDSSP=None

PDB='1PGA'
CHAIN='A'
MODE_A=set(range(13,20))
MODE_B=set(range(42,47))
HELIX=set(range(23,37))


def residue_ids(orientation):
    return {int(r['res_seq']) for r in orientation.get('contact_residues',[]) if r.get('chain')==CHAIN}


def recall(ids, target):
    return len(ids & target)/len(target)


def summarize_orientation(o):
    ids=residue_ids(o)
    ra=recall(ids,MODE_A); rb=recall(ids,MODE_B); rh=recall(ids,HELIX)
    return {
        'rank':o['rank'],
        'orientation_index':o['orientation_index'],
        'energy_kj_mol':o['total_energy_change_kj_mol'],
        'vdw_kj_mol':o['vdw_energy_kj_mol'],
        'solvation_kj_mol':o['solvation_energy_change_kj_mol'],
        'minimum_separation_A':o['minimum_separation_A'],
        'mode_A_13_19_recall':ra,
        'mode_B_42_46_recall':rb,
        'helix_23_36_recall':rh,
        'mode_A_minus_B':ra-rb,
        'contact_residues':sorted(ids),
        'n_contact_residues':len(ids),
    }


def run(n_orientations):
    work=Path(tempfile.mkdtemp(prefix='is_harrison_1pga_'))
    try:
        req=AnalyzeRequest(pdb_id=PDB,chain=CHAIN,structure_context='selected_chain_legacy',protrusion=False,
                           env=EnvParams(pH=7.0,ionic=150.0,temp=298.0))
        pdb,_,context,_=prepare_context(req,work)
        fixed=work/'forcefield_ready.pdb'
        repair=repair_for_forcefield(pdb,fixed)
        if repair.get('status')!='ok': raise RuntimeError(repair)
        struct,_,_,_=core.build_surface_residues(fixed,7.0)
        result=scan(fixed,struct,pH=7.0,n_orientations=n_orientations)
        if result.get('status')!='ok': raise RuntimeError(result)
        rows=[summarize_orientation(o) for o in result['top_orientations']]
        # First/best-ranked representatives of each published hallmark face.
        best_A=max(rows,key=lambda r:(r['mode_A_13_19_recall'],r['helix_23_36_recall'],-r['rank']))
        best_B=max(rows,key=lambda r:(r['mode_B_42_46_recall'],r['helix_23_36_recall'],-r['rank']))
        first_A_dom=next((r for r in rows if r['mode_A_13_19_recall']>r['mode_B_42_46_recall']),None)
        first_B_dom=next((r for r in rows if r['mode_B_42_46_recall']>r['mode_A_13_19_recall']),None)
        union=set()
        for r in rows[:10]: union.update(r['contact_residues'])
        return {
            'n_orientations':n_orientations,
            'structure_context':context,
            'contact_cutoff_A':result['contact_distance_A'],
            'best_energy_kj_mol':result['best_energy_change_kj_mol'],
            'lj_parameter_source':result['parameter_source'],
            'lj_parameter_diagnostics':result['lj_parameter_diagnostics'],
            'top1':rows[0],
            'best_mode_A_within_top20':best_A,
            'best_mode_B_within_top20':best_B,
            'first_mode_A_dominant_within_top20':first_A_dom,
            'first_mode_B_dominant_within_top20':first_B_dom,
            'two_distinct_hallmark_directions_present_top20':bool(first_A_dom and first_B_dom),
            'union_top10_mode_A_recall':recall(union,MODE_A),
            'union_top10_mode_B_recall':recall(union,MODE_B),
            'union_top10_helix_recall':recall(union,HELIX),
            'top20':rows,
            'structure_repair':repair,
        }
    finally:
        shutil.rmtree(work,ignore_errors=True)


def main():
    # Grid convergence is assessed without changing any energy parameter.
    r512=run(512)
    r1024=run(1024)
    def signature(r):
        return {
            'A_best_recall':r['best_mode_A_within_top20']['mode_A_13_19_recall'],
            'B_best_recall':r['best_mode_B_within_top20']['mode_B_42_46_recall'],
            'A_union10':r['union_top10_mode_A_recall'],
            'B_union10':r['union_top10_mode_B_recall'],
            'two_modes':r['two_distinct_hallmark_directions_present_top20'],
        }
    s512=signature(r512); s1024=signature(r1024)
    report={
        'status':'COMPLETED',
        'classification':'literature mechanistic reproduction / consistency validation; NOT independent held-out validation',
        'model_version':'2.0.0-dev',
        'frozen_model_sha':'220d38be3027f124ff19b7a23cf27b04d1c29374',
        'source':{
            'citation':'Harrison ET, Weidner T, Castner DG, Interlandi G. Biointerphases. 2017;12(2):02D401.',
            'doi':'10.1116/1.4971381',
            'pdb':'1PGA',
            'published_mode_A_hallmark':[13,14,15,16,17,18,19],
            'published_mode_B_hallmark':[42,43,44,45,46],
            'published_shared_helix':list(range(23,37)),
            'published_contact_cutoff_A':6.0,
        },
        'predeclared_interpretation':{
            'primary_question':'Does the low-energy orientation ensemble contain orientations dominated by each of the two published hallmark beta-strand faces?',
            'secondary':'Report continuous hallmark and helix recall, energy rank, and 512/1024 grid stability; no post-hoc score tuning.',
            'scientific_pass_threshold':'none predeclared; results are reported continuously to avoid arbitrary success criteria',
        },
        'grid_512':r512,
        'grid_1024':r1024,
        'grid_signature_512':s512,
        'grid_signature_1024':s1024,
        'grid_signature_identical':s512==s1024,
    }
    out=ROOT/'validation'/'literature_harrison_1pga_report.json'
    out.write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
