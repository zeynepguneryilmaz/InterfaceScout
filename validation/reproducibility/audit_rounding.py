"""Replay the pre-P0 kernel to attribute archived 209-map changes to rounding.

Requires repository history and the original numeric-before/numeric-after captures.
These archives are comparison evidence, never scientific reference oracles.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import types
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'backend'),str(ROOT/'validation')]
import core
import interface_engine as engine
from reproduce import checked_input

BASE='3344970751fc3778f81878d6c9d58ef6cf6fc9de'
FIELDS=('chemistry_support','mean_accessibility','patch_coherence','orientation_coherence')


def load_legacy():
    source=subprocess.check_output(['git','show',BASE+':backend/core.py'],cwd=ROOT).decode('utf-8')
    class RemoveRound(ast.NodeTransformer):
        count=0
        def visit_Call(self,node):
            node=self.generic_visit(node)
            if isinstance(node.func,ast.Name) and node.func.id=='round':
                self.count+=1
                return node.args[0]
            return node
    transform=RemoveRound()
    old=transform.visit(ast.parse(source)); new=ast.parse((ROOT/'backend/core.py').read_text(encoding='utf-8'))
    old.body=old.body[1:]; new.body=new.body[1:]  # only module documentation excluded
    assert transform.count==20
    assert ast.dump(old)==ast.dump(new),'Kernel differs beyond the 20 removed rounding calls'
    for file in ('coarse_patch.py','geometry.py','model_settings.py'):
        original=subprocess.check_output(['git','show',BASE+':backend/'+file],cwd=ROOT)
        assert original.replace(b'\r\n',b'\n')==(ROOT/'backend'/file).read_bytes().replace(b'\r\n',b'\n'),file
    legacy=types.ModuleType('legacy_rounding_kernel'); sys.modules[legacy.__name__]=legacy
    exec(compile(source,'historical/backend/core.py','exec'),legacy.__dict__)
    return legacy,hashlib.sha256(source.encode()).hexdigest()


def dominates(a,b):
    return all(a[f]>=b[f] for f in FIELDS) and any(a[f]>b[f] for f in FIELDS)


def details(before,after):
    old={p['center_key']:p for p in before['patches']}; new={p['center_key']:p for p in after['patches']}
    common=old.keys() & new.keys()
    changed=[]
    for key in sorted(old.keys()|new.keys()):
        a,b=old.get(key),new.get(key)
        if a is None or b is None or (a['display_rank'],a['pareto_front'])!=(b['display_rank'],b['pareto_front']):
            def compact(p):
                return {f:p[f] for f in ('display_rank','pareto_front','members','seed_members',*FIELDS)} if p else None
            changed.append(dict(center=key,before=compact(a),after=compact(b)))
    edges=[]
    for a in sorted(common):
        for b in sorted(common):
            if a!=b and dominates(old[a],old[b])!=dominates(new[a],new[b]):
                edges.append(dict(dominator=a,dominated=b,before=dominates(old[a],old[b]),after=dominates(new[a],new[b]),
                    before_margins={f:old[a][f]-old[b][f] for f in FIELDS},after_margins={f:new[a][f]-new[b][f] for f in FIELDS}))
    orows={r['key']:r for r in before['residues']}; nrows={r['key']:r for r in after['residues']}
    return dict(changed_patches=changed,changed_dominance_edges=edges,
        positive_scores_restored=[dict(key=k,before=orows[k]['local_score'],after=nrows[k]['local_score']) for k in orows if orows[k]['local_score']==0<nrows[k]['local_score']],
        max_residue_abs_delta={f:max((abs(orows[k][f]-nrows[k][f]) for k in orows),default=0.) for f in ('scrsa','local_score','propensity','multiscale_persistence')},
        common_center_membership_unchanged=all(old[k]['members']==new[k]['members'] and old[k]['seed_members']==new[k]['seed_members'] for k in common),
        top5_unchanged=[(p['center_key'],p['members']) for p in before['patches'][:5]]==[(p['center_key'],p['members']) for p in after['patches'][:5]])


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ('before','after','recomputed','output'):ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args()
    legacy,source_hash=load_legacy()
    manifest=json.loads(Path(__file__).with_name('inputs.json').read_text())
    cases=[c for c in manifest['cases'] if not c['id'].startswith('raw_')]
    assert len(cases)==19
    missing=[c['id'] for c in cases if not (args.recomputed/(c['id']+'.json')).is_file()]
    if missing:
        raise SystemExit('Finish the full reproduction run first; missing: '+', '.join(missing))
    report=dict(baseline_commit=BASE,legacy_core_sha256=source_hash,only_kernel_changes='20 round calls removed; module docstring excluded',cases=[],changed_maps=[])
    for case in cases:
        paths=[args.before/(case['id']+'.json'),args.after/(case['id']+'.json')]
        before,after=[json.loads(p.read_text()) for p in paths]
        assert before['input_sha256']==after['input_sha256']==case['p0_capture_sha256']
        # Historical Windows captures used CRLF; committed snapshots use LF.
        canonical=checked_input(case).encode('utf-8')
        assert hashlib.sha256(canonical.replace(b'\n',b'\r\n')).hexdigest()==case['p0_capture_sha256']
        with patch.object(engine,'core',legacy):
            replay=engine.analyze_all_maps(pdb_text=checked_input(case),chain=case['chain'],pH=case['pH'],ionic_mM=case['ionic_mM'],temp_K=case['temp_K'])
        fresh=json.loads((args.recomputed/(case['id']+'.json')).read_text())
        assert replay['maps']==before['result']['maps'],'Historical replay mismatch: '+case['id']
        assert fresh['maps']==after['result']['maps'],'P0 recomputation mismatch: '+case['id']
        for key,old in before['result']['maps'].items():
            new=after['result']['maps'][key]
            if [p['center_key'] for p in old['patches']]!=[p['center_key'] for p in new['patches']]:
                detail=dict(case=case['id'],map=key,**details(old,new))
                if case['id']=='extreme_1CRN_0.0':
                    # Tiny fixture: expose the candidate priorities behind suppression.
                    candidates={}
                    for label,kernel in [('before',legacy),('after',core)]:
                        rows=kernel._surface_residues(fresh['prepared_pdb'],case['pH'])
                        surface=[r for r in rows if r['surface_exposed']]
                        channel=kernel.chemistry_map(surface,kernel.build_distances(surface),key,case['pH'])
                        candidates[label]=[{f:c[f] for f in ('center_key','multiscale_persistence','multiscale_geomean')}
                            for c in channel['patch_centers'] if c['center_key'] in {'A:19:','A:20:','A:22:','A:24:'}]
                    detail['selection_candidates']=candidates
                report['changed_maps'].append(detail)
        report['cases'].append(dict(id=case['id'],input_sha256=case['sha256'],p0_capture_sha256=case['p0_capture_sha256'],maps=11,legacy_replay_exact=True,p0_recomputation_exact=True,
            archive_sha256={label:hashlib.sha256(p.read_bytes()).hexdigest() for label,p in zip(('before','after'),paths)}))
        print(case['id'],'209-map replay PASS',flush=True)
    assert len(report['changed_maps'])==3
    report['status']='PASS'; report['map_count']=209
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')


if __name__=='__main__':main()
