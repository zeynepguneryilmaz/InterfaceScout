"""Offline PDB-to-all-maps verification; no frozen score files are trusted."""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
HERE=Path(__file__).resolve().parent/'reproducibility'
sys.path.insert(0,str(ROOT/'backend'))
sys.path.insert(0,str(HERE))
import core
from interface_engine import analyze_all_maps
from oracle import verify_result, spec, ATOL, RTOL, BFACTOR_ATOL


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def checked_input(case):
    path=(ROOT/case['path']).resolve()
    if not path.is_relative_to(ROOT):
        raise ValueError('Input path outside repository')
    if sha(path)!=case['sha256']:
        raise ValueError('Input SHA-256 mismatch: '+case['id'])
    return path.read_text(encoding='utf-8')


def compute(case):
    captured=[]
    original=core.analyze
    def capture(request):
        result=original(request)
        captured.append(result)
        return result
    with patch.object(core,'analyze',side_effect=capture):
        result=analyze_all_maps(pdb_text=checked_input(case),chain=case['chain'],
            pH=case['pH'],ionic_mM=case['ionic_mM'],temp_K=case['temp_K'])
    if len(captured)!=1:
        raise AssertionError('Expected exactly one shared SASA/core pass')
    cfg=spec()
    if result['version']!=cfg['model_version']:
        raise AssertionError('Model version differs from declared specification')
    for k,value in cfg['parameters'].items():
        if result['method'][k]!=value:
            raise AssertionError('Parameter mismatch: '+k)
    return result,captured[0]


def verify_case(case,output):
    result,core_result=compute(case)
    checks=verify_result(result,core_result,case.get('independent_sasa',False))
    filename=output/(case['id']+'.json')
    filename.write_text(json.dumps(result,allow_nan=False),encoding='utf-8')
    run=subprocess.run(['node',str(HERE/'export_maps.cjs'),str(filename)],capture_output=True,text=True,encoding='utf-8')
    if run.returncode:
        raise AssertionError('B-factor verification failed: '+run.stderr)
    return dict(id=case['id'],input_sha256=case['sha256'],output_sha256=sha(filename),
                status='PASS',oracle=checks,bfactor=json.loads(run.stdout))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite',choices=['quick','full'],default='full')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    data=json.loads((HERE/'inputs.json').read_text())
    if sha(HERE/'model_spec.json')!=data['model_spec_sha256']:
        raise SystemExit('Declared model specification hash mismatch')
    cases=data['cases'] if args.suite=='full' else [c for c in data['cases'] if c['id'] in data['quick_case_ids']]
    def git(*arguments):
        run=subprocess.run(['git',*arguments],cwd=ROOT,capture_output=True,text=True)
        return run.stdout.strip() if run.returncode==0 else 'unavailable'
    report=dict(status='RUNNING',suite=args.suite,python=sys.version,platform=platform.platform(),
        dependencies={d.metadata['Name']:d.version for d in importlib.metadata.distributions()},
        node=subprocess.check_output(['node','--version'],text=True).strip(),
        git_commit=git('rev-parse','HEAD'),git_status=git('status','--porcelain'),
        model_version=spec()['model_version'],parameters=spec()['parameters'],
        tolerances=dict(atol=ATOL,rtol=RTOL,bfactor_atol=BFACTOR_ATOL),
        input_manifest_sha256=sha(HERE/'inputs.json'),spec_sha256=sha(HERE/'model_spec.json'),
        tolerance_contract_sha256=sha(HERE/'TOLERANCES.md'),
        source_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in
            list((ROOT/'backend').glob('*.py'))+[HERE/'oracle.py',HERE/'export_maps.cjs',Path(__file__),ROOT/'frontend/index.html']},
        cases=[])
    report_path=args.output/'verification.json'
    for case in cases:
        try:
            record=verify_case(case,args.output)
        except Exception as exc:
            record=dict(id=case['id'],status='FAIL',error=f'{type(exc).__name__}: {exc}')
        report['cases'].append(record)
        print(case['id'],record['status'],record.get('error',''),flush=True)
        report_path.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    report['status']='PASS' if all(c['status']=='PASS' for c in report['cases']) else 'FAIL'
    report_path.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(f"{report['status']}: {len(cases)} inputs, {11*len(cases)} maps; {report_path}")
    return 0 if report['status']=='PASS' else 1


if __name__=='__main__':
    raise SystemExit(main())
