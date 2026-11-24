"""Real competitor-tool preflight. Never turn an installation failure into a result.

Scientific scope:
- PEP-Patch: charged and hydrophobic protein surface patches, closest comparator.
- SURFMAP: 2D projection of hydrophobic surface feature; 2D pixels ≠ 3D patch.
- P2Rank: small-molecule pocket prediction, different task; operational only.

Each candidate uses the same frozen ISV76 normalized coordinates. We do not
calculate an invalid unified accuracy/speed figure across unlike methods.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import psutil
import shutil
import subprocess
import sys
import tempfile
import time
import traceback

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"validation"),str(ROOT/"backend")]
from isv_full_80_qa import prepare
from prepare import prepare_pdb_text

ACTIVE={r["proposed_pdb_id"]:r for r in csv.DictReader(
    (ROOT/"validation/ISV_active_candidates_76_v3.csv").open(encoding="utf-8")
)}
PIDS=("1UBQ","2OED","1HEL","1GCT")
OUT=Path("external_methods_preflight").resolve()
TOOLS=("pep_patch","surfmap","p2rank")


def save_json(p,data):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")


def run(cmd,cwd,timeout):
    cmd=list(map(str,cmd))
    start=time.monotonic()
    with tempfile.TemporaryFile(mode="w+t") as stdout, tempfile.TemporaryFile(mode="w+t") as stderr:
        env=dict(os.environ)
        if any("pep_patch_" in arg for arg in cmd):
            env["PATH"]="/opt/pep-venv/bin:"+env.get("PATH","")
        p=psutil.Popen(cmd,cwd=str(cwd),stdout=stdout,stderr=stderr,env=env)
        peak=0
        while p.poll() is None:
            try:
                children=p.children(recursive=True)
                peak=max(peak,p.memory_info().rss+sum(
                    c.memory_info().rss for c in children if c.is_running()))
            except (psutil.Error,OSError):pass
            if time.monotonic()-start>timeout:
                # APBS/PDB2PQR may spawn children. Killing only the parent
                # strands numerical processes and contaminates later cases.
                try:
                    descendants=p.children(recursive=True)
                except psutil.Error:
                    descendants=[]
                for child in reversed(descendants):
                    try: child.kill()
                    except psutil.Error: pass
                try: p.kill()
                except psutil.Error: pass
                psutil.wait_procs([p,*descendants],timeout=12)
                stdout.seek(0);stderr.seek(0)
                return {"status":"TIMEOUT","returncode":-9,"wall_s":round(time.monotonic()-start,3),
                        "peak_RSS_sampled_MiB":round(peak/1048576,3),
                        "stdout_tail":stdout.read()[-3000:],
                        "stderr_tail":stderr.read()[-5000:]}
            time.sleep(.075)
        p.wait()
        stdout.seek(0);stderr.seek(0)
        out,err=stdout.read(),stderr.read()
    return {"status":"PASS" if p.returncode==0 else "FAIL","returncode":p.returncode,
            "wall_s":round(time.monotonic()-start,3),
            "peak_RSS_sampled_MiB":round(peak/1048576,3),
            "stdout_tail":out[-3500:],"stderr_tail":err[-6000:]}


def input_pdb(pid,path):
    meta,src=prepare(ACTIVE[pid])
    text,_=prepare_pdb_text(src.read_text(encoding="utf-8"))
    path.write_text(text,encoding="utf-8")
    return {"prepared_PDB_sha256":hashlib.sha256(text.encode()).hexdigest(),
            "original_assembly_sha256":meta["assembly_cif_sha256"],
            "PDB_id":pid}


def inventory(path):
    files=[x for x in path.rglob("*") if x.is_file()]
    return {"n_files":len(files),
            "total_bytes":sum(p.stat().st_size for p in files),
            "by_extension":{e:sum(x.suffix.lower()==e for x in files)
                            for e in sorted({x.suffix.lower() for x in files})},
            "paths":[str(x.relative_to(path)) for x in files[:50]]}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--tool",choices=TOOLS,required=True)
    args=ap.parse_args()
    tool=args.tool
    runroot=OUT/tool;runroot.mkdir(parents=True,exist_ok=True)
    meta={"tool":tool,"PDBs":list(PIDS),
          "scope":{"pep_patch":"3D electrostatic/hydrophobic patches",
                   "surfmap":"2D Kyte-Doolittle hydrophobicity projection",
                   "p2rank":"small-molecule pockets, different scientific endpoint"}[tool],
          "scientific_warning":"No common binding-site accuracy or performance superiority claim",
          "measurements":[]}
    for pid in PIDS:
        d=runroot/pid;d.mkdir(parents=True,exist_ok=True)
        try:
            inp=d/(pid+".pdb")
            provenance=input_pdb(pid,inp)
            cmds=[]
            if tool=="pep_patch":
                hyd=d/"hydrophobic";hyd.mkdir(exist_ok=True)
                # PEP-Patch bundled Crippen/Eisenberg AMBER lookup does not
                # include generic HIS backbone N for these inputs. Run its
                # documented user-provided hydropathy-scale functionality,
                # with one fixed, conventional Kyte-Doolittle value per aa.
                kd={
                    "ALA":1.8,"ARG":-4.5,"ASN":-3.5,"ASP":-3.5,
                    "CYS":2.5,"GLN":-3.5,"GLU":-3.5,"GLY":-0.4,
                    "HIS":-3.2,"ILE":4.5,"LEU":3.8,"LYS":-3.9,
                    "MET":1.9,"PHE":2.8,"PRO":-1.6,"SER":-0.8,
                    "THR":-0.7,"TRP":-0.9,"TYR":-1.3,"VAL":4.2,
                }
                scale=hyd/"kyte_doolittle_scale.csv"
                scale.write_text("residue,propensity\n"+"".join(
                    k+","+str(v)+"\n" for k,v in kd.items()),encoding="utf-8")
                ele=d/"electrostatic";ele.mkdir(exist_ok=True)
                program="/opt/pep-venv/bin"
                cmds=[
                    ("hydrophobic", [program+"/pep_patch_hydrophobic",
                         inp,inp,"--scale",scale,"--potential",
                         "--patches","--out",hyd/"hydrophobic.npz",
                         "--ply_out",hyd/"hydrophobic.ply"],hyd,500),
                    ("electrostatic", [program+"/pep_patch_electrostatic",
                         inp,inp,"--apbs_dir",ele/"apbs",
                         "-o",ele/"patches.csv","-ro",ele/"residues.csv",
                         "--ply_out",ele/"electrostatic"],ele,900),
                ]
            elif tool=="surfmap":
                d2=d/"surfmap_output";d2.mkdir(exist_ok=True)
                cmds=[
                    ("kyte_doolittle", ["/opt/surf-venv/bin/surfmap","-pdb",inp,
                        "-tomap","kyte_doolittle","-d",d2,"--docker"],d,1000)
                ]
            else:
                pr="/opt/p2rank_2.5.1/prank"
                output=d/"pockets";output.mkdir(exist_ok=True)
                cmds=[
                    ("pocket_prediction",["bash",pr,"predict","-f",inp,"-o",output,
                                          "-threads","1"],d,800)
                ]
            attempts=[]
            for label,cmd,cwd,limit in cmds:
                t=run(cmd,cwd,limit)
                attempts.append({"task":label,"command":list(map(str,cmd)),
                                 "runtime":t})
                print("EXTERNAL_PREFLIGHT "+json.dumps({
                    "tool":tool,"pdb":pid,"task":label,"status":t["status"],
                    "seconds":t["wall_s"],"memory":t["peak_RSS_sampled_MiB"],
                    "error":t["stderr_tail"][-350:]}),flush=True)
            item={**provenance,"attempts":attempts,"output_inventory":inventory(d),
                  "overall_status":"PASS" if all(a["runtime"]["status"]=="PASS" for a in attempts)
                                   else "PARTIAL_OR_FAIL"}
        except Exception as exc:
            item={"PDB_id":pid,"overall_status":"FAIL","exception":str(exc),
                  "traceback":traceback.format_exc()[-3000:]}
        meta["measurements"].append(item)
        save_json(runroot/(pid+"_result.json"),item)
    # Pin actual external installation provenance, no fabricated successful
    # execution when a program or Docker daemon is absent.
    installlog=Path("/opt/external_setup_status.txt")
    if installlog.exists():
        meta["dependency_installation_log"]=installlog.read_text(encoding="utf-8")[-14000:]
    meta["n_pass"]=sum(x["overall_status"]=="PASS" for x in meta["measurements"])
    meta["n_total"]=len(PIDS)
    save_json(runroot/"SUMMARY.json",meta)
    print("EXTERNAL_PREFLIGHT_SUMMARY "+json.dumps({
        "tool":tool,"pass":meta["n_pass"],"total":meta["n_total"]}),flush=True)
    # Failure is a meaningful observation, but the workflow must retain
    # detailed output and not conflate failed cases with verified metrics.
    if meta["n_pass"]!=meta["n_total"]:
        raise SystemExit("Some external comparator cases not executable; see structured evidence")


if __name__=="__main__":
    main()
