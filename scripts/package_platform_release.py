"""Package InterfaceScout's tracked application files for a specific OS.

Release archives intentionally contain source launchers and NOT a native,
self-contained signed executable. Python 3.10-3.12 is required at installation.
No developer-only validation datasets, CI internals or runtime caches included.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import zipfile

ROOT=Path(__file__).resolve().parents[1]
VERSION="1.1.0"
COMMON_ROOT={
    "README.md","LICENSE","CITATION.cff","THIRD_PARTY_NOTICES.md",
}
KEEP_TOP={"backend","frontend","assets","launchers","docs","examples"}
OS_LAUNCHERS={
    "windows":{"run_local.bat","start.bat","Create Desktop Shortcut.bat"},
    "macos":{"run_local.sh","start.command","Create Desktop App.command"},
    "linux":{"run_local.sh","start.sh","Create Desktop Shortcut.sh"},
}
EXCLUDE_SUFFIXES={".pyc",".pyo",".log",".tmp",".bak"}
def allow(path:str,platform:str)->bool:
    p=PurePosixPath(path)
    if path in COMMON_ROOT:return True
    if path in set().union(*OS_LAUNCHERS.values()):
        return path in OS_LAUNCHERS[platform]
    if not p.parts or p.parts[0] not in KEEP_TOP:return False
    if len(p.parts)>1 and p.parts[1] in {".venv","__pycache__"}:return False
    if any(x in {".venv","__pycache__","node_modules",".DS_Store"} for x in p.parts):return False
    if p.suffix.lower() in EXCLUDE_SUFFIXES:return False
    if p.parts[0]=="examples" and "results" in p.parts and p.suffix.lower() in {".zip",".json",".csv",".xlsx",".pdb"}:return False
    return True

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--platform",choices=OS_LAUNCHERS,required=True)
    ap.add_argument("--output-dir",default="dist/release")
    a=ap.parse_args()
    target=Path(a.output_dir);target.mkdir(parents=True,exist_ok=True)
    files=subprocess.check_output(["git","ls-files","-z"],cwd=ROOT).decode("utf-8").split("\0")
    included=sorted({f for f in files if f and allow(f,a.platform)})
    required=COMMON_ROOT | OS_LAUNCHERS[a.platform] | {
        "backend/app.py","backend/requirements.txt","frontend/index.html",
        "frontend/logo.png","backend/dependency_check.py",
        f"docs/INSTALL_{a.platform.upper() if a.platform != 'macos' else 'MACOS'}.md",
    }
    if a.platform=="windows":
        required.update({"assets/InterfaceScout.ico",
                         "launchers/create_windows_shortcut.ps1"})
    missing=sorted(required-set(included))
    if missing:raise SystemExit("Missing required package files: "+", ".join(missing))
    prefix="InterfaceScout-"+VERSION
    name=f"{prefix}-{a.platform}.zip"
    archive=target/name
    report=f"""# InterfaceScout {VERSION} — {a.platform} release ZIP

This is the full local application with OS-specific setup/launch scripts.
It is NOT a compiled native installer or a signed standalone executable.
Python 3.10–3.12 and a compatible internet connection for initial pip
installation are required unless dependencies are installed separately.
The local interface runs at http://localhost:8000.

First launch:
{ {'windows':'Double-click run_local.bat','macos':'In Terminal: bash run_local.sh','linux':'In Terminal: bash run_local.sh'}[a.platform] }

Later launches:
{ {'windows':'Double-click start.bat','macos':'In Terminal: bash start.command','linux':'In Terminal: bash start.sh'}[a.platform] }

Optional desktop shortcut:
{ {'windows':'Create Desktop Shortcut.bat','macos':'bash "Create Desktop App.command"','linux':'bash "Create Desktop Shortcut.sh"'}[a.platform] }

Unzip fully before starting; do not run launchers from within a ZIP.
Optional APBS/PDB2PQR/PROPKA electrostatics requires separate installation
and is not enabled by default. Original eleven chemical channels work without it.
Please see docs/INSTALL_{a.platform.upper() if a.platform != 'macos' else 'MACOS'}.md.
"""
    with zipfile.ZipFile(archive,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=8) as z:
        for rel in included:
            path=ROOT/rel
            if not path.is_file():raise SystemExit("Tracked file missing from checkout: "+rel)
            info=zipfile.ZipInfo(prefix+"/"+rel)
            info.compress_type=zipfile.ZIP_DEFLATED
            info.external_attr=((0o100755 if path.suffix in {".sh",".command"} else 0o100644)<<16)
            z.writestr(info,path.read_bytes(),compress_type=zipfile.ZIP_DEFLATED,compresslevel=8)
        z.writestr(prefix+"/START_HERE.txt",report)
    with zipfile.ZipFile(archive) as z:
        invalid=z.testzip()
        if invalid:raise SystemExit("Archive CRC failed: "+invalid)
        members=set(z.namelist())
        for rel in required:
            assert prefix+"/"+rel in members, "Required member absent: "+rel
        for plat,launchers in OS_LAUNCHERS.items():
            if plat!=a.platform:
                for rel in launchers-OS_LAUNCHERS[a.platform]:
                    assert prefix+"/"+rel not in members, "Foreign OS launcher leaked: "+rel
    digest=hashlib.sha256(archive.read_bytes()).hexdigest()
    (target/(name+".sha256")).write_text(digest+"  "+name+"\n",encoding="utf-8")
    manifest={"version":VERSION,"platform":a.platform,"archive":name,
      "sha256":digest,"file_count":len(included),"source_commit":
      subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
      "distributable":"Python-based local application; not a compiled or signed executable",
      "python_required":"3.10-3.12","offline_wheelhouse_included":False}
    (target/(name+".json")).write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(manifest,ensure_ascii=False))

if __name__=="__main__":main()
