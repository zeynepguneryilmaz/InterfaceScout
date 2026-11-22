"""End-to-end InterfaceScout smoke QC for five retained replacements plus 1RWZ."""
from __future__ import annotations
import os, sys, json, hashlib, tempfile
from pathlib import Path
import requests
import gemmi

IDS=["1C5E","1PDB","1A1X","1STN","2CBA","1RWZ"]
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"backend"))
from interface_engine import analyze_all_maps

OUT=Path("candidate_smoke_results");OUT.mkdir(exist_ok=True)
s=requests.Session()
s.headers["User-Agent"]="InterfaceScout-panel-smoke/1.0"
items=[]
for pid in IDS:
    print("START",pid,flush=True)
    record={"pdb":pid}
    try:
        url=f"https://files.rcsb.org/download/{pid}-assembly1.cif"
        response=s.get(url,timeout=45);response.raise_for_status()
        with tempfile.TemporaryDirectory() as folder:
            cif=Path(folder)/"assembly.cif";cif.write_bytes(response.content)
            st=gemmi.read_structure(str(cif))
            model=st[0]
            names=[]
            for i,ch in enumerate(model):
                if i>=26:raise ValueError("Chain count exceeds one-character PDB chain alphabet")
                names.append(ch.name)
                ch.name=chr(65+i)
            pdb=Path(folder)/"prepared_assembly.pdb"
            st.write_pdb(str(pdb))
            lines=[ln for ln in pdb.read_text().splitlines() if ln.startswith(("ATOM  ","TER   ","MODEL ","ENDMDL","END   "))]
            pdb_text="\n".join(lines)+"\n"
            record.update({"original_assembly_chain_names":names,"selected_pdb_chains":sorted({ln[21] for ln in lines if ln.startswith("ATOM  ")}),
                          "protein_pdb_sha256":hashlib.sha256(pdb_text.encode()).hexdigest(),"pdb_ATOM_count":sum(ln.startswith("ATOM  ") for ln in lines)})
            if record["pdb_ATOM_count"]<=0:raise ValueError("No ATOM coordinate records")
            result=analyze_all_maps(pdb_text=pdb_text,chain=None,pH=7.4,ionic_mM=150.0,temp_K=298.0)
            maps=result.get("maps",{})
            record.update({"success":result.get("n_maps")==11 and len(maps)==11,
                           "n_maps":result.get("n_maps"),"n_profiles":len(result.get("residue_profiles",[])),
                           "n_patches_by_map":{name:detail.get("n_patches",0) for name,detail in maps.items()},
                           "warnings":result.get("warnings",[])})
    except Exception as e:
        record.update({"success":False,"error":f"{type(e).__name__}: {e}"})
    items.append(record)
    print("SMOKE_QA",json.dumps(record,separators=(",",":")),flush=True)
(OUT/"results.json").write_text(json.dumps(items,indent=2))
if any(not x["success"] for x in items):
    raise SystemExit("One or more selected structures failed the full 11-map smoke check.")
