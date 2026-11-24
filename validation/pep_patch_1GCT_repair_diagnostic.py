"""Diagnose missing-heavy-atom preparation for the frozen 1GCT comparator input.

The PDBFixer-prepared structure is a SECONDARY diagnostic representation.
It must never replace the frozen original PDB in the ISV76 main benchmark.
No scientific accuracy inference is made from the resulting repaired model.
"""
from __future__ import annotations
import hashlib, json, sys, traceback
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"validation"),str(ROOT/"backend")]
from isv_full_80_qa import prepare
from prepare import prepare_pdb_text
import csv
OUT=Path("pep_patch_1GCT_repair_diagnostic");OUT.mkdir(exist_ok=True)
PDB="1GCT"
record={"pdb_id":PDB,"status":"FAIL","interpretation":"Secondary repaired-structure control, not original-input benchmark"}

def write():
    (OUT/"SUMMARY.json").write_text(json.dumps(record,indent=2,ensure_ascii=False,default=str)+"\n",encoding="utf-8")

def atom_keys(text):
    result={}
    for line in text.splitlines():
        if line.startswith(("ATOM  ","HETATM")) and len(line)>=27:
            key=(line[21].strip(),line[22:26].strip(),line[26].strip())
            result.setdefault(key,set()).add(line[12:16].strip())
    return result

def main():
    try:
        manifest=list(csv.DictReader((ROOT/"validation/ISV_active_candidates_76_v3.csv").open(encoding="utf-8")))
        row=next(q for q in manifest if q["proposed_pdb_id"]==PDB)
        meta,pdb=prepare(row)
        original,_=prepare_pdb_text(pdb.read_text(encoding="utf-8"))
        (OUT/"original_input.pdb").write_text(original,encoding="utf-8")
        record["original_SHA256"]=hashlib.sha256(original.encode()).hexdigest()
        record["original_assembly_sha256"]=meta.get("assembly_cif_sha256")
        atoms=atom_keys(original)
        record["source_residue_count"]=len(atoms)
        incomplete={":".join(k):sorted({"N","CA","C","O"}-v) for k,v in atoms.items()
                    if {"N","CA","C","O"}-v}
        record["original_incomplete_backbone_residues"]=incomplete
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
        fixer=PDBFixer(filename=str(OUT/"original_input.pdb"))
        fixer.findMissingResidues()
        record["omitted_complete_residue_insertions"]={str(k):str(v) for k,v in fixer.missingResidues.items()}
        fixer.missingResidues={} # Never silently rebuild whole unresolved sequences
        fixer.findNonstandardResidues()
        record["nonstandard_residues"]=[str(v) for v in fixer.nonstandardResidues]
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        record["missing_atoms_detected"]={str(k):[str(a) for a in v] for k,v in fixer.missingAtoms.items()}
        record["missing_terminal_atoms_detected"]={str(k):[str(a) for a in v] for k,v in fixer.missingTerminals.items()}
        fixer.addMissingAtoms()
        path=OUT/"repaired_missing_atoms_only.pdb"
        with path.open("w",encoding="utf-8") as handle:
            PDBFile.writeFile(fixer.topology,fixer.positions,handle,keepIds=True)
        repaired=path.read_text(encoding="utf-8")
        record["repaired_SHA256"]=hashlib.sha256(repaired.encode()).hexdigest()
        record["repaired_residue_count"]=len(atom_keys(repaired))
        record["repaired_remaining_incomplete_backbone"]={
            ":".join(k):sorted({"N","CA","C","O"}-v)
            for k,v in atom_keys(repaired).items() if {"N","CA","C","O"}-v}
        record["status"]="REPAIRED_REPRESENTATION_CREATED"
    except Exception as exc:
        record.update(status="DIAGNOSTIC_FAILED",error=str(exc),traceback=traceback.format_exc()[-5000:])
    finally:
        write()
        print("PEP_PATCH_1GCT_REPAIR_DIAGNOSTIC "+json.dumps(record,default=str),flush=True)
    if record["status"]!="REPAIRED_REPRESENTATION_CREATED":
        raise SystemExit("1GCT repair diagnostic did not produce validated structure")

if __name__=="__main__": main()
