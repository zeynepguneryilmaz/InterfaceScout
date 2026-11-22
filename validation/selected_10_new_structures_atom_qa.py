"""Actual RCSB mmCIF atom completeness QA for new and comparison structures.
Reported flags apply to MODELED residues in first model only.
Coordinate sequence coverage, mutation identity, and wwPDB geometry are separate checks.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import re
import urllib.request
from pathlib import Path
from Bio.PDB import MMCIFParser, MMCIF2Dict

PDBS = "1AKI 1AO6 1FNF 1IGT 2LZM 4EUL 1F6S 3V83 1DEQ 1JY3".split()
AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
ATOMS = {
"ALA":"N CA C O CB","ARG":"N CA C O CB CG CD NE CZ NH1 NH2",
"ASN":"N CA C O CB CG OD1 ND2","ASP":"N CA C O CB CG OD1 OD2",
"CYS":"N CA C O CB SG","GLN":"N CA C O CB CG CD OE1 NE2",
"GLU":"N CA C O CB CG CD OE1 OE2","GLY":"N CA C O",
"HIS":"N CA C O CB CG ND1 CD2 CE1 NE2",
"ILE":"N CA C O CB CG1 CG2 CD1","LEU":"N CA C O CB CG CD1 CD2",
"LYS":"N CA C O CB CG CD CE NZ","MET":"N CA C O CB CG SD CE",
"PHE":"N CA C O CB CG CD1 CD2 CE1 CE2 CZ","PRO":"N CA C O CB CG CD",
"SER":"N CA C O CB OG","THR":"N CA C O CB OG1 CG2",
"TRP":"N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
"TYR":"N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH",
"VAL":"N CA C O CB CG1 CG2",
}
ATOMS = {k:set(v.split()) for k,v in ATOMS.items()}
OUT = Path("new_structures_atom_qa")
OUT.mkdir(exist_ok=True)

def audit(pid):
    path = OUT / f"{pid}.cif"
    if not path.exists():
        req=urllib.request.Request(f"https://files.rcsb.org/download/{pid}.cif",
                                   headers={"User-Agent":"InterfaceScout-scientific-validation/1.0"})
        with urllib.request.urlopen(req,timeout=60) as src:
            path.write_bytes(src.read())
    structure = MMCIFParser(QUIET=True,auth_chains=True,auth_residues=True).get_structure(pid,str(path))
    raw = MMCIF2Dict.MMCIF2Dict(str(path))
    rows=[]
    models=list(structure.get_models())
    for ch in models[0]:
        modeled=[r for r in ch if r.id[0]==" " and r.resname in AA]
        if not modeled:continue
        detail=[]
        for r in modeled:
            names={a.name for a in r}
            missing=ATOMS[r.resname]-names
            bb=sorted(missing & {"N","CA","C","O"})
            side=sorted(missing-{"N","CA","C","O"})
            alt=any(a.is_disordered() for a in r)
            if bb or side or alt:
                detail.append({"chain":ch.id,"auth_residue_number":r.id[1],"icode":r.id[2].strip(),
                               "resname":r.resname,"missing_backbone_atoms":bb,
                               "missing_sidechain_atoms":side,"alternate_conformation":bool(alt)})
        rows.append({"chain":ch.id,"modeled_standard_residues":len(modeled),
                     "missing_backbone_residues":sum(bool(x["missing_backbone_atoms"]) for x in detail),
                     "missing_sidechain_residues":sum(bool(x["missing_sidechain_atoms"]) for x in detail),
                     "altloc_residues":sum(bool(x["alternate_conformation"]) for x in detail),
                     "flagged_positions":detail})
    entities=raw.get("_entity_poly.pdbx_seq_one_letter_code_can",[])
    expected=[len(re.sub("[^A-Z]","",s.replace("\n","").upper())) for s in entities]
    return {"pdb":pid,"status":"COORDINATES_AUDITED","source_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
            "number_of_models":len(models),"chain_count":len(rows),"standard_modeled_residues":sum(x["modeled_standard_residues"] for x in rows),
            "missing_backbone_residues":sum(x["missing_backbone_residues"] for x in rows),
            "missing_sidechain_residues":sum(x["missing_sidechain_residues"] for x in rows),
            "altloc_residues":sum(x["altloc_residues"] for x in rows),
            "deposited_polymer_sequence_lengths":expected,"chains":rows,
            "method_note":"First model; modeled canonical amino acids only; NOT whole-chain sequence/wwPDB local fit validation"}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--shard",type=int,default=0)
    parser.add_argument("--shards",type=int,default=2)
    args=parser.parse_args()
    assert args.shards>0 and 0<=args.shard<args.shards
    results=[]
    for pid in PDBS[args.shard::args.shards]:
        try: rec=audit(pid)
        except Exception as exc: rec={"pdb":pid,"status":"ERROR","error":str(exc)}
        results.append(rec)
        print("ATOM_QA",pid,rec["status"],rec.get("standard_modeled_residues"),rec.get("error",""),flush=True)
    target=OUT / f"shard_{args.shard}"
    target.mkdir(exist_ok=True)
    (target/"atom_quality.json").write_text(json.dumps(results,indent=2))
    keys=["pdb","status","source_sha256","number_of_models","chain_count","standard_modeled_residues",
          "missing_backbone_residues","missing_sidechain_residues","altloc_residues",
          "deposited_polymer_sequence_lengths","chains","method_note","error"]
    with (target/"atom_quality.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,keys)
        writer.writeheader()
        for rec in results:
            writer.writerow({k:json.dumps(rec.get(k)) if isinstance(rec.get(k),(list,dict)) else rec.get(k,"")
                             for k in keys})
    if any(x.get("status")!="COORDINATES_AUDITED" for x in results):
        raise RuntimeError("Some selected structures could not be QA audited")

if __name__=="__main__":
    main()
