"""Context audit: flag protein-nucleic acid complexes, bound ligands and metals.

A model can execute successfully yet misrepresent biological context when
protein-only preparation removes a bound cofactor or interaction partner.
"""
import csv,json,requests,time
from pathlib import Path
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import gemmi

ROOT=Path(__file__).resolve().parents[1]
IDS=[x["pdb_id"] for x in csv.DictReader(open(ROOT/"validation"/"strict_48_selected_candidates_v1.csv",encoding="utf-8"))]
P=Path("strict_context_48");P.mkdir(exist_ok=True)
STD={"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
NUCLEIC={"A","U","C","G","I","DA","DT","DC","DG","DI","DU"}
METALS={"ZN","MN","MG","CA","FE","FE2","CU","CO","CD","NI","MO","W","K","NA","SR","BA","HG","PB"}
EXCLUDE_WATER={"HOH","DOD","WAT"}
def fetch(pid):
    e=None
    for tries in range(3):
        try:
            q=requests.get(f"https://files.rcsb.org/download/{pid}-assembly1.cif",timeout=50)
            q.raise_for_status();return q.content
        except Exception as err:e=err;time.sleep(tries+1)
    raise e
def run(pid):
    try:
        p=P/(pid+"-assembly1.cif")
        p.write_bytes(fetch(pid))
        st=gemmi.read_structure(str(p))
        aa=0;nucl=0;metals=[];ligands=[];small=[]
        for chain in st[0]:
            for res in chain:
                name=res.name.upper().strip()
                if name in STD:aa+=1;continue
                if name in EXCLUDE_WATER:continue
                atoms=[a for a in res if a.element.name not in ("H","D")]
                if name in NUCLEIC or (name.startswith("D") and name[1:] in NUCLEIC):
                    nucl+=1;continue
                item={"name":name,"chain":chain.name,"number":str(res.seqid),"heavy_atoms":len(atoms)}
                if name in METALS or (len(atoms)==1 and atoms[0].element.name.upper() in METALS):metals.append(item)
                elif len(atoms)>=6:ligands.append(item)
                else:small.append(item)
        flags=[]
        if nucl:flags.append("BOUND_DNA_RNA")
        if metals:flags.append("BOUND_METAL_ION")
        if ligands:flags.append("NONTRIVIAL_LIGAND")
        return {"pdb":pid,"protein_residues":aa,"nucleic_residues":nucl,
          "metal_count":len(metals),"metals":metals[:20],
          "large_ligand_count":len(ligands),"large_ligands":ligands[:25],
          "small_hetero_count":len(small),"small_heteros":small[:10],
          "flags":flags}
    except Exception as e:return {"pdb":pid,"error":str(e),"flags":["FETCH_OR_PARSE_FAILED"]}
with ThreadPoolExecutor(max_workers=12) as ex:rows=list(ex.map(run,IDS))
for x in rows:print("CONTEXT_RECORD "+json.dumps({k:v for k,v in x.items() if k not in ("small_heteros","metals","large_ligands")},separators=(",",":")),flush=True)
(P/"context_48.json").write_text(json.dumps(rows,indent=2))
with (P/"context_48.csv").open("w",newline="") as file:
 w=csv.DictWriter(file,fieldnames=["pdb","nucleic_residues","metal_count","large_ligand_count","small_hetero_count","flags","error"])
 w.writeheader()
 for r in rows:w.writerow({k:(";".join(v) if isinstance(v,list) else v) for k,v in r.items() if k in w.fieldnames})
summary={"n":len(rows),"n_clean":sum(not q["flags"] for q in rows),
 "n_bound_nucleic":sum("BOUND_DNA_RNA" in q["flags"] for q in rows),
 "n_bound_metals":sum("BOUND_METAL_ION" in q["flags"] for q in rows),
 "n_large_ligands":sum("NONTRIVIAL_LIGAND" in q["flags"] for q in rows),
 "flagged":[{"pdb":q["pdb"],"flags":q["flags"]} for q in rows if q["flags"]]}
(P/"context_summary.json").write_text(json.dumps(summary,indent=2))
print("CONTEXT_SUMMARY "+json.dumps(summary,separators=(",",":")),flush=True)
