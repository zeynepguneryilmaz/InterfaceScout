"""Independent residue/side-chain bookkeeping QA for InterfaceScout.
Executes the real core._surface_residues, compares atom-level SASA summed
by a separately written reference routine. Does not change production code.
"""
import sys, json, math
from io import StringIO
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"backend"))
import core
from Bio.PDB import Structure, Model, Chain, Residue, Atom, PDBIO, PDBParser
from Bio.PDB.SASA import ShrakeRupley

def construct():
    st=Structure.Structure("QA"); model=Model.Model(0);st.add(model)
    # include unusual numbering/insertion code and duplicate residue index across chains
    entries=[("A",10," ","ALA",["N","CA","C","O","CB"]),
             ("A",11,"A","SER",["N","CA","C","O","CB","OG"]),
             ("A",12," ","GLY",["N","CA","C","O"]),
             ("B",10," ","PHE",["N","CA","C","O","CB","CG","CD1","CD2","CE1","CE2","CZ"]),
             ("B",13," ","GLY",["N","CA","C","O"])]
    chains={}
    for chainid,seq,icode,rname,names in entries:
        if chainid not in chains:
            ch=Chain.Chain(chainid);model.add(ch);chains[chainid]=ch
        rr=Residue.Residue((" ",seq,icode),rname," ");chains[chainid].add(rr)
        for i,name in enumerate(names):
            # crowd some neighbors but keep all coordinates unique
            center={"A":0,"B":17}[chainid]+(seq-10)*4
            coords=[center+i*.46, (i%3)*.7, (i%2)*.6]
            element=name[0]
            atom=Atom.Atom(name,coords,10.,1.," ",name,i+1,element=element)
            rr.add(atom)
    io=PDBIO();io.set_structure(st);buffer=StringIO();io.save(buffer)
    return buffer.getvalue()

pdb=construct()
rows=core._surface_residues(pdb,7.4)
ref=PDBParser(QUIET=True).get_structure("ref",StringIO(pdb))
ShrakeRupley(probe_radius=core.SASA_PROBE_A,n_points=core.SASA_POINTS).compute(ref,level="A")
expected={}
for chain in next(ref.get_models()):
    for residue in chain:
        rn=residue.resname
        key=f"{chain.id}:{residue.id[1]}:{residue.id[2].strip()}"
        if rn=="GLY":
            selected=[residue["CA"]]
        else:
            selected=[atom for atom in residue if atom.get_name().strip() not in ("N","CA","C","O","OXT") and atom.element!="H"]
        area=sum(float(a.sasa) for a in selected)
        raw=area/core.SIDECHAIN_REF_ASA[rn]
        expected[key]=(rn,area,raw,raw>=core.SC_RSA_THRESHOLD)
found={r["key"]:r for r in rows}
assert len(found)==len(rows)==len(expected),(len(found),len(expected))
records=[]
for key,(rn,area,raw,exposed) in expected.items():
    r=found[key]
    assert r["res_name"]==rn,(key,"type")
    assert abs(r["sidechain_sasa"]-area)<=.000501,(key,"area",r["sidechain_sasa"],area)
    assert abs(r["scrsa_raw"]-raw)<=.00000501,(key,"relative",r["scrsa_raw"],raw)
    assert r["surface_exposed"]==exposed,(key,"threshold")
    assert r["scrsa"]==round(min(max(raw,0.),1.),5),(key,"clamping")
    records.append({"key":key,"residue":rn,"sidechain_A2":round(area,3),"scrsa_raw":round(raw,5),"surface_exposed":exposed})
assert "A:11:A" in found and "B:10:" in found and "A:10:" in found
# Explicitly document that GLY uses CA pseudo-sidechain, not chemically defined sidechain.
gly=[x for x in records if x["residue"]=="GLY"]
assert len(gly)==2
print("RESIDUE_QA "+json.dumps({"status":"PASS","n_residues":len(rows),"n_chains":2,
    "insertion_code":"PASS","duplicate_number_different_chain":"PASS",
    "sidechain_atom_selection":"PASS","raw_and_clamped_scRSA":"PASS",
    "surface_threshold":"PASS","glycine_CA_surrogate":gly,
    "scope":"Synthetic multi-chain structure; independent aggregation of same Biopython SASA kernel"}))
