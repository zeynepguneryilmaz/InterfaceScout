"""ISV structure-level and production-engine QA for all 80 proposed PDB structures.

This is a candidate screening workflow, NOT experimental-site validation.
It assesses assembly processing, standard-residue heavy-atom sidechains, and the
unchanged production 11-map engine. All limitations are output, never repaired
by silently adding/mutating coordinates.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,subprocess,sys,time,traceback
from pathlib import Path
from collections import Counter
import gemmi
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"validation/ISV_protein_candidate_pool_80_v2.csv"
ROOTOUT=Path("isv_80_full_qa");ROOTOUT.mkdir(exist_ok=True)
CANON={"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
SC={"ALA":"CB","ARG":"CB CG CD NE CZ NH1 NH2","ASN":"CB CG OD1 ND2","ASP":"CB CG OD1 OD2","CYS":"CB SG","GLN":"CB CG CD OE1 NE2","GLU":"CB CG CD OE1 OE2","GLY":"","HIS":"CB CG ND1 CD2 CE1 NE2","ILE":"CB CG1 CG2 CD1","LEU":"CB CG CD CE NZ","LYS":"CB CG CD CE NZ","MET":"CB CG SD CE","PHE":"CB CG CD1 CD2 CE1 CE2 CZ","PRO":"CB CG CD","SER":"CB OG","THR":"CB OG1 CG2","TRP":"CB CG CD1 NE1 CE2 CE3 CZ2 CZ3 CH2","TYR":"CB CG CD1 CD2 CE1 CE2 CZ OH","VAL":"CB CG1 CG2"}
ASCII="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
def prepare(row):
    import urllib.request
    pid=row["proposed_pdb_id"]
    cif=ROOTOUT/(pid+"-assembly1.cif")
    if not cif.exists():
        url=f"https://files.rcsb.org/download/{pid}-assembly1.cif"
        with urllib.request.urlopen(url,timeout=80) as f: cif.write_bytes(f.read())
    raw=cif.read_bytes()
    st=gemmi.read_structure(str(cif))
    if len(st)==0:raise ValueError("no structural models")
    st.remove_alternative_conformations();st.remove_hydrogens()
    model=st[0]
    chain_info=[]
    for chain in model:
        aa=[r for r in chain.get_polymer() if r.name in CANON]
        if not aa:continue
        sc_missing=0; backbone_missing=0; missing_atoms=0
        for residue in aa:
            atoms={atom.name for atom in residue if atom.element.name not in ("H","D")}
            backbone_missing+=not {"N","CA","C","O"}.issubset(atoms)
            absent=set(SC[residue.name].split())-atoms
            sc_missing+=bool(absent);missing_atoms+=len(absent)
        chain_info.append({"source_chain":chain.name,"residues":len(aa),
           "backbone_incomplete_residues":backbone_missing,
           "sidechain_incomplete_residues":sc_missing,
           "missing_sidechain_atom_count":missing_atoms})
    if not chain_info:raise ValueError("no protein polymer chains")
    if len(chain_info)>len(ASCII):raise ValueError("more than 62 protein chains; PDB chain-ID limit")
    clean=gemmi.Structure();clean.name=pid;clean.cell=st.cell;clean.spacegroup_hm=st.spacegroup_hm
    clean_model=gemmi.Model(1)
    for i,old in enumerate([c for c in model if any(r.name in CANON for r in c.get_polymer())]):
        cp=old.clone();cp.name=ASCII[i];clean_model.add_chain(cp)
        chain_info[i]["prepared_chain"]=ASCII[i]
    clean.add_model(clean_model);clean.remove_ligands_and_waters()
    out=ROOTOUT/(pid+"_protein_only.pdb")
    pdb=clean.make_pdb_string()
    if "ATOM  " not in pdb:raise ValueError("PDB conversion produced zero ATOM lines")
    out.write_text(pdb,encoding="utf8")
    summary={"pdb":pid,"candidate_id":row["candidate_id"],"protein":row["protein_name"],
      "category":row["category"],"assembly_cif_sha256":hashlib.sha256(raw).hexdigest(),
      "prepared_pdb_sha256":hashlib.sha256(pdb.encode()).hexdigest(),"n_chains":len(chain_info),
      "n_residues":sum(x["residues"] for x in chain_info),
      "sidechain_incomplete_residues":sum(x["sidechain_incomplete_residues"] for x in chain_info),
      "backbone_incomplete_residues":sum(x["backbone_incomplete_residues"] for x in chain_info),
      "n_missing_sidechain_heavy_atoms":sum(x["missing_sidechain_atom_count"] for x in chain_info),
      "chains":chain_info,"source_cif_path":str(cif),"protein_only_path":str(out)}
    return summary,out
def worker(filename,pid):
    sys.path.insert(0,str(ROOT/"backend"))
    from interface_engine import analyze_all_maps
    from model_settings import MODEL_VERSION
    r=analyze_all_maps(pdb_text=Path(filename).read_text(),pH=7.4,chain=None)
    assert r["version"]==MODEL_VERSION
    assert r["n_maps"]==11 and len(r["maps"])==11 and len(r["map_order"])==11
    profiles=r["residue_profiles"]
    keys=[p["key"] for p in profiles]
    assert len(keys)==len(set(keys)) and len(keys)>0
    exposed=[p for p in profiles if p["surface_exposed"]]
    assert all(p["surface_exposed"]==(p["trsa_raw"]>=0.05) for p in profiles)
    assert all(p["scrsa"]==0 and p["sidechain_sasa"]==0 for p in profiles if p["res_name"]=="GLY")
    meta={}
    for mapname in r["map_order"]:
        d=r["maps"][mapname];patches=d.get("patches",[])
        assert all(int(p["pareto_front"])>=1 for p in patches)
        meta[mapname]={"patches":len(patches),"local_residues":len(d["residues"])}
    print("ISV_ENGINE_RESULT "+json.dumps({"pdb":pid,"model":r["version"],
      "n_residues":len(profiles),"n_exposed":len(exposed),"n_maps":11,
      "n_nonempty_patch_maps":sum(v["patches"]>0 for v in meta.values()),
      "map_counts":meta}),flush=True)
def do_shard(shard,shards,timeout=480):
    rows=list(csv.DictReader(DATA.open()))
    assert len(rows)==80 and len({x["proposed_pdb_id"] for x in rows})==80
    records=[]
    for i,row in enumerate(rows):
        if i%shards!=shard:continue
        pid=row["proposed_pdb_id"]
        try:
            summary,filename=prepare(row)
            flags=[]
            if summary["backbone_incomplete_residues"]:flags.append("BACKBONE_INCOMPLETE")
            if summary["sidechain_incomplete_residues"]:flags.append("SIDECHAIN_INCOMPLETE")
            if summary["n_chains"]>12:flags.append("LARGE_MULTICHAIN_ASSEMBLY")
            if summary["n_residues"]>2500:flags.append("LARGE_PROTEIN_SYSTEM")
            if summary["n_residues"]<50:flags.append("SMALL_PROTEIN")
            proc=subprocess.run([sys.executable,__file__,"--worker",str(filename),pid],
                capture_output=True,text=True,timeout=timeout)
            payload=None
            if proc.returncode==0 and "ISV_ENGINE_RESULT " in proc.stdout:
                payload=json.loads(proc.stdout.split("ISV_ENGINE_RESULT ",1)[1].splitlines()[0])
                engine="PASS"
            else:
                engine="FAIL";flags.append("ENGINE_FAILURE")
            summary.update({"engine_status":engine,"engine":payload,
                "engine_error":proc.stderr[-1400:],"quality_flags":flags,
                "screening_decision":"ENGINE_PASSED_REVIEW_REQUIRED" if engine=="PASS" else "NOT_QUALIFIED_YET"})
        except subprocess.TimeoutExpired:
            summary={"pdb":pid,"candidate_id":row["candidate_id"],"protein":row["protein_name"],"category":row["category"],
              "engine_status":"TIMEOUT","quality_flags":["ENGINE_TIMEOUT"],
              "screening_decision":"NOT_QUALIFIED_YET"}
        except Exception as e:
            summary={"pdb":pid,"candidate_id":row["candidate_id"],"protein":row["protein_name"],"category":row["category"],
              "engine_status":"PREPARATION_FAILURE","quality_flags":["PREPARATION_FAILURE"],
              "error":str(e),"trace":traceback.format_exc()[-1800:],
              "screening_decision":"NOT_QUALIFIED_YET"}
        records.append(summary)
        print("ISV80_QA "+json.dumps({"pdb":pid,"engine":summary["engine_status"],
          "n_residues":summary.get("n_residues"),"flags":summary["quality_flags"]}),flush=True)
    (ROOTOUT/f"details_{shard}.json").write_text(json.dumps(records,indent=2))
    with (ROOTOUT/f"screen_{shard}.csv").open("w",newline="",encoding="utf8") as f:
        fields=["pdb","candidate_id","protein","category","engine_status","screening_decision",
            "n_chains","n_residues","n_maps","n_exposed","n_nonempty_patch_maps",
            "sidechain_incomplete_residues","backbone_incomplete_residues",
            "assembly_cif_sha256","prepared_pdb_sha256","quality_flags","error"]
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for x in records:
            w.writerow({**{k:(json.dumps(x[k]) if isinstance(x.get(k),(list,dict)) else x.get(k,"")) for k in fields},
                "n_maps":(x.get("engine") or {}).get("n_maps",""),
                "n_exposed":(x.get("engine") or {}).get("n_exposed",""),
                "n_nonempty_patch_maps":(x.get("engine") or {}).get("n_nonempty_patch_maps","")})
    print("ISV80_SHARD_SUMMARY "+json.dumps({"shard":shard,"n":len(records),
      "engine_pass":sum(x["engine_status"]=="PASS" for x in records),
      "failed":sum(x["engine_status"]!="PASS" for x in records)}),flush=True)
def aggregate(path):
    detail=[]
    for i in range(8):
        d=path/f"details_{i}.json"
        if not d.exists():raise ValueError(f"Missing shard {i}")
        detail+=json.loads(d.read_text())
    if len(detail)!=80 or len({x["pdb"] for x in detail})!=80:raise ValueError("Incomplete 80-protein evidence")
    order={x["pdb"]:i for i,x in enumerate(csv.DictReader(DATA.open()))}
    detail.sort(key=lambda x:order[x["pdb"]])
    with (path/"ISV80_FULL_SCREEN.csv").open("w",newline="",encoding="utf8") as f:
        fields=["pdb","candidate_id","protein","category","engine_status","screening_decision",
          "n_chains","n_residues","n_maps","n_exposed","n_nonempty_patch_maps",
          "sidechain_incomplete_residues","backbone_incomplete_residues",
          "assembly_cif_sha256","prepared_pdb_sha256","quality_flags","error"]
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for x in detail:
            w.writerow({**{k:json.dumps(x.get(k)) if isinstance(x.get(k),(list,dict)) else x.get(k,"") for k in fields},
                "n_maps":(x.get("engine") or {}).get("n_maps",""),
                "n_exposed":(x.get("engine") or {}).get("n_exposed",""),
                "n_nonempty_patch_maps":(x.get("engine") or {}).get("n_nonempty_patch_maps","")})
    counts=dict(Counter(x["engine_status"] for x in detail))
    summary={"n_candidates":80,"engine_status_counts":counts,
        "n_engine_pass":counts.get("PASS",0),
        "n_with_zero_missing_backbone":sum(x.get("backbone_incomplete_residues",-1)==0 for x in detail),
        "n_without_missing_sidechains":sum(x.get("sidechain_incomplete_residues",-1)==0 for x in detail),
        "n_biofilm_engine_pass":sum(x["category"].startswith("biofilm_") and x["engine_status"]=="PASS" for x in detail),
        "statuses":"PASS indicates computational compatibility only; no final structural identity/sequence family acceptance or experimental adsorption validation",
        "next":"Manual domain/assembly/entity check; source evidence optional; wwPDB reports; family redundancy; then accept every suitable protein and freeze manifest."}
    (path/"SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print("ISV80_FULL_SUMMARY "+json.dumps(summary),flush=True)
if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--shard",type=int)
    parser.add_argument("--shards",type=int,default=8)
    parser.add_argument("--worker",nargs=2)
    parser.add_argument("--aggregate",type=Path)
    args=parser.parse_args()
    if args.worker:worker(*args.worker)
    elif args.aggregate:aggregate(args.aggregate)
    elif args.shard is not None:do_shard(args.shard,args.shards)
    else:parser.error("specify --shard or --aggregate")
