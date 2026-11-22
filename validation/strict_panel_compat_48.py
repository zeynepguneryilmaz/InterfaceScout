"""Actual all-map InterfaceScout compatibility and assembly preparation.
Runs every PDB in the frozen 48-ID v3 list. Does NOT certify experimental validity.
"""
import csv, json, hashlib, sys, traceback, subprocess, time, string
from pathlib import Path
from collections import Counter
import gemmi
ROOT=Path(__file__).resolve().parents[1]
OUT=Path("strict_panel_48")
OUT.mkdir(exist_ok=True)
IDS=[r["pdb_id"].upper() for r in csv.DictReader(open(ROOT/"validation"/"strict_48_selected_candidates_v1.csv",encoding="utf-8"))]
assert len(IDS)==48 and len(set(IDS))==48

def prepare(pid):
    file=OUT/f"{pid}-assembly1.cif"
    if not file.exists():
        import urllib.request
        url=f"https://files.rcsb.org/download/{pid}-assembly1.cif"
        with urllib.request.urlopen(url,timeout=60) as u:file.write_bytes(u.read())
    st=gemmi.read_structure(str(file))
    if not len(st):raise RuntimeError("No structural model")
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.remove_ligands_and_waters()
    model=st[0]
    alph=string.ascii_uppercase+string.ascii_lowercase+string.digits
    mapping=[]
    kept=[]
    for c in model:
        cnt=sum(1 for res in c.get_polymer() if res.name in {
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET",
            "PHE","PRO","SER","THR","TRP","TYR","VAL"})
        if cnt:kept.append(c)
    if len(kept)>len(alph):raise RuntimeError("Too many chains to write single-character PDB IDs")
    # Construct an entirely NEW protein-only model. In-place removal from a
    # Gemmi Model can leave unrelated DNA or generated chains in some assemblies.
    clean_structure=gemmi.Structure()
    clean_structure.name=pid
    clean_structure.cell=st.cell
    clean_structure.spacegroup_hm=st.spacegroup_hm
    clean_model=gemmi.Model(1)
    for i,c in enumerate(kept):
        orig=c.name
        copied=c.clone()
        copied.name=alph[i]
        clean_model.add_chain(copied)
        mapping.append({"source_chain":orig,"prepared_chain":copied.name,"standard_residues":sum(
            1 for r in copied.get_polymer() if r.name in {
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET",
            "PHE","PRO","SER","THR","TRP","TYR","VAL"})})
    clean_structure.add_model(clean_model)
    clean_structure.remove_ligands_and_waters()
    pdb=clean_structure.make_pdb_string()
    if not pdb.count("\nATOM"):raise RuntimeError("Empty ATOM PDB conversion")
    out=OUT/f"{pid}_prepared.pdb"
    out.write_text(pdb,encoding="utf-8")
    return {"pdb":pid,"chain_mapping":mapping,
            "n_chains":len(mapping),"n_expected_standard_residues":sum(x["standard_residues"] for x in mapping),
            "pdb_sha256":hashlib.sha256(pdb.encode("utf-8")).hexdigest(),"pdb_bytes":len(pdb.encode("utf-8"))}

def worker(pid):
    sys.path.insert(0,str(ROOT/"backend"))
    from interface_engine import analyze_all_maps, PUBLIC_MAP_SPECS
    file=OUT/f"{pid}_prepared.pdb"
    pdb=file.read_text(encoding="utf-8")
    result=analyze_all_maps(pdb_text=pdb,chain=None,pH=7.4,ionic_mM=150.,temp_K=298.)
    n=result.get("n_maps")
    if n!=11:raise ValueError(f"expected 11, got {n}")
    order=result.get("map_order")
    if len(order)!=11 or len(set(order))!=11:raise ValueError("map order invalid")
    profiles=result.get("residue_profiles",[])
    if len(profiles)<1:raise ValueError("No residue profiles")
    keys=[str(p["key"]) for p in profiles]
    if len(set(keys))!=len(keys):raise ValueError("Duplicate residue key (chain collision)")
    maps=result["maps"]
    if set(maps)!=set(order):raise ValueError("Map keys differ from map_order")
    checks={}
    for key in order:
        data=maps[key]
        patches=data.get("patches",[])
        residues=data.get("residues",[])
        if len(residues)>len(profiles):raise ValueError(f"invalid map size for {key}")
        if any(int(p.get("pareto_front",0))<1 for p in patches):raise ValueError("Invalid Pareto front")
        checks[key]={"n_patches":len(patches),"n_rows":len(residues),"n_front1":sum(p.get("pareto_front")==1 for p in patches)}
    print("WORKER_OK "+json.dumps({"pdb":pid,"n_maps":n,"n_profiles":len(profiles),
        "n_surface":result["diagnostics"]["n_surface_residues"],
        "maps":checks},separators=(",",":")),flush=True)

def main():
    mode=sys.argv[1] if len(sys.argv)>1 else "all"
    if mode=="--worker":
        try:worker(sys.argv[2])
        except BaseException:
            traceback.print_exc();sys.exit(2)
        return
    rows=[]
    for idx,pid in enumerate(IDS,1):
        item={"pdb":pid}
        try:
            item.update(prepare(pid))
            res=subprocess.run([sys.executable,__file__,"--worker",pid],capture_output=True,text=True,timeout=420)
            item["returncode"]=res.returncode
            item["stdout_tail"]=res.stdout[-3000:]
            item["stderr_tail"]=res.stderr[-3000:]
            if res.returncode==0 and "WORKER_OK " in res.stdout:
                payload=json.loads(res.stdout.split("WORKER_OK ",1)[1].strip().splitlines()[0])
                item["engine"]="PASS"
                item["n_maps"]=payload["n_maps"]
                item["n_residue_profiles"]=payload["n_profiles"]
                item["n_surface_residues"]=payload["n_surface"]
                item["map_details"]=payload["maps"]
            else:item["engine"]="FAIL"
        except subprocess.TimeoutExpired:
            item["engine"]="TIMEOUT_420S"
        except Exception:
            item["engine"]="PREP_FAIL"
            item["error"]=traceback.format_exc()[-3000:]
        rows.append(item)
        print("COMPAT "+json.dumps({"index":idx,"pdb":pid,"engine":item["engine"],
            "residues":item.get("n_residue_profiles"),"n_maps":item.get("n_maps"),
            "error":item.get("error",item.get("stderr_tail",""))[:350]},separators=(",",":")),flush=True)
    (OUT/"compatibility_48.json").write_text(json.dumps(rows,indent=2))
    with (OUT/"compatibility_48.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["pdb","engine","n_chains","n_expected_standard_residues",
            "n_residue_profiles","n_surface_residues","n_maps","pdb_sha256","error"])
        w.writeheader()
        for r in rows:w.writerow({k:r.get(k,"") for k in w.fieldnames})
    print("COMPAT_SUMMARY "+json.dumps({"n":len(rows),
      "engine_pass":sum(x["engine"]=="PASS" for x in rows),
      "engine_fail":[{"pdb":x["pdb"],"status":x["engine"],"error":x.get("error",x.get("stderr_tail",""))[:500]} for x in rows if x["engine"]!="PASS"]},separators=(",",":")),flush=True)
if __name__=="__main__":main()
