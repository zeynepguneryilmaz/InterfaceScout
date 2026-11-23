"""Surface-localized missing-atom and biological assembly audit, ISV80.

All scores are from the unchanged 2.5.0 production engine. Defects are
associated with actual *surface and top-five candidate* residue keys, not
global missing-atom totals. This does not impute missing atoms or claim their
unmeasured energetic effect. Cofactor removal is logged, not assumed harmless.
"""
from __future__ import annotations
import csv,json,hashlib,sys,urllib.request,time,collections,traceback
from pathlib import Path
import gemmi
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"validation"))
sys.path.insert(0,str(ROOT/"backend"))
from isv_full_80_qa import prepare,CANON,SC,ROOTOUT
from interface_engine import analyze_all_maps
from model_settings import MODEL_VERSION
from concurrent.futures import ThreadPoolExecutor
PANEL=list(csv.DictReader((ROOT/"validation/ISV_protein_candidate_pool_80_v2.csv").open()))
FOCUS={"3GHG","1QJS","1STP","4DN0","6LQH"}
OUT=Path("isv80_scientific_audit");OUT.mkdir(exist_ok=True)
IGNORE={"HOH","DOD"}
COMMON_SALTS={"NA","CL","SO4","PO4","GOL","EDO","PEG","ACT","ACY","MES","TRS","HEP","BME","MPD","FMT","IOD","NO3"}
METALS={"ZN","CA","CU","FE","MN","NI","MG","CO","CD","K","NA"}
COFACTORS={"HEM","HEC","HEA","FAD","FMN","NAD","NAP","NDP","ADP","ATP","PLP","PQQ","SAM","SAH"}
def getjson(url):
    for k in range(3):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"InterfaceScout-scientific-structure-QA/1.1","Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=25) as r:return json.load(r)
        except Exception:
            if k==2:raise
            time.sleep(k+1)
def source_info(pid):
    e=getjson(f"https://data.rcsb.org/rest/v1/core/entry/{pid}")
    assembly_ids=(e.get("rcsb_entry_container_identifiers") or {}).get("assembly_ids") or []
    assemblies=[]
    for aid in assembly_ids[:10]:
        a=getjson(f"https://data.rcsb.org/rest/v1/core/assembly/{pid}/{aid}")
        info=a.get("rcsb_assembly_info") or {}
        assemblies.append({"id":str(aid),
          "assembly_name":(a.get("pdbx_struct_assembly") or {}).get("details",""),
          "oligomeric_count":(a.get("pdbx_struct_assembly") or {}).get("oligomeric_count",""),
          "polymer_entity_instances":info.get("polymer_entity_instance_count"),
          "polymer_composition":info.get("polymer_composition"),
          "selected_assembly_1":str(aid)=="1"})
    return {"pdb":pid,"assembly_ids":assembly_ids,"assemblies":assemblies,
        "rcsb_assembly1_reported":any(x["id"]=="1" for x in assemblies),
        "pdb_description":(e.get("struct") or {}).get("title",""),
        "experimental_method":";".join(x.get("method","") for x in e.get("exptl") or [])}
def structure_notes(cif_path,chain_mapping):
    st=gemmi.read_structure(str(cif_path))
    st.remove_alternative_conformations()
    notes=[]
    ligands=collections.Counter()
    per_chain={}
    for i,chain in enumerate(st[0]):
        aa=[r for r in chain.get_polymer() if r.name in CANON]
        if not aa:continue
        prepared=chain_mapping[i]["prepared_chain"]
        per_chain[prepared]={"source_chain":chain.name,"protein_residues":len(aa),
          "sequence_min":min(r.seqid.num for r in aa),
          "sequence_max":max(r.seqid.num for r in aa)}
        for r in chain:
            if r.name in CANON or r.name in IGNORE:continue
            ligands[r.name]+=1
    for name,count in ligands.most_common():
        notes.append({"component":name,"count":count,
            "class":"metal_or_ion" if name in METALS else "cofactor" if name in COFACTORS else "small_molecule_or_noncanonical",
            "screening_relevance":"recorded and omitted from protein-only model; interaction effect not quantified"})
    return notes,per_chain
def check_one(row):
    pid=row["proposed_pdb_id"]
    prep,pdbfile=prepare(row)
    model=analyze_all_maps(pdb_text=pdbfile.read_text(),pH=7.4,chain=None)
    assert model["version"]==MODEL_VERSION and model["n_maps"]==11
    profiles=model["residue_profiles"]
    key_map={p["key"]:p for p in profiles}
    exposed={p["key"] for p in profiles if p.get("surface_exposed")}
    top5_by_map={}
    top5_members=set()
    for ch in model["map_order"]:
        patches=model["maps"][ch]["patches"][:5]
        members=set()
        for patch in patches:members.update(patch.get("members",[]))
        top5_by_map[ch]=members
        top5_members.update(members)
    # Same chain-remapping order as prepare(): standard-polymer-bearing
    # chains only, but do not silently drop all source chains with no standard AA.
    cifpath=ROOTOUT/(pid+"-assembly1.cif")
    st=gemmi.read_structure(str(cifpath))
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    issue=[]
    atom_summary=collections.Counter()
    for orig,chain in zip(prep["chains"],
        [c for c in st[0] if any(r.name in CANON for r in c.get_polymer())]):
        newchain=orig["prepared_chain"]
        for res in chain.get_polymer():
            if res.name not in CANON:continue
            atoms={a.name for a in res if a.element.name not in ("H","D")}
            missing_side=sorted(set(SC[res.name].split())-atoms)
            missing_backbone=sorted({"N","CA","C","O"}-atoms)
            if not missing_side and not missing_backbone:continue
            icode=str(res.seqid.icode).strip()
            key=f"{newchain}:{res.seqid.num}:{icode}"
            # profile numbering may use converted author numbering;
            # verify mapping rather than treating all unmatched as buried.
            profile=key_map.get(key)
            if profile is None:
                atom_summary["missing_residue_from_engine_profile"]+=1
                issue.append({"pdb":pid,"prepared_key":key,"source_chain":chain.name,
                    "res_name":res.name,"missing_sidechain_heavy_atoms":";".join(missing_side),
                    "missing_backbone_heavy_atoms":";".join(missing_backbone),
                    "engine_surface_exposed":"UNRESOLVED_PROFILE_ABSENT",
                    "engine_trsa":"","engine_scrsa":"","n_top5_maps":0,"top5_maps":""})
                continue
            is_surface=key in exposed
            maps=[k for k,v in top5_by_map.items() if key in v]
            if missing_side:atom_summary["incomplete_sidechain_residues"]+=1
            if missing_backbone:atom_summary["incomplete_backbone_residues"]+=1
            if is_surface and missing_side:atom_summary["exposed_incomplete_sidechain_residues"]+=1
            if maps and missing_side:atom_summary["top5_incomplete_sidechain_residues"]+=1
            if maps and missing_backbone:atom_summary["top5_incomplete_backbone_residues"]+=1
            issue.append({"pdb":pid,"prepared_key":key,"source_chain":chain.name,
                "res_name":res.name,"missing_sidechain_heavy_atoms":";".join(missing_side),
                "missing_backbone_heavy_atoms":";".join(missing_backbone),
                "engine_surface_exposed":int(is_surface),
                "engine_trsa":profile.get("trsa_raw"),
                "engine_scrsa":profile.get("scrsa"),
                "n_top5_maps":len(maps),"top5_maps":";".join(maps)})
    notes,chainmeta=structure_notes(cifpath,prep["chains"])
    src=source_info(pid)
    if not src["rcsb_assembly1_reported"]:raise ValueError("Source assembly 1 not listed in RCSB")
    if src["assemblies"][0]["id"]!="1":raise ValueError("Selected assembly unexpectedly non-first")
    n_modified=sum(x["engine_surface_exposed"] for x in issue)
    res={"pdb":pid,"protein":row["protein_name"],"category":row["category"],
       "model_version":MODEL_VERSION,"assembly_cif_sha256":prep["assembly_cif_sha256"],
       "protein_only_sha256":prep["prepared_pdb_sha256"],
       "n_residues":len(profiles),"n_exposed":len(exposed),"n_chains":prep["n_chains"],
       "missing_sidechain_residues":atom_summary["incomplete_sidechain_residues"],
       "missing_residue_from_engine_profile":atom_summary["missing_residue_from_engine_profile"],
       "exposed_missing_sidechain_residues":atom_summary["exposed_incomplete_sidechain_residues"],
       "top5_missing_sidechain_residues":atom_summary["top5_incomplete_sidechain_residues"],
       "missing_backbone_residues":atom_summary["incomplete_backbone_residues"],
       "top5_missing_backbone_residues":atom_summary["top5_incomplete_backbone_residues"],
       "n_top5_unique_members":len(top5_members),
       "missing_sidechain_fraction_of_exposed":atom_summary["exposed_incomplete_sidechain_residues"]/len(exposed) if exposed else 0,
       "source_assembly_count":len(src["assembly_ids"]),
       "source_assembly_ids":";".join(src["assembly_ids"]),
       "assembly_1_protein_chains":prep["n_chains"],
       "source_bioassembly_details":src["assemblies"],
       "component_names":";".join(x["component"] for x in notes),
       "n_distinct_removed_components":len(notes),
       "removed_components":notes,"chain_info":chainmeta,
       "warnings":[]}
    if res["missing_residue_from_engine_profile"]:res["warnings"].append("UNMAPPED_MISSING_ATOM_RESIDUE_REQUIRES_REVIEW")
    if res["exposed_missing_sidechain_residues"]:res["warnings"].append("MISSING_SURFACE_CHEMISTRY_UNCERTAINTY")
    if res["top5_missing_sidechain_residues"]:res["warnings"].append("TOP5_REGION_ATOM_INCOMPLETE")
    if res["n_distinct_removed_components"]:res["warnings"].append("OMITTED_NONSOLVENT_COMPONENTS")
    if res["source_assembly_count"]>1:res["warnings"].append("MULTIPLE_SOURCE_BIOASSEMBLIES_REVIEW")
    if pid in FOCUS:
        (OUT/f"focus_{pid}_missing_atoms.json").write_text(json.dumps(issue,indent=2))
        (OUT/f"focus_{pid}_assembly_components.json").write_text(json.dumps({"source":src,"removed":notes,"chains":chainmeta},indent=2))
    return res,issue
def run(shard,shards):
    records=[];issues=[];errors=[]
    for i,row in enumerate(PANEL):
        if i%shards!=shard:continue
        try:
            result,detail=check_one(row)
            records.append(result);issues.extend(detail)
            print("ISV_SCIENCE_OK",row["proposed_pdb_id"],
                result["n_exposed"],result["exposed_missing_sidechain_residues"],
                result["top5_missing_sidechain_residues"],result["component_names"],flush=True)
        except Exception as ex:
            errors.append({"pdb":row["proposed_pdb_id"],"error":str(ex),"trace":traceback.format_exc()[-1700:]})
            print("ISV_SCIENCE_FAIL",row["proposed_pdb_id"],str(ex),flush=True)
    (OUT/f"science_{shard}.json").write_text(json.dumps(records,indent=2))
    (OUT/f"errors_{shard}.json").write_text(json.dumps(errors,indent=2))
    cols=["pdb","prepared_key","source_chain","res_name","missing_sidechain_heavy_atoms","missing_backbone_heavy_atoms","engine_surface_exposed","engine_trsa","engine_scrsa","n_top5_maps","top5_maps"]
    with (OUT/f"issues_{shard}.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(issues)
    if errors:raise SystemExit(2)
def aggregate():
    records=[];errors=[]
    for shard in range(8):
        records+=json.loads((OUT/f"science_{shard}.json").read_text())
        errors+=json.loads((OUT/f"errors_{shard}.json").read_text())
    if errors or len(records)!=80:raise ValueError(f"Incomplete audit {len(records)} records {len(errors)} errors")
    byid={r["pdb"]:r for r in records}
    if len(byid)!=80 or not FOCUS.issubset(byid):raise ValueError("Missing records")
    src={r["proposed_pdb_id"]:r for r in PANEL}
    counts=collections.Counter()
    for item in records:
        if item["exposed_missing_sidechain_residues"]>0:counts["proteins_with_exposed_missing_sidechain"]+=1
        if item["top5_missing_sidechain_residues"]>0:counts["proteins_with_top5_missing_sidechain"]+=1
        if item["n_distinct_removed_components"]>0:counts["proteins_with_omitted_nonsolvent_components"]+=1
        if item["source_assembly_count"]>1:counts["proteins_with_alternative_assembly_record"]+=1
        counts["total_exposed_missing_sidechain_residues"]+=item["exposed_missing_sidechain_residues"]
        counts["total_top5_missing_sidechain_residues"]+=item["top5_missing_sidechain_residues"]
    allcols=["pdb","protein","category","n_chains","n_residues","n_exposed","missing_sidechain_residues",
      "exposed_missing_sidechain_residues","top5_missing_sidechain_residues","missing_backbone_residues",
      "top5_missing_backbone_residues","missing_residue_from_engine_profile","missing_sidechain_fraction_of_exposed",
      "n_distinct_removed_components","component_names","source_assembly_count","warnings"]
    with (OUT/"ISV80_SURFACE_QA.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=allcols);w.writeheader()
        for r in sorted(records,key=lambda v:v["pdb"]):
            w.writerow({k:";".join(r[k]) if isinstance(r.get(k),list) else r.get(k,"") for k in allcols})
    (OUT/"ISV80_SURFACE_QA_DETAILS.json").write_text(json.dumps(records,indent=2))
    result={"status":"PASS","n_proteins":80,"model":MODEL_VERSION,**counts,
      "focused_five":{p:{k:v for k,v in byid[p].items() if k not in ("source_bioassembly_details","chain_info","removed_components")} for p in sorted(FOCUS)},
      "interpretation":"Audits missing side-chain heavy atoms and whether their residue keys fall in geometric surface/top5 mapped regions. Actual corrected physicochemical maps or cofactors' energetic impacts were NOT calculated."}
    (OUT/"SUMMARY.json").write_text(json.dumps(result,indent=2))
    print("ISV_SURFACE_AUDIT_SUMMARY",json.dumps(result),flush=True)
if __name__=="__main__":
    if len(sys.argv)==3:run(int(sys.argv[1]),int(sys.argv[2]))
    elif len(sys.argv)==2 and sys.argv[1]=="aggregate":aggregate()
    else:raise SystemExit("Usage: script.py shard 8 | aggregate")
