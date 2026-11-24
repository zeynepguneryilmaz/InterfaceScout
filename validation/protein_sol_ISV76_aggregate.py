"""Strict ISV76 Protein-sol full-panel audit and result completeness.

Outputs 76 row-level results (PASS, partial, failure). Actual server PDB
contents are required: HTTP 200 and a results page alone never count as PASS.
Provider B factors are not asserted to share physical units with InterfaceScout
APBS values. Two source-exported PDBs must both be nonempty and parseable.
"""
from __future__ import annotations
import csv,json,statistics,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
DATA=Path("protein_sol_ISV76")
OUT=DATA/"final";OUT.mkdir(parents=True,exist_ok=True)
manifest=list(csv.DictReader((ROOT/"validation"/"ISV_active_candidates_76_v3.csv").open(encoding="utf-8",newline="")))
ids=[q["proposed_pdb_id"] for q in manifest]
assert len(ids)==len(set(ids))==76
rows=[]
def keyset(path):
    ids=set();n=0
    for line in path.read_text(encoding="utf-8",errors="replace").splitlines():
        if not line.startswith(("ATOM  ","HETATM")) or len(line)<27:continue
        try:
            key=(line[21].strip(),int(line[22:26].strip()),line[26].strip())
            ids.add(key);n+=1
        except ValueError:pass
    return ids,n
for q in manifest:
    pid=q["proposed_pdb_id"]
    home=DATA/"proteins"/pid
    summary=home/"SUMMARY.json"
    item={"pdb_id":pid,"category":q["category"],"status":"MISSING",
          "job_id":None,"elapsed_s":None,"input_sha256":None,
          "potential_PDB_valid":False,"ratio_PDB_valid":False,
          "n_input_residues":None,"n_potential_residues":None,
          "n_ratio_residues":None,"potential_input_residue_coverage":None,
          "ratio_input_residue_coverage":None,"error":""}
    if summary.exists():
        obj=json.loads(summary.read_text(encoding="utf-8"))
        item.update({"status":obj["status"],
                     "job_id":obj.get("job_id"),
                     "elapsed_s":obj.get("elapsed_s"),
                     "input_sha256":obj.get("input_sha256"),
                     "error":obj.get("error","")})
        if (home/"input_prepared.pdb").exists():
            orig,natom=keyset(home/"input_prepared.pdb")
            item["n_input_residues"]=len(orig)
            for mode,label in [("potlfield","potential"),("ratio","ratio")]:
                path=home/(mode+".pdb")
                if not path.exists():continue
                residue,atoms=keyset(path)
                item[label+"_PDB_valid"]=atoms>=4
                item["n_"+label+"_residues"]=len(residue)
                item[label+"_input_residue_coverage"]=round(
                    len(orig & residue)/max(1,len(orig)),7)
        if item["status"]=="SUCCESS_BOTH_PDB" and not (
            item["potential_PDB_valid"] and item["ratio_PDB_valid"]):
            item["status"]="INVALID_BOTH_PDB_CLAIM"
    rows.append(item)
outpath=OUT/"PROTEIN_SOL_76_ALL_CASES.csv"
with outpath.open("w",encoding="utf-8",newline="") as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]))
    writer.writeheader();writer.writerows(rows)
complete=[q for q in rows if q["status"]=="SUCCESS_BOTH_PDB"]
times=[float(q["elapsed_s"]) for q in rows if q["elapsed_s"] is not None]
stats={"expected":76,"observed":len(rows),"n_with_manifest_and_result":sum(q["status"]!="MISSING" for q in rows),
       "n_with_both_valid_scientific_PDBs":len(complete),
       "n_incomplete":76-len(complete),
       "n_with_job_id":sum(bool(q["job_id"]) for q in rows),
       "status_counts":{k:sum(q["status"]==k for q in rows) for k in sorted(set(q["status"] for q in rows))},
       "median_client_end_to_end_elapsed_s":statistics.median(times) if times else None,
       "source_pdb_provenance":"Frozen ISV76, normalized PDB",
       "server_memory":"Cannot be measured externally",
       "comparison_scope":"Descriptive provider PDB coverage and surface region attributes, not adsorption-site predictive accuracy",
       "success_definition":"Both genuine downloaded potlfield.pdb and ratio.pdb are valid atomic PDBs",
       "cannot_compare_raw_values":"Protein-sol FDPB/NPP model and InterfaceScout APBS/chemical descriptors differ; do not assert identical units",
       "all_76_results_accounted":all(q["status"]!="MISSING" for q in rows)}
(OUT/"SUMMARY.json").write_text(json.dumps(stats,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
print("PROTEIN_SOL_ALL76_FINAL "+json.dumps(stats,ensure_ascii=False),flush=True)
if not stats["all_76_results_accounted"]:
    sys.exit("Incomplete 76-member result provenance. Missing items are retained.")
if stats["n_incomplete"]:
    sys.exit("Some Protein-sol jobs did not yield both actual output PDBs. See CSV.")
