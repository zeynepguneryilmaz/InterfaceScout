"""RCSB screening of proposed ISV PDBs. No InterfaceScout performance labels used.

Outputs independently downloaded entry metadata, experimental method, resolution,
polymer entity identity/length, entity/chain mapping, and review flags.
A valid PDB ID alone is NOT a pass for protein identity or coordinate quality.
"""
from __future__ import annotations
import csv,json,time,urllib.request,urllib.error
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PANEL=ROOT/"validation/ISV_protein_candidate_pool_68_v1.csv"
TV=ROOT/"validation/strict_apo_48_selected_candidates_v2_final.csv"
OUT=Path("isv_rcsb_screen");OUT.mkdir(exist_ok=True)
def fetch(url):
    for i in range(3):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"InterfaceScout-ISV-Screen/1.0 (research metadata audit)","Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=30) as f:return json.load(f)
        except Exception as e:
            if i==2:raise
            time.sleep(2*(i+1))
def names_from_entity(record):
    p=record.get("entity_poly",{})
    ref=record.get("rcsb_polymer_entity",{})
    meta=record.get("rcsb_polymer_entity_container_identifiers",{})
    descr=record.get("pdbx_description","")
    if not descr:
        descr=(record.get("rcsb_polymer_entity",{}) or {}).get("pdbx_description","")
    return {"description":str(descr),
        "polymer_type":str(p.get("type","")),
        "sequence_len":len(str(p.get("pdbx_seq_one_letter_code_can","")).replace("\n","")),
        "entity_id":str(meta.get("entity_id","")),
        "asym_ids":";".join(meta.get("asym_ids") or []),
        "auth_asym_ids":";".join(meta.get("auth_asym_ids") or [])}
def main():
    rows=list(csv.DictReader(PANEL.open()))
    tv={r["pdb_id"] for r in csv.DictReader(TV.open())}
    assert len(rows)==68 and len({r["protein_name"] for r in rows})==68
    results=[]
    for r in rows:
        pdb=r["proposed_pdb_id"].upper()
        item={**r,"rcsb_url":f"https://www.rcsb.org/structure/{pdb}"}
        try:
            entry=fetch(f"https://data.rcsb.org/rest/v1/core/entry/{pdb}")
            ids=(entry.get("rcsb_entry_container_identifiers") or {}).get("polymer_entity_ids",[])
            entities=[]
            for eid in ids:
                try:
                    entity=fetch(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{eid}")
                    row=names_from_entity(entity)
                    if row["polymer_type"].lower().find("polypeptide")>=0:
                        entities.append(row)
                except Exception as e:
                    entities.append({"entity_id":eid,"entity_error":str(e)})
            method=";".join(str(x.get("method","")) for x in entry.get("exptl") or [])
            info=entry.get("rcsb_entry_info") or {}
            resolutions=info.get("resolution_combined") or []
            resolution=float(resolutions[0]) if resolutions else None
            complete=all("entity_error" not in x for x in entities)
            protein=bool(entities)
            item.update({"rcsb_metadata":"FOUND",
                 "experimental_method":method,"resolution_A":resolution,
                 "experimental_structure":int(bool(method) and (("X-RAY" in method.upper()) or ("ELECTRON" in method.upper()) or ("NMR" in method.upper()))),
                 "protein_entity_count":len(entities),
                 "protein_entities_json":json.dumps(entities,ensure_ascii=False),
                 "metadata_complete":int(complete),
                 "pdb_id_overlaps_TV48":int(pdb in tv)})
            # This is metadata screening only; identity/construct and coordinates
            # require subsequent independent per-chain manual inspection.
            flags=[]
            if pdb in tv:flags.append("EXACT_TV48_OVERLAP")
            if not protein:flags.append("NO_PROTEIN_POLYMER")
            if not complete:flags.append("ENTITY_METADATA_INCOMPLETE")
            if method.upper().find("X-RAY")!=-1 and resolution is not None and resolution>2.5:flags.append("XRAY_RESOLUTION_GT_2P5_REVIEW")
            if not method:flags.append("EXPERIMENT_METHOD_UNAVAILABLE")
            item["metadata_status"]="METADATA_REVIEW" if flags else "METADATA_PASS_ONLY"
            item["metadata_flags"]=";".join(flags)
            print("ISV_RCSB",pdb,item["metadata_status"],method,resolution,"entities",len(entities),flush=True)
        except Exception as e:
            item.update({"rcsb_metadata":"NOT_VERIFIED","experimental_method":"",
                 "resolution_A":"","experimental_structure":0,
                 "protein_entity_count":0,"protein_entities_json":"[]",
                 "metadata_complete":0,"pdb_id_overlaps_TV48":int(pdb in tv),
                 "metadata_status":"UNVERIFIED","metadata_flags":str(e)})
            print("ISV_RCSB_FAIL",pdb,str(e),flush=True)
        results.append(item)
    with (OUT/"isv_rcsb_metadata_68.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    totals={s:sum(x["metadata_status"]==s for x in results) for s in ("METADATA_PASS_ONLY","METADATA_REVIEW","UNVERIFIED")}
    report={"number_of_proteins":len(results),"metadata_counts":totals,
        "experimental_method_counts":{},
        "warning":"METADATA_PASS_ONLY is not final protein identity, experimental interface evidence, domain completeness, or coordinate quality acceptance. Manual sequence/construct and per-chain coordinate checks required."}
    for x in results:
        key=x["experimental_method"]
        report["experimental_method_counts"][key]=report["experimental_method_counts"].get(key,0)+1
    (OUT/"SUMMARY.json").write_text(json.dumps(report,indent=2,ensure_ascii=False))
    print("ISV_SCREEN_SUMMARY",json.dumps(report),flush=True)
    if totals["UNVERIFIED"]>0:raise SystemExit(2)
if __name__=="__main__":main()
