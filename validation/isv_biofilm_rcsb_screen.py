"""RCSB metadata audit of biofilm-related additions to ISV dataset 2.

Checks archival experimental methods, resolution, and polymer descriptions.
Does not certify direct binding contacts, missing atoms, biological assembly,
or methodological suitability of non-globular biofilm fibers.
"""
import csv,json,time,urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
IN=ROOT/"validation/ISV_protein_candidate_pool_80_v2.csv"
OUT=Path("isv_biofilm_audit");OUT.mkdir(exist_ok=True)
def get(url):
    for retry in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers={"User-Agent":"InterfaceScout-ISV-biofilm/1.0"}),timeout=25) as fd:
                return json.load(fd)
        except Exception:
            if retry==2:raise
            time.sleep(retry+1)
def main():
    source=list(csv.DictReader(IN.open()))
    entries=[r for r in source if r["category"].startswith("biofilm_")]
    assert len(source)==80 and len(entries)==12
    results=[]
    for r in entries:
        pdb=r["proposed_pdb_id"]
        entry=get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb}")
        polymer_ids=(entry.get("rcsb_entry_container_identifiers") or {}).get("polymer_entity_ids",[])
        entities=[]
        for entity_id in polymer_ids:
            entity=get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{entity_id}")
            pol=entity.get("entity_poly",{})
            if "polypeptide" not in str(pol.get("type","")).lower():continue
            description=(entity.get("rcsb_polymer_entity") or {}).get("pdbx_description") or entity.get("pdbx_description") or ""
            entities.append({"entity_id":entity_id,"pdb_description":description,
              "residues":len((pol.get("pdbx_seq_one_letter_code_can") or "").replace("\\n","").replace(" ",""))})
        info=entry.get("rcsb_entry_info") or {}
        resolution=(info.get("resolution_combined") or [None])[0]
        method=";".join(x.get("method","") for x in entry.get("exptl") or [])
        row={**r,"experimental_method":method,"resolution_A":resolution,
             "polymer_entities":json.dumps(entities,ensure_ascii=False),
             "metadata_status":"VERIFIED_PDB_ENTRY_ONLY" if entities and method else "REVIEW"}
        results.append(row)
        print("BIOFILM_RCSB",pdb,method,resolution,entities,flush=True)
    with (OUT/"biofilm_12_rcsb_records.csv").open("w",newline="",encoding="utf-8") as fd:
        writer=csv.DictWriter(fd,fieldnames=list(results[0]));writer.writeheader();writer.writerows(results)
    summary={"n_candidates":12,"n_rcsb_verified":sum(r["metadata_status"]=="VERIFIED_PDB_ENTRY_ONLY" for r in results),
      "not_final_accepted":12,"caution":"No chain-level coordinate completeness, sidechain completeness, species/construct/manual biology, or 11-map InterfaceScout verification performed"}
    (OUT/"SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print("BIOFILM_SCREEN_SUMMARY",json.dumps(summary),flush=True)
    assert summary["n_rcsb_verified"]==12
if __name__=="__main__":main()
