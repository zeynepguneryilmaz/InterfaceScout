"""RCSB mmCIF coordinate-level pre-screen for ISV candidate pool (not final QA).

For each deposited polymer chain in model 1, count unique CA positions divided
by the deposited polymer sequence length. Does not rebuild missing sidechains,
nor establish biological assembly or source publication-to-PDB correspondence.
"""
from __future__ import annotations
import csv,json,urllib.request,time,collections
from pathlib import Path
import gemmi
ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/"validation/ISV_protein_candidate_pool_80_v2.csv"
OUT=Path("isv_coordinate_screen");OUT.mkdir(exist_ok=True)
def fetch(url):
    for attempt in range(3):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"InterfaceScout-ISV-Coordinate-Audit/1.0"})
            with urllib.request.urlopen(req,timeout=50) as response:return response.read()
        except Exception:
            if attempt==2:raise
            time.sleep(attempt+1)
def get_entity(pdb,eid):
    return json.loads(fetch(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{eid}"))
def scan(entry):
    pdb=entry["proposed_pdb_id"]
    cif_text=fetch(f"https://files.rcsb.org/download/{pdb}.cif").decode("utf8")
    block=gemmi.cif.read_string(cif_text).sole_block()
    def v(k):return list(block.find_values("_atom_site."+k))
    ents=v("label_entity_id")
    seq=v("label_seq_id");atoms=v("label_atom_id");chains=v("label_asym_id")
    models=v("pdbx_PDB_model_num")
    if not (len(ents)==len(seq)==len(atoms)==len(chains)==len(models)):
        raise ValueError("atom-site column length mismatch")
    ids=sorted(set(en for en,chain in by_chain))
    if not ids:raise ValueError("no numeric protein entity labels")
    by_chain=collections.defaultdict(set)
    for en,se,atom,chain,model in zip(ents,seq,atoms,chains,models):
        if model=="1" and atom=="CA" and se.isdigit():
            by_chain[(en,chain)].add(int(se))
    methods=block.find_values("_exptl.method")
    method=";".join(str(x) for x in methods)
    res_values=[str(x) for x in block.find_values("_refine.ls_d_res_high") if str(x) not in ("?",".")]
    resolution=next((float(x) for x in res_values if x),None)
    ids=sorted({en for en,_ in by_chain})
    records=[]
    for eid in ids:
        ent=get_entity(pdb,eid)
        ep=ent.get("entity_poly",{})
        if "polypeptide" not in str(ep.get("type","")).lower():continue
        length=len(str(ep.get("pdbx_seq_one_letter_code_can","")).replace("\n","").replace(" ",""))
        label=(ent.get("rcsb_polymer_entity") or {}).get("pdbx_description","")
        for (en,chain),observed in by_chain.items():
            if en!=eid:continue
            records.append({"entity_id":eid,"label_asym_chain":chain,"protein_name_rcsb":label,
                "seq_length":length,"n_observed_CA":len(observed),
                "CA_coverage":round(len(observed)/length,5) if length else None,
                "n_internal_missing_seq_positions":sum(
                    k not in observed for k in range(min(observed),max(observed)+1)) if observed else 0})
    if not records:raise ValueError("no polymer chain CA coordinates")
    best=max(records,key=lambda r:(r["CA_coverage"] or 0,r["n_observed_CA"]))
    method=(method or "NMR/unknown-see-record")
    return records,best,method,resolution
def main():
    rows=list(csv.DictReader(INPUT.open()))
    assert len(rows)==80
    out=[];errors=[]
    for r in rows:
        pid=r["proposed_pdb_id"]
        try:
            chains,best,method,res=scan(r)
            coverage=float(best["CA_coverage"] or 0)
            # Coordinate-only screen: do not conflate with structural full QA.
            coord_pass=(coverage>=.9 and best["n_observed_CA"]>=50)
            row={**r,"coordinate_status":"CA_COVERAGE_PASS_PROVISIONAL" if coord_pass else "CA_COVERAGE_REVIEW",
                 "best_chain":best["label_asym_chain"],"best_chain_CA_coverage":coverage,
                 "best_chain_n_observed_CA":best["n_observed_CA"],"best_chain_seq_length":best["seq_length"],
                 "best_chain_n_internal_missing":best["n_internal_missing_seq_positions"],
                 "all_protein_chains_json":json.dumps(chains),"error":""}
            print("ISV_CA",pid,row["coordinate_status"],best["n_observed_CA"],best["seq_length"],coverage,flush=True)
        except Exception as e:
            row={**r,"coordinate_status":"CA_UNVERIFIED","best_chain":"",
                 "best_chain_CA_coverage":"","best_chain_n_observed_CA":"",
                 "best_chain_seq_length":"","best_chain_n_internal_missing":"",
                 "all_protein_chains_json":"[]","error":str(e)}
            errors.append((pid,str(e)))
            print("ISV_CA_ERROR",pid,str(e),flush=True)
        out.append(row)
    with (OUT/"isv_coordinate_80.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
    counts=dict(collections.Counter(x["coordinate_status"] for x in out))
    summary={"n_pdb":80,"counts":counts,"errors":errors,
      "scope":"Deposited mmCIF model-1 C-alpha coverage relative to deposited polymer construct only. Not full UniProt coverage, experimental quality, sidechain completeness, assembly selection, or acceptance."}
    (OUT/"SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print("ISV_CA_SUMMARY",json.dumps(summary),flush=True)
    if errors:raise SystemExit(2)
if __name__=="__main__":main()
