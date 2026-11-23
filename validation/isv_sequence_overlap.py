"""ISV80 protein-sequence redundancy audit against the frozen TV48 panel.

Fetches original RCSB polymer entity sequences, prefilters with shared 4-mers,
then aligns suspected overlaps with Biopython. Reports close sequence similarity,
not complete structural domain-family homology. No score/performance selection.
"""
from __future__ import annotations
import csv,json,re,time,urllib.request,hashlib
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
from itertools import combinations
from Bio.Align import PairwiseAligner
ROOT=Path(__file__).resolve().parents[1]
ISV=list(csv.DictReader((ROOT/"validation/ISV_protein_candidate_pool_80_v2.csv").open()))
TV=list(csv.DictReader((ROOT/"validation/strict_apo_48_selected_candidates_v2_final.csv").open()))
OUT=Path("isv_sequence_similarity");OUT.mkdir(exist_ok=True)
def fetch(url):
    for n in range(3):
        try:
            req=urllib.request.Request(url,headers={"Accept":"application/json","User-Agent":"InterfaceScout-ISV80-sequence-audit/1.0"})
            with urllib.request.urlopen(req,timeout=35) as r:return json.load(r)
        except Exception:
            if n==2:raise
            time.sleep(1+n)
def get(pid):
    entry=fetch(f"https://data.rcsb.org/rest/v1/core/entry/{pid}")
    ids=(entry.get("rcsb_entry_container_identifiers") or {}).get("polymer_entity_ids") or []
    entities=[]
    for id in ids:
        e=fetch(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pid}/{id}")
        poly=e.get("entity_poly") or {}
        if "polypeptide" not in str(poly.get("type","")).lower():continue
        seq=re.sub("[^A-Z]","",str(poly.get("pdbx_seq_one_letter_code_can","")).upper())
        if len(seq)<30:continue
        entities.append({"pdb":pid,"entity_id":str(id),"seq":seq,
          "n_residues":len(seq),"pdb_entity_description":(e.get("rcsb_polymer_entity") or {}).get("pdbx_description","")})
    if not entities:raise ValueError("No >=30 aa polymer entity for "+pid)
    return entities
def grams(s,k=4):return {s[i:i+k] for i in range(len(s)-k+1)}
def candidate_pairs(data):
    # Same entry, same protein: exclude; only distinct PDB records.
    targets=set(r["proposed_pdb_id"] for r in ISV)
    all_pairs=[]
    for i,a in enumerate(data):
        if a["pdb"] not in targets:continue
        ga=a["kmer"]
        for b in data[i+1:]:
            if a["pdb"]==b["pdb"]:continue
            if b["pdb"] not in targets and a["pdb"] not in targets:continue
            if min(len(a["seq"]),len(b["seq"]))<30:continue
            if max(len(a["seq"]),len(b["seq"]))>10000:continue
            gb=b["kmer"]
            score=len(ga&gb)/max(1,len(ga|gb))
            if score>=.01 or a["seq"]==b["seq"]:
                all_pairs.append((a,b,score))
    return all_pairs
def evaluate(item):
    a,b,score=item
    aligner=PairwiseAligner(mode="local",match_score=2,mismatch_score=-1,open_gap_score=-9,extend_gap_score=-1)
    alignment=aligner.align(a["seq"],b["seq"])[0]
    counts=alignment.counts()
    matching=int(counts.identities)
    mismatching=int(counts.mismatches)
    aligned=matching+mismatching
    cov=aligned/min(len(a["seq"]),len(b["seq"])) if min(len(a["seq"]),len(b["seq"])) else 0
    ident=matching/max(1,aligned)
    return {"pdb_a":a["pdb"],"entity_a":a["entity_id"],"description_a":a["pdb_entity_description"],
      "pdb_b":b["pdb"],"entity_b":b["entity_id"],"description_b":b["pdb_entity_description"],
      "a_is_ISV":int(a["pdb"] in {r["proposed_pdb_id"] for r in ISV}),
      "b_is_ISV":int(b["pdb"] in {r["proposed_pdb_id"] for r in ISV}),
      "length_a":len(a["seq"]),"length_b":len(b["seq"]),"fourmer_jaccard":round(score,5),
      "local_identical":matching,"local_mismatches":mismatching,
      "local_aligned":aligned,"local_identity":round(ident,5),
      "coverage_of_shorter":round(cov,5),
      "high_similarity_flag":int((ident>=.40 and cov>=.70 and aligned>=40) or
          (ident>=.80 and cov>=.40 and aligned>=30))}
def main():
    ids=list(dict.fromkeys([r["proposed_pdb_id"] for r in ISV]+[r["pdb_id"] for r in TV]))
    assert len(ids)==128
    found={};errors=[]
    with ThreadPoolExecutor(max_workers=10) as exe:
        pending={exe.submit(get,pdb):pdb for pdb in ids}
        for fut in as_completed(pending):
            pid=pending[fut]
            try:
                result=fut.result();found[pid]=result
                print("ISV_SEQ_SOURCE",pid,len(result),flush=True)
            except Exception as e:
                errors.append({"pdb":pid,"error":str(e)})
                print("ISV_SEQ_ERROR",pid,str(e),flush=True)
    entities=[dict(row,kmer=grams(row["seq"])) for pid in ids for row in found.get(pid,[])]
    (OUT/"source_entity_sequences.json").write_text(json.dumps([
      {k:v for k,v in x.items() if k!="kmer"} for x in entities],indent=2))
    if errors:
        (OUT/"download_errors.json").write_text(json.dumps(errors,indent=2))
        raise RuntimeError(f"{len(errors)} RCSB sequence downloads failed; no final redundancy decision")
    pairs=candidate_pairs(entities)
    print("ISV_SEQ_PREFILTER",len(pairs),flush=True)
    results=[]
    for pair in pairs:
        try:results.append(evaluate(pair))
        except Exception as e:
            raise RuntimeError(f"Alignment failed for {pair[0]['pdb']} and {pair[1]['pdb']}: {e}") from e
    with (OUT/"candidate_sequence_similarity.csv").open("w",newline="",encoding="utf-8") as f:
        fields=["pdb_a","entity_a","description_a","pdb_b","entity_b","description_b",
        "a_is_ISV","b_is_ISV","length_a","length_b","fourmer_jaccard",
        "local_identical","local_mismatches","local_aligned","local_identity",
        "coverage_of_shorter","high_similarity_flag"]
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(results)
    strong=[x for x in results if x["high_similarity_flag"]]
    ids_isv={x["proposed_pdb_id"] for x in ISV}
    ids_tv={x["pdb_id"] for x in TV}
    overlaps_tv={x["pdb_a"] for x in strong if x["pdb_b"] in ids_tv and x["pdb_a"] in ids_isv}
    overlaps_tv|={x["pdb_b"] for x in strong if x["pdb_a"] in ids_tv and x["pdb_b"] in ids_isv}
    overlaps_isv=set()
    for x in strong:
        if x["pdb_a"] in ids_isv and x["pdb_b"] in ids_isv:
            overlaps_isv|={x["pdb_a"],x["pdb_b"]}
    summary={"n_isv":80,"n_tv48":48,"n_pdb_records":len(found),
        "n_protein_entities":len(entities),"fourmer_screened_pairs":len(pairs),
        "high_similarity_entity_pairs":len(strong),
        "isv_ids_with_high_similarity_to_TV48":sorted(overlaps_tv),
        "isv_ids_with_high_similarity_to_ISV":sorted(overlaps_isv),
        "interpretation":"High-similarity flags require source/construct/domain context and do not automatically reject proteins. Low 4-mer similarity does not exclude remote structural homology."}
    (OUT/"SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print("ISV_SEQ_SUMMARY",json.dumps(summary),flush=True)
if __name__=="__main__":main()
