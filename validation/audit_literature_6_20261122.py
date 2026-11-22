"""Independent candidate QC for experimentally localized InterfaceScout literature cases.
PASS means strict structure/sequence/CATH criteria only, NEVER experimental validation.
No structure repair. Retain rejected candidates, untested fields, and accession provenance.
"""
import csv, hashlib, json, re, sys, time, traceback
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import requests
from Bio.Align import PairwiseAligner
import rebuild_strict_apo_48 as base

ROOT=Path(__file__).resolve().parents[1]
OUT=Path("literature_candidate_qa_20261122"); OUT.mkdir(exist_ok=True)
LITERATURE_23="4F5S 2VUF 2CBA 1UBQ 1JNJ 3ECA 2HAV 1TGU 3NWV 7RSA 3LN2 1LYZ 1HRC 1MBN 2HHB 2PTN 1A5R 2Y0G 256L 1AML 2OED 2GB1 7KWI".split()
CANDIDATES={
"2VH7":{"protein":"human acylphosphatase","experimental_doi":"10.1021/la401985d","original":True},
"6O5E":{"protein":"vitronectin HX","experimental_doi":"10.1073/pnas.2007699117","original":True},
"7TXR":{"protein":"vitronectin HX alternate wildtype calcium","experimental_doi":"10.1073/pnas.2007699117","original":False},
"7U68":{"protein":"vitronectin HX alternate wildtype sodium","experimental_doi":"10.1073/pnas.2007699117","original":False},
"5J2R":{"protein":"BRAF RBD drug-bound","experimental_doi":"10.1021/jacs.4c18380","original":True},
"5J17":{"protein":"BRAF RBD unbound NMR","experimental_doi":"10.1021/jacs.4c18380","original":False},
"1F6R":{"protein":"apo bovine alpha-lactalbumin","experimental_doi":"10.1073/pnas.0401603101","original":True},
"3STN":{"protein":"apo human liver FABP","experimental_doi":"10.1021/ja507310m","original":True},
"2L67":{"protein":"human liver FABP apo NMR alternate","experimental_doi":"10.1021/ja507310m","original":False},
"8ONG":{"protein":"human EMAP II","experimental_doi":"10.3390/ijms27020605","original":True}
}
DEV=list(csv.DictReader(open(ROOT/"validation"/"strict_apo_48_selected_candidates_v2_final.csv")))
DEV_IDS=[r["pdb_id"] for r in DEV]
assert len(DEV)==48 and len(set(DEV_IDS))==48 and len(LITERATURE_23)==23
HEADERS={"User-Agent":"InterfaceScout research reproducibility audit 2026-11-22"}
def get(pid):
    e=base.entry_data(pid)
    return e
def load_entity(pid):
    e=base.jget(f"https://data.rcsb.org/rest/v1/core/entry/{pid}")
    seqs=[]
    for n in e.get("rcsb_entry_container_identifiers",{}).get("polymer_entity_ids",[]):
        p=base.jget(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pid}/{n}")
        z=p.get("entity_poly",{})
        if "polypeptide" not in z.get("type","").lower():continue
        seq=re.sub("[^A-Z]","",z.get("pdbx_seq_one_letter_code_can","").upper())
        if not seq:continue
        seqs.append({"entity":str(n),"sequence":seq,"length":len(seq),
          "mutation_count":p.get("rcsb_polymer_entity",{}).get("pdbx_mutation",None),
          "description":p.get("rcsb_polymer_entity",{}).get("pdbx_description",None),
          "auth_chains":p.get("rcsb_polymer_entity_container_identifiers",{}).get("auth_asym_ids",[])})
    return {"pdb":pid,"entities":seqs,
      "method":[x.get("method") for x in e.get("exptl",[])],
      "resolution":(e.get("rcsb_entry_info",{}).get("resolution_combined") or [None])[0],
      "title":e.get("struct",{}).get("title",None),
      "n_polymer_entities":len(seqs),
      "rcsb_entry":f"https://www.rcsb.org/structure/{pid}"}
def similarity(a,b,aligner):
    best=None
    for ea in a["entities"]:
        for eb in b["entities"]:
            sa,sb=ea["sequence"],eb["sequence"]
            if not sa or not sb:continue
            al=aligner.align(sa,sb)[0]
            paired=[(int(i),int(j)) for i,j in zip(*al.indices) if i>=0 and j>=0]
            n=len(paired)
            if not n:continue
            same=sum(sa[i]==sb[j] for i,j in paired)
            ident=same/n
            cov=n/min(len(sa),len(sb))
            metric={"identity":round(ident,5),"shorter_coverage":round(cov,5),"aligned":n,
                "identical":same,"candidate_entity":ea["entity"],"reference_entity":eb["entity"],
                "candidate_len":len(sa),"reference_len":len(sb)}
            if best is None or (ident*cov, n)>(best["identity"]*best["shorter_coverage"],best["aligned"]):
                best=metric
    return best
def main():
    # Use same CATH release/wwPDB thresholds and coordinate parser as frozen development audit.
    cath,cat_url,cat_hash=base.cath_map()
    ids=list(dict.fromkeys(DEV_IDS+LITERATURE_23+list(CANDIDATES)))
    with ThreadPoolExecutor(max_workers=10) as pool:
        answers=list(pool.map(lambda pid:(pid,load_entity(pid)),ids))
    data={pid:x for pid,x in answers}
    aligner=PairwiseAligner(mode="local",match_score=2,mismatch_score=-1,open_gap_score=-8,extend_gap_score=-1)
    results=[]
    for k,(pid,meta) in enumerate(CANDIDATES.items(),1):
        r={"pdb":pid,**meta,**data[pid],"CATH_H":sorted(cath.get(pid,[])),
           "candidate_CATH_known":bool(cath.get(pid))}
        res=base.inspect(pid,cath.get(pid,[]))
        r["strict_xray_qa"]=res.get("decision")
        r["qa_reason"]=res.get("reason")
        for field in ["wwpdb","rfree","resolution_A","assembly_sha256","assembly_chains","n_chains","n_residues"]:
            if field in res:r[field]=res[field]
        def topmatch(refs,n=8):
            arr=[]
            for other in refs:
                z=similarity(data[pid],data[other],aligner)
                if z:arr.append({"pdb":other,**z,
                    "shared_CATH_H":sorted(set(cath.get(pid,[])) & set(cath.get(other,[]))),
                    "shared_CATH_T":sorted({".".join(c.split(".")[:3]) for c in cath.get(pid,[])} & {".".join(c.split(".")[:3]) for c in cath.get(other,[])})})
            arr.sort(key=lambda x:(x["identity"]*x["shorter_coverage"],x["aligned"]),reverse=True)
            return arr[:n],arr
        r["closest_development"],all_dev=topmatch(DEV_IDS)
        r["closest_existing_literature"],all_lit=topmatch(LITERATURE_23)
        r["development_sequence_conflicts"]=[x for x in all_dev if x["identity"]>=.30 and x["shorter_coverage"]>=.70 and x["aligned"]>=40]
        r["literature_sequence_conflicts"]=[x for x in all_lit if x["identity"]>=.30 and x["shorter_coverage"]>=.70 and x["aligned"]>=40]
        r["development_CATH_H_conflicts"]=[x["pdb"] for x in all_dev if x["shared_CATH_H"]]
        r["literature_CATH_H_conflicts"]=[x["pdb"] for x in all_lit if x["shared_CATH_H"]]
        r["development_CATH_T_conflicts"]=[x["pdb"] for x in all_dev if x["shared_CATH_T"]]
        r["literature_CATH_T_conflicts"]=[x["pdb"] for x in all_lit if x["shared_CATH_T"]]
        r["unique_PDB_vs_48_and_23"]=pid not in DEV_IDS and pid not in LITERATURE_23
        r["passed_independence_sequence_to_development"]=not r["development_sequence_conflicts"]
        r["passed_independence_sequence_to_old_23"]=not r["literature_sequence_conflicts"]
        r["passed_CATH_H_to_development"]=r["candidate_CATH_known"] and not r["development_CATH_H_conflicts"]
        r["passed_CATH_T_to_development"]=r["candidate_CATH_known"] and not r["development_CATH_T_conflicts"]
        r["passed_structure_and_development_independence"]=(res.get("decision")=="QUALIFIED"
              and r["passed_independence_sequence_to_development"] and r["passed_CATH_H_to_development"])
        results.append(r)
        print(f"CHECK {k}/{len(CANDIDATES)} {pid} qa={r['strict_xray_qa']} reason={r['qa_reason']} "
          f"dev_seq_conflicts={len(r['development_sequence_conflicts'])} old23_conflicts={len(r['literature_sequence_conflicts'])}"
          f" CATH_dev={r['development_CATH_H_conflicts'][:4]}",flush=True)
    (OUT/"full_candidate_qa.json").write_text(json.dumps(results,indent=2,ensure_ascii=False))
    fields=["pdb","protein","experimental_doi","original","method","resolution","rfree","strict_xray_qa","qa_reason",
      "candidate_CATH_known","CATH_H","n_chains","n_residues","unique_PDB_vs_48_and_23",
      "passed_independence_sequence_to_development","passed_independence_sequence_to_old_23",
      "passed_CATH_H_to_development","passed_CATH_T_to_development","passed_structure_and_development_independence",
      "development_CATH_H_conflicts","literature_CATH_H_conflicts",
      "development_CATH_T_conflicts","literature_CATH_T_conflicts",
      "development_sequence_conflicts","literature_sequence_conflicts",
      "closest_development","closest_existing_literature","assembly_sha256"]
    with (OUT/"candidate_qa_summary.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for r in results:
            x={k:r.get(k) for k in fields}
            for k in fields:
                if isinstance(x[k],(dict,list)):x[k]=json.dumps(x[k],ensure_ascii=False)
            w.writerow(x)
    (OUT/"methods.json").write_text(json.dumps({
      "development_n":len(DEV_IDS),"literature_n":len(LITERATURE_23),"candidate_n":len(CANDIDATES),
      "cath_source":cat_url,"cath_sha256":cat_hash,
      "strict_gates":base.GATES,
      "sequence_alignment":{"method":"Biopython PairwiseAligner local","match":2,"mismatch":-1,"gap_open":-8,"gap_extend":-1,
             "flag":"identity >= 0.30 and alignment coverage >= 0.70 of shorter sequence and aligned >= 40"},
      "caution":"Geometry pass is not a paper-level ground-truth pass; experimental residue mapping must be audited separately.",
      "reference_48_ids":DEV_IDS,"reference_23_ids":LITERATURE_23},indent=2))
    print("RESULTS "+json.dumps({"candidates":len(results),
        "xray_qa_pass":[x["pdb"] for x in results if x["strict_xray_qa"]=="QUALIFIED"],
        "structural_development_independence_pass":[x["pdb"] for x in results if x["passed_structure_and_development_independence"]],
        "n_errors":sum(bool(x.get("qa_reason") and "connection" in str(x["qa_reason"]).lower()) for x in results)}),flush=True)
if __name__=="__main__":main()
