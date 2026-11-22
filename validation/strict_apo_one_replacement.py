"""Find one *new* strict-quality, CATH-unique PDB to replace 1H99.
Does not use InterfaceScout map scores; reuses deterministic screened-pool order.
"""
import csv, json, sys, time, hashlib
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import gemmi
from Bio.Align import PairwiseAligner
from rebuild_strict_apo_48 import cath_map,pick_candidates,inspect,entry_data,homolog,OUT,MAX_WORKERS

ROOT=Path(__file__).resolve().parents[1]
P=Path("strict_apo_one_replacement");P.mkdir(exist_ok=True)
oldids=[r["pdb_id"] for r in csv.DictReader(open(ROOT/"validation"/"strict_48_selected_candidates_v1.csv",encoding="utf-8"))]
rows=list(csv.DictReader(open(ROOT/"validation"/"strict_apo_48_selected_candidates_v1.csv",encoding="utf-8")))
keep=[x for x in rows if x["pdb_id"]!="1H99"]
assert len(keep)==47 and len({x["pdb_id"] for x in keep})==47
used_H={c for x in keep for c in x["CATH_H"].split(";")}
used_T={".".join(c.split(".")[:3]) for c in used_H}
print("KEEP "+json.dumps({"n":len(keep),"nH":len(used_H),"nT":len(used_T)}),flush=True)
cath,src,cathhash=cath_map()
EXPECTED="2f87544bb8844ff37ad1d9348dfe0a23a7c9e4837e0d207733935148e088cc81"
if cathhash!=EXPECTED:raise SystemExit("CATH source changed: screened candidate ordering must be rebuilt")
old=oldids
ordered=list(dict.fromkeys(old+pick_candidates(cath,set(old))))[:34000]
ordered=[pid for pid in ordered if pid in cath]
previous_tested=11477
new_candidates=[pid for pid in ordered[previous_tested:] if pid not in {x["pdb_id"] for x in rows}]
print("POOL "+json.dumps({"new_candidates":len(new_candidates),"previous_tested":previous_tested,"cath_sha256":cathhash}),flush=True)

aligner=PairwiseAligner(mode="local",match_score=2,mismatch_score=-1,open_gap_score=-8,extend_gap_score=-1)
refs=[]
for x in keep:
    pid=x["pdb_id"]
    _,entities,_=entry_data(pid)
    refs.append({"pdb":pid,"entities":entities})
print("REFERENCE_SEQUENCES "+str(len(refs)),flush=True)

def residue_gap_audit(pid,ent):
    raw=(OUT/(pid+"-assembly1.cif")).read_bytes()
    st=gemmi.read_structure(str(OUT/(pid+"-assembly1.cif")))
    chains=[]
    for c in st[0]:
        polymer=list(c.get_polymer())
        if not polymer:continue
        main=c.name.split("-")[0]
        choices=[e for e in ent if main in e["auth_chains"]]
        if len(choices)!=1:raise ValueError("ambiguous entity/chain")
        seq=choices[0]["sequence"]
        observed=[r for r in polymer if r.name in {
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET",
            "PHE","PRO","SER","THR","TRP","TYR","VAL"}]
        labels=[int(r.label_seq) for r in observed]
        if not labels or any(x<=0 for x in labels):raise ValueError("invalid label_seq")
        if any(b-a!=1 for a,b in zip(labels,labels[1:])):raise ValueError("internal missing label_seq residues")
        for r,label in zip(observed,labels):
            if label>len(seq):raise ValueError("observed label_seq beyond declared polymer")
            if gemmi.find_tabulated_residue(r.name).one_letter_code!=seq[label-1]:
                raise ValueError("observed sequence mismatch")
        chains.append(c.name)
    if not chains:raise ValueError("no valid protein chains")
    return {"protein_chains":chains}

counters=Counter();tested=0;selected=None;qualified=0
for k in range(0,len(new_candidates),128):
    batch=new_candidates[k:k+128]
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        tests=list(executor.map(lambda pid:inspect(pid,cath[pid]),batch))
    for r in tests:
        tested+=1
        if r["decision"]!="QUALIFIED":
            counters[r["reason"].split(":")[0] if r.get("reason") else "UNKNOWN"]+=1
            continue
        qualified+=1
        if set(r["CATH_H"]) & used_H:
            counters["CATH_H_COLLISION"]+=1
            continue
        topologies={".".join(x.split(".")[:3]) for x in r["CATH_H"]}
        if topologies & used_T:
            counters["CATH_T_COLLISION"]+=1
            continue
        if any(homolog(r,previous,aligner) for previous in refs):
            counters["SEQUENCE_SIMILARITY"]+=1
            continue
        try:
            r["residue_level_check"]=residue_gap_audit(r["pdb"],r["entities"])
        except Exception as e:
            counters["RESIDUE_AUDIT_"+str(e)[:50]]+=1
            continue
        selected=r
        print("REPLACEMENT_FOUND "+json.dumps({"pdb":r["pdb"],"resolution":r["resolution_A"],
              "rfree":r["rfree"],"CATH_H":r["CATH_H"],"n_chains":r["n_chains"],
              "CATH_T":sorted(topologies),"assembly_sha256":r["assembly_sha256"]}),flush=True)
        break
    print("REPLACEMENT_PROGRESS "+json.dumps({"tested":tested,"qualified":qualified,"failures":counters.most_common(8)}),flush=True)
    if selected is not None:break
result={"replacing":"1H99","candidate":selected,"new_candidates_tested":tested,"qualified_tested":qualified,
        "selection_without_interfacescout_scores":True,"CATH_hash":cathhash,"failures":dict(counters)}
(P/"replacement.json").write_text(json.dumps(result,indent=2))
if selected is None:
    print("REPLACEMENT_SUMMARY "+json.dumps({"status":"NOT_FOUND","tested":tested,"qualified":qualified}),flush=True)
    sys.exit(1)
print("REPLACEMENT_SUMMARY "+json.dumps({"status":"FOUND","pdb":selected["pdb"],"tested":tested,"qualified":qualified}),flush=True)
