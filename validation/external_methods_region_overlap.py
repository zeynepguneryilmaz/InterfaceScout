"""Descriptive PEP-Patch vs InterfaceScout region-set agreement, no accuracy labels.

Only compare SINGLE-chain structures with uniquely mappable residue IDs.
Multi-chain PEP-Patch residue reports omit chain identity and cannot be safely
paired to original InterfaceScout keys; record exclusion transparently.
"""
from pathlib import Path
import json
import re
import pandas as pd
import numpy as np

IS=Path("isv76_final_model")
EXT=Path("external_methods_full_panel/pep_patch")
OUT=Path("external_methods_full_panel/final")
OUT.mkdir(parents=True,exist_ok=True)
summary=Path("external_methods_full_panel/final/SUMMARY.json")
if not summary.exists():
    raise ValueError("Run full-panel strict aggregation first")

case=[]
for p in sorted(IS.joinpath("proteins").glob("*/canonical_maps_complete.json")):
    pid=p.parent.name
    own=json.loads(p.read_text())
    profiles=own["residue_profiles"]
    exposed={r["key"] for r in profiles if r.get("surface_exposed",False)}
    chains={r["chain"] for r in profiles}
    result={"pdb_id":pid,"n_exposed":len(exposed),"n_chains":len(chains)}
    if len(chains)!=1:
        result.update(status="NOT_COMPARABLE_MULTICHAIN",reason="PEP-Patch patch-members CSV has no chain labels")
        case.append(result);continue
    prepared=EXT/(pid+"_result.json")
    if not prepared.exists():
        result.update(status="MISSING_PEP_RESULT");case.append(result);continue
    prep=json.loads(prepared.read_text())
    mode={v["task"]:v["runtime"]["status"] for v in prep.get("attempts",[])}
    if mode.get("electrostatic")!="PASS":
        result.update(status="PEP_ELECTROSTATIC_FAILED");case.append(result);continue
    patches_file=EXT/pid/"electrostatic"/"patches.csv"
    residues_file=EXT/pid/"electrostatic"/"residues.csv"
    if not patches_file.is_file() or not residues_file.is_file():
        result.update(status="MISSING_PEP_PATCH_EXPORT");case.append(result);continue
    patches=pd.read_csv(patches_file)
    residues=pd.read_csv(residues_file)
    lookup=dict(zip(residues["nr"],residues["residues"]))
    chain=next(iter(chains))
    for sign,channel in (("positive","cationic"),("negative","anionic")):
        pep=set();unmatched=[]
        for nr in patches.loc[patches.type==sign].sort_values("area",ascending=False).head(5)["nr"]:
            for token in str(lookup.get(nr,"")).split():
                m=re.fullmatch(r"([A-Z]{3})(-?\d+)([A-Za-z]?)",token)
                if not m:
                    unmatched.append(token);continue
                key=f"{chain}:{int(m.group(2))}:{m.group(3)}"
                if key not in {q["key"] for q in profiles} or next(
                    q["res_name"] for q in profiles if q["key"]==key)!=m.group(1):
                    unmatched.append(token);continue
                pep.add(key)
        ism={str(x) for patch in own["maps"][channel]["patches"][:5]
             for x in patch.get("members",[])}
        den=pep|ism
        result.update({
            sign+"_n_PEP_members":len(pep),
            sign+"_n_InterfaceScout_members":len(ism),
            sign+"_n_common_members":len(pep&ism),
            sign+"_jaccard":len(pep&ism)/len(den) if den else None,
            sign+"_PEP_exposed_coverage":len(pep&exposed)/max(len(exposed),1),
            sign+"_InterfaceScout_exposed_coverage":len(ism&exposed)/max(len(exposed),1),
            sign+"_unmatched_tokens":";".join(unmatched[:30]),
        })
    if result["positive_unmatched_tokens"] or result["negative_unmatched_tokens"]:
        result["status"]="PARTIAL_RESIDUE_IDENTITY"
    else:
        result["status"]="VALID_DESCRIPTIVE_COMPARISON"
    case.append(result)

if len(case)!=76 or len({x["pdb_id"] for x in case})!=76:
    raise ValueError("Not all active protein identities accounted for")
df=pd.DataFrame(case)
df.to_csv(OUT/"PEP_PATCH_VS_INTERFACESCOUT_CHARGED_TOP5_76_QA.csv",index=False)
eligible=df[df.status=="VALID_DESCRIPTIVE_COMPARISON"]
report={
  "n_total":len(case),
  "status_counts":{str(k):int(v) for k,v in df.status.value_counts().items()},
  "n_valid":len(eligible),
  "median_positive_Jaccard":(float(eligible.positive_jaccard.median()) if len(eligible) else None),
  "median_negative_Jaccard":(float(eligible.negative_jaccard.median()) if len(eligible) else None),
  "limitation":"PEP mesh electrostatic positive/negative regions and InterfaceScout residue patches are different definitions; no experimental truth labels",
  "PDB_matching":"Only verified unique single-chain residue names/numbering; missing insertion/chain IDs not guessed",
}
(OUT/"CHARGED_REGION_OVERLAP_SUMMARY.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
print("PEP_PATCH_REGION_OVERLAP "+json.dumps(report),flush=True)
