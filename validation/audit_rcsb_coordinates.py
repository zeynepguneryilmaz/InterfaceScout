"""Coordinate-aware RCSB audit; metadata screening, not final scientific acceptance."""
import csv, json, os, re, sys, tempfile, time
from pathlib import Path
from itertools import combinations
import requests
import gemmi
from Bio.Align import PairwiseAligner

IDS="1CRN 1R69 1SHG 2PPN 4AKE 1OMP 1TIM 5CSC 1UBQ 1PGA 1BPI 1HEL 1A6M 1HRC 1GFL 1TEN 1BRS 1A3N 1ALD 1GKY 1CDG 1CSE 1A4Y 1TCA 1ATP 1PPB 1HVR 1LDM 1J8U 1CLL 1FNF 1B0L 1PPA 1BTL 1RNH 1E9I 1MSO 1CSP 1DPS 1PLC 1MUP 1RBL 2TRX 1SMD 1HML 1THW 1BEB 1CON".split()
ADDITIONAL="1LYZ 1AKI 1AON 1FAT 1F8E 2PTC 1G6X 1A2P 1TUP 1AKE 1AKE 1BTA 1HZH 1GCI 1CYO 1A0I 1GHA 1TGS 1C5E 1PDB 1ECA 1PPE 1RGS 1A1X".split()
# Additional codes are a search pool ONLY; no automatic approval.
ALL=list(dict.fromkeys(IDS+ADDITIONAL))
s=requests.Session();s.headers["User-Agent"]="InterfaceScout academic validation 1.0 contact via GitHub"
base="https://data.rcsb.org/rest/v1/core"
out=Path("audit_results");out.mkdir(exist_ok=True)
def get(url):
    for i in range(3):
        try:
            r=s.get(url,timeout=35)
            if r.status_code==404:return None
            r.raise_for_status();return r
        except requests.RequestException as e:
            if i==2:return {"error":str(e)}
            time.sleep(i+1)
def api(path):
    x=get(base+path)
    if isinstance(x,dict) or x is None:return x
    return x.json()
def first(x,*keys):
    for k in keys:
        if isinstance(x,dict):x=x.get(k)
        else:return None
    return x
def cif_summary(pid,assembly=False):
    suffix="-assembly1.cif" if assembly else ".cif"
    resp=get("https://files.rcsb.org/download/"+pid+suffix)
    if not resp or isinstance(resp,dict):return {"error":resp}
    name=out/(pid+suffix)
    name.write_bytes(resp.content)
    try:
        st=gemmi.read_structure(str(name))
        if len(st)==0:return {"error":"no model"}
        m=st[0]; chains=[]; models=len(st)
        for ch in m:
            polymer=ch.get_polymer()
            if len(polymer)==0:continue
            residues=[x for x in polymer if x.name in ("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL")]
            count_ca=sum(any(a.name=="CA" for a in residue) for residue in residues)
            chains.append({"chain":ch.name,"residues":len(residues),"ca":count_ca})
        return {"models":models,"polymer_chains":len(chains),"chains":chains,"total_ca":sum(c["ca"] for c in chains),"bytes":len(resp.content)}
    except Exception as e:return {"error":str(e)}
    finally:
        if name.exists():name.unlink()
def entity_sequence(pid,ent):
    q=api("/polymer_entity/"+pid+"/"+str(ent))
    if not q or "error" in q:return None
    seq=first(q,"entity_poly","pdbx_seq_one_letter_code_can") or ""
    seq=re.sub(r"[^A-Z]","",seq.upper())
    return {"entity":str(ent),"seq":seq,"length":len(seq),"auth_chains":first(q,"rcsb_polymer_entity_container_identifiers","auth_asym_ids") or [],"description":first(q,"rcsb_polymer_entity","pdbx_description") or ""}
records=[]; seqs=[]
for k,pid in enumerate(ALL,1):
    e=api("/entry/"+pid)
    if not e or "error" in e:
        print("FAILED_ENTRY",pid,str(e),flush=True)
        records.append({"pdb":pid,"in_original_48":pid in IDS,"error":str(e)});continue
    rr=first(e,"rcsb_entry_info","resolution_combined") or []
    res=rr[0] if rr else None
    meth=";".join(x.get("method","") for x in e.get("exptl",[]))
    refine=e.get("refine") or []
    rfree=next((x.get("ls_R_factor_R_free") for x in refine if x.get("ls_R_factor_R_free") is not None),None)
    assembly=api("/assembly/"+pid+"/1")
    assembler=first(assembly,"pdbx_struct_assembly","oligomeric_details") or "unknown"
    peids=first(e,"rcsb_entry_container_identifiers","polymer_entity_ids") or []
    q=[entity_sequence(pid,ent) for ent in peids]
    q=[x for x in q if x and x["length"]]
    for z in q:seqs.append((pid,z))
    orig=cif_summary(pid)
    ass=cif_summary(pid,True)
    record={"pdb":pid,"in_original_48":pid in IDS,"method":meth,"resolution_A":res,"r_free":rfree,
        "assembly_annotation":assembler,"entities":q,"asymmetric_unit":orig,"assembly1":ass,
        "review_flags":[]}
    flags=record["review_flags"]
    if meth!="X-RAY DIFFRACTION":flags.append("NON_XRAY")
    if res is None or res>2.5:flags.append("RESOLUTION_GT_2_5")
    if rfree is None:flags.append("RFREE_NOT_REPORTED")
    elif rfree>0.30:flags.append("RFREE_GT_0_30")
    elif rfree>0.28:flags.append("RFREE_GT_0_28")
    if any(z["length"]<40 for z in q):flags.append("SHORT_POLYMER_ENTITY")
    if len(q)>1:flags.append("MULTI_PROTEIN_ENTITY")
    if any("error" in x for x in [orig,ass]):flags.append("COORDINATE_FETCH_OR_PARSE_FAILURE")
    if "error" not in orig and "error" not in ass:
        if orig["total_ca"]==0 or ass["total_ca"]==0:flags.append("NO_CA")
        if ass["polymer_chains"]>4:flags.append("LARGE_OLIGOMER")
        elif ass["polymer_chains"]>1:flags.append("OLIGOMER")
        # Unresolved residue count only interpretable where one entity/chain and no domain fusion
        if len(q)==1 and orig["polymer_chains"]==1:
            ratio=orig["total_ca"]/q[0]["length"]
            record["observed_ca_fraction"]=round(ratio,3)
            if ratio<0.9:flags.append("LOW_COORDINATE_COVERAGE")
    records.append(record)
    print("RECORD",json.dumps({x:record[x] for x in ["pdb","in_original_48","resolution_A","r_free","assembly_annotation","observed_ca_fraction","review_flags"] if x in record},separators=(",",":")),flush=True)
print("PAIRWISE_SEQUENCE_START",len(seqs),flush=True)
aligner=PairwiseAligner(mode="local",match_score=2,mismatch_score=-1,open_gap_score=-8,extend_gap_score=-1)
hits=[]
for (id1,a),(id2,b) in combinations(seqs,2):
    if id1==id2:continue
    if not a["seq"] or not b["seq"]:continue
    alignment=aligner.align(a["seq"],b["seq"])[0]
    aa,bb=alignment.indices
    aligned=[(int(i),int(j)) for i,j in zip(aa,bb) if i>=0 and j>=0]
    if not aligned:continue
    ident=sum(a["seq"][i]==b["seq"][j] for i,j in aligned)/len(aligned)
    coverage=len(aligned)/min(len(a["seq"]),len(b["seq"]))
    if ident>=0.30 and coverage>=0.70:
        hits.append({"pdb1":id1,"pdb2":id2,"entity1":a["entity"],"entity2":b["entity"],"identity":round(ident,4),"coverage_shorter":round(coverage,4),"aligned_positions":len(aligned)})
for h in sorted(hits,key=lambda x:(-x["identity"],-x["coverage_shorter"])):
    print("HOMOLOGY_REVIEW",json.dumps(h,separators=(",",":")),flush=True)
summary={"original48":len(IDS),"audited_original":sum(r.get("in_original_48") and "error" not in r for r in records),
 "original_resolution_le_2_5":sum(r.get("in_original_48") and isinstance(r.get("resolution_A"),(int,float)) and r["resolution_A"]<=2.5 for r in records),
 "original_coordinate_fetch_success":sum(r.get("in_original_48") and "error" not in r.get("asymmetric_unit",{}) and "error" not in r.get("assembly1",{}) for r in records if "error" not in r),
 "pairwise_homology_review_count":sum(h["pdb1"] in IDS and h["pdb2"] in IDS for h in hits),
 "original48_flags":{r["pdb"]:r.get("review_flags",[]) for r in records if r.get("in_original_48")},
 "not_final_validation":True}
(out/"full_audit.json").write_text(json.dumps(records,indent=2))
(out/"homology_review.json").write_text(json.dumps(hits,indent=2))
(out/"summary.json").write_text(json.dumps(summary,indent=2))
print("FINAL_SUMMARY",json.dumps(summary,separators=(",",":")),flush=True)
