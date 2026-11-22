"""Select quality-clean, CATH-nonredundant PDB entries without InterfaceScout scores.

This program is prospective selection, not model fitting. It reads RCSB, CATH
and wwPDB primary records; fails closed if standards are not met. No atom repair.
"""
import csv, gzip, hashlib, json, random, re, string, sys, time, traceback, requests
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from xml.etree import ElementTree as ET
import gemmi
from Bio.Align import PairwiseAligner

ROOT=Path(__file__).resolve().parents[1]
OUT=Path("strict_panel_48")
OUT.mkdir(exist_ok=True)
PROTEIN={"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
SIDE={"ALA":"CB","ARG":"CB CG CD NE CZ NH1 NH2","ASN":"CB CG OD1 ND2","ASP":"CB CG OD1 OD2","CYS":"CB SG","GLN":"CB CG CD OE1 NE2","GLU":"CB CG CD OE1 OE2","GLY":"","HIS":"CB CG ND1 CD2 CE1 NE2","ILE":"CB CG1 CG2 CD1","LEU":"CB CG CD1 CD2","LYS":"CB CG CD CE NZ","MET":"CB CG SD CE","PHE":"CB CG CD1 CD2 CE1 CE2 CZ","PRO":"CB CG CD","SER":"CB OG","THR":"CB OG1 CG2","TRP":"CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2","TYR":"CB CG CD1 CD2 CE1 CE2 CZ OH","VAL":"CB CG1 CG2"}
SEED=20261122
TARGET=48
DESIRED_PASS_POOL=65
BATCH=64
MAX_CANDIDATES=2800
MAX_WORKERS=18
# Predeclared scientific-quality gates (all satisfied; no repair)
GATES={"max_resolution_A":2.50,"max_rfree":0.280,"max_clashscore":15.0,
       "max_rama_outliers_pct":1.0,"max_rota_outliers_pct":5.0,
       "max_rsrz_outliers_pct":5.0,"min_CA_coverage":0.98,
       "max_missing_sidechain_residues":0,"max_altloc_residues":0,
       "min_protein_length":70,"max_protein_length":650,
       "max_assembly_chains":4,"max_assembly_residues":1500}
HEADERS={"User-Agent":"InterfaceScout nonredundant structural QA scholarly software (2026-11)"}
def fetch(url,tries=2,timeout=25):
    exc=None
    for n in range(tries):
        try:
            r=requests.get(url,timeout=timeout,headers=HEADERS)
            r.raise_for_status()
            return r.content
        except Exception as e:
            exc=e
            time.sleep(0.5+n)
    raise RuntimeError(str(exc))
def jget(url):return json.loads(fetch(url))
def parse_number(v):
    try:
        f=float(v)
        return f if f==f and abs(f)<1e10 else None
    except (TypeError,ValueError):return None
def cath_map():
    url="https://download.cathdb.info/cath/releases/latest-release/cath-classification-data/cath-domain-list.txt"
    content=fetch(url,tries=3,timeout=100)
    mapping=defaultdict(set)
    for line in content.decode("utf-8","replace").splitlines():
        fields=line.split()
        if len(fields)<5 or len(fields[0])<7 or not fields[0][0].isdigit():continue
        mapping[fields[0][:4].upper()].add(".".join(fields[1:5]))
    return dict(mapping),url,hashlib.sha256(content).hexdigest()
def wwpdb(pid):
    sub=pid.lower()
    url=f"https://files.rcsb.org/pub/pdb/validation_reports/{sub[1:3]}/{sub}/{sub}_validation.xml.gz"
    raw=fetch(url)
    if raw[:2]==b"\x1f\x8b":raw=gzip.decompress(raw)
    root=ET.fromstring(raw)
    en=next((e for e in root.iter() if e.tag.split("}")[-1]=="Entry"),None)
    if en is None:raise ValueError("wwPDB entry node absent")
    d=en.attrib
    metrics={k:parse_number(d.get(k)) for k in ("clashscore","percent-rama-outliers","percent-rota-outliers","percent-RSRZ-outliers","DCC_Rfree","PDB-Rfree")}
    rfree=metrics["DCC_Rfree"] if metrics["DCC_Rfree"] is not None else metrics["PDB-Rfree"]
    if rfree is None:raise ValueError("Rfree missing")
    for key,lim in [("clashscore",GATES["max_clashscore"]),("percent-rama-outliers",GATES["max_rama_outliers_pct"]),("percent-rota-outliers",GATES["max_rota_outliers_pct"])]:
        if metrics[key] is None:raise ValueError("wwPDB metric missing "+key)
        if metrics[key]>lim:raise ValueError("wwPDB threshold "+key)
    if metrics["percent-RSRZ-outliers"] is not None and metrics["percent-RSRZ-outliers"]>GATES["max_rsrz_outliers_pct"]:
        raise ValueError("wwPDB threshold RSRZ")
    if rfree>GATES["max_rfree"]:raise ValueError("wwPDB threshold Rfree")
    return metrics,rfree,url
def entry_data(pid):
    e=jget(f"https://data.rcsb.org/rest/v1/core/entry/{pid}")
    if [x.get("method") for x in e.get("exptl",[])]!=["X-RAY DIFFRACTION"]:raise ValueError("not X-ray")
    resolution=e.get("rcsb_entry_info",{}).get("resolution_combined") or []
    if not resolution or resolution[0]>GATES["max_resolution_A"]:raise ValueError("resolution")
    title=e.get("struct",{}).get("title","").lower()
    blockers=("membrane protein","transmembrane","rhodopsin","photosystem","integral membrane","membrane-bound","membrane embedded","gpcr","ion channel")
    if any(x in title for x in blockers):raise ValueError("membrane-related entry")
    entities=[]
    for ent in e.get("rcsb_entry_container_identifiers",{}).get("polymer_entity_ids",[]):
        obj=jget(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pid}/{ent}")
        pp=obj.get("entity_poly",{})
        if "polypeptide" not in str(pp.get("type","")).lower():continue
        seq=re.sub(r"[^A-Z]","",pp.get("pdbx_seq_one_letter_code_can","").upper())
        auth=obj.get("rcsb_polymer_entity_container_identifiers",{}).get("auth_asym_ids") or []
        description=obj.get("rcsb_polymer_entity",{}).get("pdbx_description","").lower()
        if any(x in description for x in blockers):raise ValueError("membrane-related entity")
        entities.append({"entity_id":ent,"sequence":seq,"length":len(seq),"auth_chains":auth,"description":description[:160]})
    if not entities:raise ValueError("no protein entity")
    if any(x["length"]<GATES["min_protein_length"] or x["length"]>GATES["max_protein_length"] for x in entities):
        raise ValueError("entity length gate")
    return float(resolution[0]),entities,title[:140]
def coordinate_info(pid,entities):
    raw=fetch(f"https://files.rcsb.org/download/{pid}-assembly1.cif",tries=3,timeout=35)
    cifpath=OUT/(pid+"-assembly1.cif")
    cifpath.write_bytes(raw)
    st=gemmi.read_structure(str(cifpath))
    chains=[]
    for ch in st[0]:
        residues=list(ch.get_polymer())
        aa=[r for r in residues if r.name in PROTEIN]
        if not aa:continue
        if len(residues)!=len(aa):raise ValueError("noncanonical polymer residue")
        original_name=ch.name.split("-")[0]
        refs=[e for e in entities if original_name in e["auth_chains"]]
        if len(refs)!=1:raise ValueError("unmapped or ambiguous chain "+ch.name)
        declared=refs[0]["length"]
        missing=sum(not any(a.name=="CA" for a in r) for r in aa)
        observed=len(aa)-missing
        cov=observed/declared
        if cov<GATES["min_CA_coverage"]:raise ValueError(f"CA coverage {cov:.3f}")
        scmiss=0;alt=0
        for res in aa:
            atoms={a.name for a in res if a.element.name not in ("H","D")}
            scmiss+=int(bool(set(SIDE[res.name].split())-atoms))
            alt+=int(any(a.altloc not in (chr(0)," ",".","?") for a in res))
        if scmiss>GATES["max_missing_sidechain_residues"]:raise ValueError("missing sidechain")
        if alt>GATES["max_altloc_residues"]:raise ValueError("altloc exists")
        chains.append({"chain":ch.name,"length":declared,"n_CA":observed,"coverage":round(cov,4),"sequence":refs[0]["sequence"],"sidechain_missing":scmiss,"altloc":alt})
    if not chains:raise ValueError("no protein chain in bioassembly")
    if len(chains)>GATES["max_assembly_chains"]:raise ValueError("large assembly")
    if sum(x["n_CA"] for x in chains)>GATES["max_assembly_residues"]:raise ValueError("too many residues")
    return {"assembly_sha256":hashlib.sha256(raw).hexdigest(),"assembly_chains":chains,"n_chains":len(chains),"n_residues":sum(x["n_CA"] for x in chains)}
def inspect(pid,cath):
    result={"pdb":pid,"CATH_H":sorted(cath),"decision":"REJECT","reason":None}
    try:
        metrics,rfree,url=wwpdb(pid)
        res,entities,title=entry_data(pid)
        z=coordinate_info(pid,entities)
        result.update({"decision":"QUALIFIED","wwpdb":metrics,"rfree":rfree,"resolution_A":res,"entities":entities,"title":title,"validation_url":url,**z})
    except Exception as e:
        result["reason"]=str(e)[:180]
    return result

def homolog(a,b,aligner):
    for e in a["entities"]:
        for f in b["entities"]:
            seq1=e["sequence"];seq2=f["sequence"]
            if not seq1 or not seq2:continue
            al=aligner.align(seq1,seq2)[0]
            pairs=[(int(i),int(j)) for i,j in zip(*al.indices) if i>=0 and j>=0]
            if len(pairs)<40:continue
            identity=sum(seq1[i]==seq2[j] for i,j in pairs)/len(pairs)
            cov=len(pairs)/min(len(seq1),len(seq2))
            if identity>=0.30 and cov>=0.70:return True
    return False
def pick_candidates(mapping,existing_ids):
    # diversify homologous superfamilies before quality screening.
    group=defaultdict(list)
    for pid,families in mapping.items():
        if pid in existing_ids:continue
        if len(families)>2 or len(families)==0:continue
        if pid[0] not in "3456789":continue
        key=sorted(families)[0]
        group[key].append(pid)
    families=sorted(group,key=lambda x:hashlib.sha256((str(SEED)+x).encode()).hexdigest())
    candidates=[]
    for roundno in range(6):
        for fam in families:
            entries=sorted(group[fam],key=lambda x:(x[0] not in "6789",hashlib.sha256((str(SEED)+x).encode()).hexdigest()))
            if roundno<len(entries):
                pid=entries[roundno]
                if pid not in existing_ids:candidates.append(pid)
    return list(dict.fromkeys(candidates))
def main():
    cath,source,cathhash=cath_map()
    oldpath=ROOT/"validation"/"SCIENTIFIC_QA_48_V4_2026-11-22.csv"
    old=list(csv.DictReader(open(oldpath,encoding="utf-8")))
    oldids=[r["PDB"] for r in old]
    # Pre-screen historical structures by the exactly same strict criteria. 
    # We do not grandfather lower-quality pilot structures.
    candidates=list(dict.fromkeys(oldids+pick_candidates(cath,set(oldids))))[:MAX_CANDIDATES]
    # all candidates must possess domain annotation, otherwise not independent on CATH basis
    candidates=[pid for pid in candidates if pid in cath]
    print("DISCOVERY "+json.dumps({"candidate_count":len(candidates),"cath_entries":len(cath),"cath_source":source,"cath_sha256":cathhash,"gates":GATES}),flush=True)
    selected=[]; qualified=[]; audited=[]; reasons=Counter();used_families=set()
    aligner=PairwiseAligner(mode="local",match_score=2,mismatch_score=-1,open_gap_score=-8,extend_gap_score=-1)
    for k in range(0,len(candidates),BATCH):
        chunk=candidates[k:k+BATCH]
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
            checks=list(ex.map(lambda pid:inspect(pid,cath[pid]),chunk))
        for r in checks:
            audited.append({"pdb":r["pdb"],"decision":r["decision"],"reason":r.get("reason"),"cath":r["CATH_H"]})
            if r["decision"]!="QUALIFIED":
                reasons[r["reason"].split(":")[0] if r.get("reason") else "unknown"]+=1
                continue
            qualified.append(r)
            if set(r["CATH_H"])&used_families:continue
            if any(homolog(r,accepted,aligner) for accepted in selected):continue
            if len(selected)<TARGET:
                selected.append(r);used_families.update(r["CATH_H"])
                print("ACCEPT "+json.dumps({"index":len(selected),"pdb":r["pdb"],"resolution":r["resolution_A"],"rfree":r["rfree"],
                    "CATH":r["CATH_H"],"chains":r["n_chains"],"n_residues":r["n_residues"]}),flush=True)
            if len(selected)>=TARGET:
                break
        print("BATCH_SUMMARY "+json.dumps({"checked":min(k+BATCH,len(candidates)),"pass_qa":len(qualified),"selected_nonhomologous":len(selected),"top_reasons":reasons.most_common(7)}),flush=True)
        if len(selected)>=TARGET:break
    final={"n_selected":len(selected),"required":TARGET,"n_tested":len(audited),
       "n_qualified_before_homology":len(qualified),"pass_48":len(selected)==TARGET,
       "cath_sha256":cathhash,"gates":GATES,"failure_reasons":reasons.most_common(35),
       "original_ids_retained":[x["pdb"] for x in selected if x["pdb"] in oldids],
       "selection_without_interfacescout_scores":True}
    (OUT/"screened_candidates.json").write_text(json.dumps(audited,indent=2))
    (OUT/"quality_qualified.json").write_text(json.dumps(qualified,indent=2))
    (OUT/"selected_48_structures.json").write_text(json.dumps(selected,indent=2))
    (OUT/"selection_summary.json").write_text(json.dumps(final,indent=2))
    with open(OUT/"selected_48.csv","w",newline="") as f:
        w=csv.writer(f);w.writerow(["pdb_id","resolution_A","wwpdb_Rfree","clashscore","rama_pct","rotamer_pct","rsrz_pct","CATH_H","assembly_chains","total_residues","minimum_CA_coverage","assembly_sha256"])
        for x in selected:
            m=x["wwpdb"];w.writerow([x["pdb"],x["resolution_A"],x["rfree"],m["clashscore"],m["percent-rama-outliers"],m["percent-rota-outliers"],m["percent-RSRZ-outliers"],";".join(x["CATH_H"]),x["n_chains"],x["n_residues"],min(y["coverage"] for y in x["assembly_chains"]),x["assembly_sha256"]])
    print("STRICT_SUMMARY "+json.dumps(final,separators=(",",":")),flush=True)
    # Insufficient quality-qualified options => fail closed; no pretend 48.
    if len(selected)!=TARGET:sys.exit(1)
if __name__=="__main__":main()
