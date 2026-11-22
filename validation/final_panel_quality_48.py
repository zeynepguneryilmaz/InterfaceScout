"""Auditable wwPDB+CATH+coordinate screening of the frozen 48-ID development panel.
Scientific gate: missing data stay missing. This does NOT certify chemical validity.
"""
import csv, gzip, hashlib, io, json, re, sys, time
from collections import defaultdict
from pathlib import Path
import requests, gemmi

ROOT=Path(__file__).resolve().parents[1]
OUT=Path("final_panel_qa");OUT.mkdir(exist_ok=True)
manifest=list(csv.DictReader(open(ROOT/"validation"/"development_panel_48_candidates_v4_REVIEW.csv", encoding="utf-8")))
IDS=[r["pdb_id"].upper() for r in manifest]
assert len(IDS)==48 and len(set(IDS))==48
session=requests.Session()
session.headers["User-Agent"]="InterfaceScout 48 scientific quality audit / publication validation"
CANON={"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
SC={"ALA":"CB","ARG":"CB CG CD NE CZ NH1 NH2","ASN":"CB CG OD1 ND2","ASP":"CB CG OD1 OD2","CYS":"CB SG","GLN":"CB CG CD OE1 NE2","GLU":"CB CG CD OE1 OE2","GLY":"","HIS":"CB CG ND1 CD2 CE1 NE2","ILE":"CB CG1 CG2 CD1","LEU":"CB CG CD1 CD2","LYS":"CB CG CD CE NZ","MET":"CB CG SD CE","PHE":"CB CG CD1 CD2 CE1 CE2 CZ","PRO":"CB CG CD","SER":"CB OG","THR":"CB OG1 CG2","TRP":"CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2","TYR":"CB CG CD1 CD2 CE1 CE2 CZ OH","VAL":"CB CG1 CG2"}

def download(url, tries=3):
    err=None
    for attempt in range(tries):
        try:
            r=session.get(url,timeout=40)
            r.raise_for_status()
            return r.content
        except Exception as e:
            err=e;time.sleep(attempt+1)
    raise RuntimeError(f"{url}: {err}")

def getjson(url):
    return json.loads(download(url))

def get_validation(pid):
    # current public wwPDB PDB archive validation report
    stem=pid.lower()
    url=f"https://files.rcsb.org/pub/pdb/validation_reports/{stem[1:3]}/{stem}/{stem}_validation.xml.gz"
    try:
        raw=download(url)
        if raw[:2]==b"\x1f\x8b":raw=gzip.decompress(raw)
        import xml.etree.ElementTree as ET
        root=ET.fromstring(raw)
        attrs={}
        entry_metrics={}
        for elem in root.iter():
            tag=elem.tag.split("}")[-1].lower()
            if tag in ("entry","overall","validationreport","experiment"):
                attrs.setdefault(tag,[]).append(dict(elem.attrib))
            if tag=="entry" and not entry_metrics:
                entry_metrics=dict(elem.attrib)
        # extract validated residue outlier counts from modelled subgroup reporting
        counts=defaultdict(int)
        for elem in root.iter():
            if elem.tag.split("}")[-1].lower()=="modelledsubgroup":
                counts["modelled_subgroup"]+=1
                d={k.lower():str(v).lower() for k,v in elem.attrib.items()}
                if d.get("rama")=="outlier":counts["rama_outlier"]+=1
                if d.get("rota")=="outlier":counts["rota_outlier"]+=1
                if d.get("rsrz") not in (None,"", "none"):
                    try:
                        if float(d["rsrz"])>2:counts["rsrz_gt2"]+=1
                    except ValueError:pass
        return {"available":True,"source":url,"sha256":hashlib.sha256(raw).hexdigest(),
                "root":root.tag.split("}")[-1],"summary_attributes":attrs,
                "metrics":{k:entry_metrics.get(k) for k in ("clashscore","percent-rama-outliers","percent-rota-outliers","percent-RSRZ-outliers","DCC_Rfree","PDB-Rfree")},
                "residue_outlier_counts":dict(counts)}
    except Exception as e:
        return {"available":False,"source":url,"error":str(e)}

def get_cath():
    urls=[
      "https://download.cathdb.info/cath/releases/latest-release/cath-classification-data/cath-domain-list.txt",
      "https://download.cathdb.info/cath/releases/all-releases/v4_4_0/cath-classification-data/cath-domain-list-v4_4_0.txt",
    ]
    err=[]
    for url in urls:
        try:
            raw=download(url,2)
            if raw[:2]==b"\x1f\x8b":raw=gzip.decompress(raw)
            mapping=defaultdict(list)
            for line in raw.decode("utf-8","replace").splitlines():
                if line.startswith("#"):continue
                p=line.split()
                if len(p)<5 or len(p[0])<7:continue
                pid=p[0][:4].upper()
                if pid in IDS:
                    mapping[pid].append({"domain_id":p[0],"cath_superfamily":".".join(p[1:5]),"class":".".join(p[1:4])})
            print("CATH_SOURCE "+json.dumps({"url":url,"bytes":len(raw),"matched_pdb":len(mapping),"assignments":sum(map(len,mapping.values()))}),flush=True)
            return dict(mapping),url,[]
        except Exception as e:
            err.append(str(e))
    return {},None,err

def chain_assessment(pid, entities):
    raw=download(f"https://files.rcsb.org/download/{pid}-assembly1.cif")
    path=OUT/f"{pid}-assembly1.cif"
    path.write_bytes(raw)
    st=gemmi.read_structure(str(path))
    chains=[]
    for ch in st[0]:
        p=[r for r in ch.get_polymer() if r.name in CANON]
        if not p:continue
        missing=0;alt=0;no_ca=0
        for r in p:
            atoms={a.name for a in r if a.element.name not in ("H","D")}
            missing+=bool(set(SC[r.name].split())-atoms)
            alt+=any(a.altloc not in (chr(0)," ",".","?") for a in r)
            no_ca+="CA" not in atoms
        # match author asym ID (including symmetry-generated A-2)
        auth=ch.name.split("-")[0]
        possible=[e["seq_len"] for e in entities if auth in e["auth_chains"]]
        decl=max(possible) if possible else None
        coverage=round((len(p)-no_ca)/decl,4) if decl else None
        chains.append({"chain":ch.name,"residues":len(p),"ca_missing":no_ca,
          "seq_declared":decl,"ca_coverage":coverage,"sidechain_missing_residues":missing,
          "alternate_location_residues":alt})
    return {"assembly1_cif_sha256":hashlib.sha256(raw).hexdigest(),"chain_count":len(chains),
            "residue_count":sum(x["residues"] for x in chains),"chains":chains}

def entries(pid):
    e=getjson(f"https://data.rcsb.org/rest/v1/core/entry/{pid}")
    ents=[]
    for eid in e.get("rcsb_entry_container_identifiers",{}).get("polymer_entity_ids",[]):
        p=getjson(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pid}/{eid}")
        poly=p.get("entity_poly",{})
        if "polypeptide" not in str(poly.get("type","")).lower():continue
        seq=re.sub(r"[^A-Z]","",poly.get("pdbx_seq_one_letter_code_can","").upper())
        chains=p.get("rcsb_polymer_entity_container_identifiers",{}).get("auth_asym_ids") or []
        annotations=p.get("rcsb_polymer_entity_annotation") or []
        ents.append({"entity_id":str(eid),"seq_len":len(seq),"sequence":seq,"auth_chains":chains,
          "annotation":[{"type":q.get("type"),"id":q.get("annotation_id"),"name":q.get("name"),"source":q.get("provenance_source")} for q in annotations]})
    ref=e.get("refine") or []
    rf=next((x.get("ls_R_factor_R_free") for x in ref if x.get("ls_R_factor_R_free") is not None),None)
    rr=e.get("rcsb_entry_info",{}).get("resolution_combined") or []
    exp=[x.get("method") for x in e.get("exptl",[])]
    return {"resolution":rr[0] if rr else None,"rfree":rf,"method":exp,"entities":ents}

def main():
    cath,cath_src,cath_errors=get_cath()
    result=[]
    for i,pid in enumerate(IDS,1):
        row={"pdb":pid}
        try:
            row.update(entries(pid))
            row["assembly"]=chain_assessment(pid,row["entities"])
            row["wwpdb"]=get_validation(pid)
            row["cath"]=cath.get(pid,[])
            flags=[]
            if row["method"]!=["X-RAY DIFFRACTION"]:flags.append("NON_XRAY")
            if row["resolution"] is None or row["resolution"]>2.5:flags.append("RESOLUTION_EXCEPTION_OR_FAIL")
            if row["rfree"] is not None and row["rfree"]>0.30:flags.append("RFREE_GT_0.30")
            if row["rfree"] is None:flags.append("RFREE_MISSING")
            if not row["wwpdb"]["available"]:flags.append("WWPDB_REPORT_MISSING")
            else:
                m=row["wwpdb"].get("metrics",{})
                for key,limit,flag in [("clashscore",25,"WWPDB_CLASHSCORE_GT25"),("percent-rama-outliers",2,"WWPDB_RAMA_OUTLIERS_GT2PCT"),("percent-rota-outliers",10,"WWPDB_ROTAMER_OUTLIERS_GT10PCT"),("percent-RSRZ-outliers",10,"WWPDB_RSRZ_OUTLIERS_GT10PCT")]:
                    try:
                        if float(m[key])>limit:flags.append(flag)
                    except (TypeError,ValueError,KeyError):pass
            if not row["cath"]:flags.append("CATH_UNANNOTATED")
            if not row["assembly"]["chains"]:flags.append("NO_PROTEIN_CHAINS")
            for c in row["assembly"]["chains"]:
                if c["ca_coverage"] is None:flags.append("COVERAGE_UNMAPPED_"+c["chain"])
                elif c["ca_coverage"]<.90:flags.append("COVERAGE_LT90_"+c["chain"])
                if c["sidechain_missing_residues"]:flags.append("MISSING_SIDECHAIN_"+c["chain"])
                if c["alternate_location_residues"]:flags.append("ALTLOC_"+c["chain"])
            row["flags"]=flags
        except Exception as ex:
            row["error"]=str(ex);row["flags"]=["FETCH_OR_PARSE_FAILED"]
        result.append(row)
        print(f"QARECORD {i}/48 "+json.dumps({"pdb":pid,"resolution":row.get("resolution"),"rfree":row.get("rfree"),
            "wwpdb":row.get("wwpdb",{}).get("available"),"wwpdb_metrics":row.get("wwpdb",{}).get("metrics",{}),"cath":[c["cath_superfamily"] for c in row.get("cath",[])],
            "chains":row.get("assembly",{}).get("chains",[]),"flags":row["flags"]},separators=(",",":")),flush=True)
    cls=defaultdict(set)
    for r in result:
        for d in r.get("cath",[]):cls[d["cath_superfamily"]].add(r["pdb"])
    clashes={k:sorted(v) for k,v in cls.items() if len(v)>1}
    for r in result:
        r["shared_cath_superfamilies"]=[{"superfamily":k,"other_pdb":sorted(set(v)-{r["pdb"]})} for k,v in clashes.items() if r["pdb"] in v]
        if r["shared_cath_superfamilies"]:r["flags"].append("CATH_HOMOLOGOUS_SUPERFAMILY_SHARED")
    (OUT/"quality_48.json").write_text(json.dumps(result,indent=2))
    (OUT/"quality_summary.json").write_text(json.dumps({
      "n_ids":len(IDS),"n_entries_metadata":sum("error" not in r for r in result),
      "n_wwpdb_xml_fetched":sum(r.get("wwpdb",{}).get("available",False) for r in result),
      "n_cath_annotated":sum(bool(r.get("cath")) for r in result),
      "shared_cath_superfamilies":clashes,"cath_source":cath_src,"cath_errors":cath_errors,
      "not_equivalent_to_final_approval":True},indent=2))
    with (OUT/"quality_48.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["pdb","resolution_A","rfree","assembly_chains","min_CA_coverage","missing_sidechain_residues","wwpdb_report","clashscore","rama_pct","rota_pct","rsrz_pct","cath_superfamilies","shared_cath","flags","error"])
        w.writeheader()
        for r in result:
            cs=r.get("assembly",{}).get("chains",[])
            cov=[c["ca_coverage"] for c in cs if c["ca_coverage"] is not None]
            w.writerow({"pdb":r["pdb"],"resolution_A":r.get("resolution"),"rfree":r.get("rfree"),
              "assembly_chains":len(cs),"min_CA_coverage":min(cov) if cov else "",
              "missing_sidechain_residues":sum(c["sidechain_missing_residues"] for c in cs),
              "wwpdb_report":r.get("wwpdb",{}).get("available",False),
              "clashscore":r.get("wwpdb",{}).get("metrics",{}).get("clashscore"),
              "rama_pct":r.get("wwpdb",{}).get("metrics",{}).get("percent-rama-outliers"),
              "rota_pct":r.get("wwpdb",{}).get("metrics",{}).get("percent-rota-outliers"),
              "rsrz_pct":r.get("wwpdb",{}).get("metrics",{}).get("percent-RSRZ-outliers"),
              "cath_superfamilies":";".join(d["cath_superfamily"] for d in r.get("cath",[])),
              "shared_cath":";".join(d["superfamily"] for d in r.get("shared_cath_superfamilies",[])),
              "flags":";".join(r["flags"]),"error":r.get("error","")})
    summary=json.load(open(OUT/"quality_summary.json"))
    print("QUALITY_SUMMARY "+json.dumps(summary,separators=(",",":")),flush=True)
if __name__=="__main__":main()
