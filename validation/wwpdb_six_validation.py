"""Audit independent wwPDB validation summary metrics for six selected candidates."""
import gzip,json,requests,xml.etree.ElementTree as ET,sys
from pathlib import Path
IDS=["1C5E","1PDB","1A1X","1STN","2CBA","1RWZ"]
out=Path("wwpdb_six_validation");out.mkdir(exist_ok=True)
items=[]
for pid in IDS:
    report={"pdb":pid}
    url=f"https://files.rcsb.org/pub/pdb/validation_reports/{pid[1:3].lower()}/{pid.lower()}/{pid.lower()}_validation.xml.gz"
    try:
        r=requests.get(url,timeout=45)
        r.raise_for_status()
        root=ET.fromstring(gzip.decompress(r.content))
        en=next((x for x in root.iter() if x.tag.rsplit("}",1)[-1]=="Entry"),None)
        if en is None:raise ValueError("No Entry node in validation XML")
        metrics=en.attrib
        wanted=("DCC_Rfree","clashscore","percent-rama-outliers","percent-rota-outliers","percent-RSRZ-outliers","bonds_rmsz","angles_rmsz","high-resol","low-resol")
        report["metrics"]={k:metrics.get(k) for k in wanted}
        report["raw_available_keys"]=sorted(metrics.keys())
        residues=[x for x in root.iter() if x.tag.rsplit("}",1)[-1]=="ModelledSubgroup" and x.attrib.get("resname","") in ("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL")]
        rama=[x for x in residues if x.attrib.get("rama","").lower().startswith("out")]
        rota=[x for x in residues if x.attrib.get("rota","").lower().startswith("out")]
        report["modelled_protein_rows"]=len(residues)
        report["observed_rama_outlier_rows"]=len(rama)
        report["observed_rota_outlier_rows"]=len(rota)
        report["rama_outlier_examples"]=[{"chain":x.attrib.get("chain"),"resnum":x.attrib.get("resnum"),"name":x.attrib.get("resname")} for x in rama[:10]]
        report["rota_outlier_examples"]=[{"chain":x.attrib.get("chain"),"resnum":x.attrib.get("resnum"),"name":x.attrib.get("resname")} for x in rota[:10]]
        report["report_status"]="FETCHED"
        report["source"]=url
    except Exception as e:
        report["report_status"]="UNAVAILABLE"
        report["reason"]=str(e)
    items.append(report)
    print("WWPDB_QA",json.dumps(report,separators=(",",":")),flush=True)
(out/"metrics.json").write_text(json.dumps(items,indent=2))
if any(x.get("report_status")!="FETCHED" for x in items):raise SystemExit("One or more report documents unavailable")
