from __future__ import annotations
import io, json, re, sys, urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

ROOT=Path(__file__).resolve().parents[2]
HERE=ROOT/"validation"/"benchmark_v1_2026-10-06"
OUT=HERE/"direct_adsorption_bad_gsc_primary_2026-10-07"
OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/"backend"))
import core
from prepare import prepare_pdb_text

UA="InterfaceScout-BAD-GSC/1.0"
MAP_TO_INTERNAL={
"anionic":"anionic","cationic":"cationic","hydrophobic":"hydrophobic","pi/aromatic":"pi_carbon",
"H-bond donor":"hbond_donor","H-bond acceptor":"hbond_acceptor","metal-oxide":"oxide",
"calcium/phosphate charged sites":"hydroxyapatite","transition-metal coordination":"metal_coord",
"soft-metal sulfur affinity":"gold","phosphate-rich":"phosphate"}

def ns(s): return re.sub(r"\s+"," ",str(s or "").strip().lower())
def map_surface(surface):
    s=ns(surface)
    if re.search(r"(carboxyl|\bcooh\b)",s): return ["anionic","H-bond acceptor"],"carboxyl"
    if re.search(r"(amine|amino|\bnh2\b)",s): return ["cationic","H-bond donor"],"amine"
    if re.search(r"(hydroxy|hydroxyl|\boh[- ]?terminated|hydroxy-thiol)",s): return ["H-bond donor","H-bond acceptor"],"hydroxyl"
    if re.search(r"(phosphonate|phosphate[- ]?terminated)",s): return ["phosphate-rich","anionic","H-bond acceptor"],"phosphate"
    if re.search(r"(polyethylene glycol|\bpeg\b)",s): return ["H-bond acceptor"],"PEG"
    if re.search(r"(methyl[- ]?thiol|methyl[- ]?terminated|\bch3\b)",s): return ["hydrophobic"],"methyl"
    if s in {"gold","au","bare gold","unmodified gold"}: return ["soft-metal sulfur affinity"],"bare gold"
    if re.search(r"(silica|silicon oxide|quartz|\bglass\b)",s): return ["metal-oxide"],"silica"
    if re.search(r"(titanium oxide|titania|tio2)",s): return ["metal-oxide"],"titania"
    if re.search(r"(iron oxide|magnetite|fe3o4)",s): return ["metal-oxide"],"iron oxide"
    if re.search(r"(alumina|aluminum oxide|aluminium oxide)",s): return ["metal-oxide"],"alumina"
    if re.search(r"(hydroxyapatite|calcium phosphate)",s): return ["calcium/phosphate charged sites"],"CaP"
    if re.search(r"(poly\(styrene\)|polystyrene|\bps\b)",s): return ["hydrophobic","pi/aromatic"],"polystyrene"
    if re.search(r"(poly\(tetrafluoroethylene\)|tetrafluorethylene|ptfe)",s): return ["hydrophobic"],"PTFE"
    if re.search(r"(graphite|graphene|carbon nanotube|\bcnt\b)",s): return ["hydrophobic","pi/aromatic"],"graphitic carbon"
    if re.search(r"(diamond-like carbon|\bdlc\b)",s): return ["hydrophobic","pi/aromatic"],"DLC"
    if re.search(r"\bmica\b",s): return ["anionic","H-bond acceptor"],"mica"
    return [],"unmapped"

def fnum(x):
    try:
        if pd.isna(x): return None
        s=str(x).strip()
        return None if not s or s.lower() in {"nan","none","na","n/a"} else float(s)
    except: return None

def fetch_bad():
    url="https://molecularsense.com/bad-2/"
    req=urllib.request.Request(url,headers={"User-Agent":UA})
    with urllib.request.urlopen(req,timeout=120) as r: html=r.read().decode("utf-8",errors="replace")
    (OUT/"bad2_snapshot.html").write_text(html,encoding="utf-8")
    for df in pd.read_html(io.StringIO(html)):
        df.columns=[str(c).strip() for c in df.columns]
        if {"Surface Concentration","Solution Concentration","Adsorbing Surface","PDB"}.issubset(df.columns):
            df.to_csv(OUT/"bad2_snapshot.csv",index=False)
            return df,url
    raise RuntimeError("BAD table not found")

def compute_one(key):
    pdb,ph,channels=key
    req=urllib.request.Request(f"https://files.rcsb.org/download/{pdb}.pdb",headers={"User-Agent":UA})
    with urllib.request.urlopen(req,timeout=60) as r: raw=r.read().decode("utf-8",errors="replace")
    prepared,prep=prepare_pdb_text(raw,chain=None)
    res=core.analyze(core.AnalyzeRequest(pdb_text=prepared,chain=None,env=core.EnvParams(pH=ph,ionic=150.0,temp=298.0)))
    if hasattr(res,"body"): raise RuntimeError("core HTTP response")
    surface=res.get("surface_residues",[])
    denom=sum(float(r.get("scrsa") or 0) for r in surface)
    mx={str(r.get("key")):0.0 for r in surface}
    for pub in channels:
        ch=res.get("chemistries",{}).get(MAP_TO_INTERNAL[pub],{})
        for rr in ch.get("residues",[]):
            k=str(rr.get("key"))
            if k in mx: mx[k]=max(mx[k],float(rr.get("local_score") or 0))
    return {"GSC":sum(mx.values())/denom if denom>0 else np.nan,
            "n_surface_residues":len(surface),"surface_scrsa_sum":denom,
            "selected_chain":prep.get("selected_chain","ALL")}

def main():
    df,url=fetch_bad()
    prelim=[]; excl=[]
    needed=set()
    for _,r in df.iterrows():
        rid=str(r.get("ID","")).strip(); pdb=str(r.get("PDB","")).strip().upper()
        protein=str(r.get("Protein","")).strip(); surf=str(r.get("Adsorbing Surface","")).strip()
        gamma=fnum(r.get("Surface Concentration")); csol=fnum(r.get("Solution Concentration"))
        ph=fnum(r.get("pH")); ionic=fnum(r.get("Ionic Strength")); temp=fnum(r.get("Temperature"))
        method=str(r.get("Measurement Method","")).strip(); etype=str(r.get("Experiment Type","")).strip()
        val=str(r.get("Data Validation",r.get("Data validation",""))).strip(); notes=str(r.get("Explanatory Notes","")).strip()
        chans,rule=map_surface(surf)
        reasons=[]
        if not re.fullmatch(r"[0-9A-Za-z]{4}",pdb): reasons.append("invalid_or_missing_pdb")
        if gamma is None or gamma<=0: reasons.append("bad_surface_concentration")
        if csol is None or csol<=0: reasons.append("bad_solution_concentration")
        if ph is None: reasons.append("missing_pH")
        if not chans: reasons.append("surface_unmapped")
        if re.search(r"(uncertain|invalid)",val+" "+notes,re.I): reasons.append("source_flag_uncertain_or_invalid")
        base={"source_row_id":rid,"protein":protein,"pdb_id":pdb,"surface_concentration_mg_m2":gamma,
              "solution_concentration_mg_ml":csol,"adsorbing_surface":surf,"pH":ph,"ionic_strength":ionic,
              "temperature_C":temp,"measurement_method":method,"experiment_type":etype,
              "reference":str(r.get("Reference","")).strip(),"DOI":str(r.get("DOI","")).strip(),
              "mapped_channels":";".join(chans),"surface_mapping_rule":rule}
        if reasons: excl.append({**base,"exclusion_reasons":";".join(reasons)}); continue
        key=(pdb,float(ph),tuple(chans)); needed.add(key); prelim.append((base,key))
    scores={}
    with ThreadPoolExecutor(max_workers=8) as ex:
        fut={ex.submit(compute_one,k):k for k in needed}
        for f in as_completed(fut):
            k=fut[f]
            try:scores[k]=f.result()
            except Exception as e:scores[k]={"error":repr(e)}
    rows=[]
    for base,key in prelim:
        sc=scores[key]
        if "error" in sc:
            excl.append({**base,"exclusion_reasons":"v1_or_structure_failure","score_error":sc["error"]})
        else: rows.append({**base,**sc})
    pd.DataFrame(rows).to_csv(OUT/"bad2_gsc_scored_rows.csv",index=False)
    pd.DataFrame(excl).to_csv(OUT/"bad2_gsc_exclusions.csv",index=False)

    def ex(v):
        if v is None:return ""
        return f"{v:.12g}" if isinstance(v,float) else str(v).strip()
    groups=defaultdict(list)
    for r in rows:
        k=(ns(r["adsorbing_surface"]),ex(r["solution_concentration_mg_ml"]),ex(r["pH"]),
           ex(r["ionic_strength"]),ex(r["temperature_C"]),r["measurement_method"].lower(),r["experiment_type"].lower())
        groups[k].append(r)
    strata=[]
    for k,arr in groups.items():
        if len(set(x["pdb_id"] for x in arr))<3: continue
        x=[float(z["GSC"]) for z in arr]; y=[float(z["surface_concentration_mg_m2"]) for z in arr]
        if len(set(x))<2 or len(set(y))<2: continue
        rho,p=spearmanr(x,y)
        if not np.isfinite(rho):continue
        strata.append({"surface":k[0],"solution_concentration_mg_ml":k[1],"pH":k[2],
                       "ionic_strength":k[3],"temperature_C":k[4],"measurement_method":k[5],
                       "experiment_type":k[6],"n_rows":len(arr),"n_unique_pdb":len(set(z["pdb_id"] for z in arr)),
                       "spearman_rho":float(rho),"spearman_p":float(p) if np.isfinite(p) else "",
                       "protein_pdbs":";".join(sorted(set(z["pdb_id"] for z in arr)))})
    pd.DataFrame(strata).to_csv(OUT/"bad2_gsc_primary_strata.csv",index=False)
    rhos=np.array([z["spearman_rho"] for z in strata],float)
    summ={"n_strata":len(strata),"median_rho":float(np.median(rhos)) if len(rhos) else None,
          "q1_rho":float(np.quantile(rhos,.25)) if len(rhos) else None,
          "q3_rho":float(np.quantile(rhos,.75)) if len(rhos) else None,
          "fraction_positive":float(np.mean(rhos>0)) if len(rhos) else None,
          "wilcoxon_p_vs_zero":None}
    if len(rhos) and np.any(np.abs(rhos)>1e-15):
        try:summ["wilcoxon_p_vs_zero"]=float(wilcoxon(rhos,zero_method="wilcox").pvalue)
        except:pass
    audit={"source_url":url,"live_rows":len(df),"eligible_scored_rows":len(rows),"excluded_rows":len(excl),
           "unique_scored_pdb":len(set(r["pdb_id"] for r in rows)),"unique_scored_surfaces":len(set(r["adsorbing_surface"] for r in rows)),
           "unique_gsc_calculations_success":sum("error" not in v for v in scores.values()),
           "unique_gsc_calculations_failed":sum("error" in v for v in scores.values()),
           "exact_comparable_groups_with_3plus_pdb":sum(len(set(x["pdb_id"] for x in a))>=3 for a in groups.values()),
           "primary_evaluable_strata":len(strata),"experimental_adsorption_values_used_in_GSC":False,"v1_refit":False}
    (OUT/"audit.json").write_text(json.dumps(audit,indent=2))
    (OUT/"summary.json").write_text(json.dumps(summ,indent=2))
    print(json.dumps(audit,indent=2));print(json.dumps(summ,indent=2))
if __name__=="__main__":main()
