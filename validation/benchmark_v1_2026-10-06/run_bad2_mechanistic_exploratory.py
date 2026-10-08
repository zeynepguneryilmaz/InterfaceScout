from __future__ import annotations
import io, json, math, re, sys, urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon
from Bio.PDB import PDBParser
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import seq1

ROOT=Path(__file__).resolve().parents[2]
HERE=ROOT/"validation"/"benchmark_v1_2026-10-06"
OUT=HERE/"direct_adsorption_mechanistic_2026-10-08"
OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/"backend"))

import core
from prepare import prepare_pdb_text
from interface_engine import analyze_all_maps

UA="InterfaceScout-BAD-mechanistic/1.0"
MAP_TO_INTERNAL={
    "anionic":"anionic","cationic":"cationic","hydrophobic":"hydrophobic","pi/aromatic":"pi_carbon",
    "H-bond donor":"hbond_donor","H-bond acceptor":"hbond_acceptor","metal-oxide":"oxide",
    "calcium/phosphate charged sites":"hydroxyapatite","transition-metal coordination":"metal_coord",
    "soft-metal sulfur affinity":"gold","phosphate-rich":"phosphate",
}
PUBLIC_KEY={
    "anionic":"anionic","cationic":"cationic","hydrophobic":"hydrophobic","pi/aromatic":"pi_carbon",
    "H-bond donor":"hbond_donor","H-bond acceptor":"hbond_acceptor","metal-oxide":"oxide",
    "calcium/phosphate charged sites":"calcium_phosphate","transition-metal coordination":"metal_coord",
    "soft-metal sulfur affinity":"soft_metal_sulfur","phosphate-rich":"phosphate",
}

def ns(s): return re.sub(r"\s+"," ",str(s or "").strip().lower())

def map_surface(surface):
    s=ns(surface)
    if re.search(r"(carboxyl|\bcooh\b)",s): return ["anionic","H-bond acceptor"],"carboxyl"
    if re.search(r"(amine|amino|\bnh2\b)",s): return ["cationic","H-bond donor"],"amine"
    if re.search(r"(hydroxy|hydroxyl|\boh[- ]?terminated|hydroxy-thiol)",s): return ["H-bond donor","H-bond acceptor"],"hydroxyl"
    if re.search(r"(phosphonate|phosphate[- ]?terminated)",s): return ["phosphate-rich","anionic","H-bond acceptor"],"phosphate"
    if re.search(r"(polyethylene glycol|\bpeg\b)",s): return ["H-bond acceptor"],"PEG"
    if re.search(r"(methylated|methyl[- ]?thiol|methyl[- ]?terminated|\bch3\b)",s): return ["hydrophobic"],"methyl"
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
    with urllib.request.urlopen(req,timeout=120) as r:
        html=r.read().decode("utf-8",errors="replace")
    (OUT/"bad2_snapshot.html").write_text(html,encoding="utf-8")
    for df in pd.read_html(io.StringIO(html)):
        df.columns=[str(c).strip() for c in df.columns]
        if {"Surface Concentration","Solution Concentration","Adsorbing Surface","PDB"}.issubset(df.columns):
            df.to_csv(OUT/"bad2_snapshot.csv",index=False)
            return df,url
    raise RuntimeError("BAD table not found")

def structure_mw_from_prepared(prepared):
    parser=PDBParser(QUIET=True)
    st=parser.get_structure("x",io.StringIO(prepared))
    aas=[]
    for res in st.get_residues():
        if res.id[0]!=" ":
            continue
        rn=res.get_resname().strip().upper()
        try:
            aa=seq1(rn)
        except:
            continue
        if aa and len(aa)==1 and aa in "ACDEFGHIKLMNPQRSTVWY":
            aas.append(aa)
    if not aas:
        return None
    try:
        return float(ProteinAnalysis("".join(aas)).molecular_weight())
    except:
        return None

def compute_descriptors(key):
    pdb,ph,channels=key
    req=urllib.request.Request(f"https://files.rcsb.org/download/{pdb}.pdb",headers={"User-Agent":UA})
    with urllib.request.urlopen(req,timeout=60) as r:
        raw=r.read().decode("utf-8",errors="replace")
    prepared,prep=prepare_pdb_text(raw,chain=None)
    cres=core.analyze(core.AnalyzeRequest(pdb_text=prepared,chain=None,env=core.EnvParams(pH=ph,ionic=150.0,temp=298.0)))
    if hasattr(cres,"body"): raise RuntimeError("core HTTP response")
    surface=cres.get("surface_residues",[])
    denom=sum(float(r.get("scrsa") or 0.0) for r in surface)
    scrsa={str(r["key"]):float(r.get("scrsa") or 0.0) for r in surface}
    fav={k:0.0 for k in scrsa}
    rep={k:0.0 for k in scrsa}
    for pub in channels:
        ch=cres.get("chemistries",{}).get(MAP_TO_INTERNAL[pub],{})
        for rr in ch.get("residues",[]):
            k=str(rr.get("key"))
            if k in fav: fav[k]=max(fav[k],float(rr.get("local_score") or 0.0))
        for rr in ch.get("repulsive_residues",[]):
            k=str(rr.get("key"))
            if k in rep:
                # repulsion_propensity is normalized; reconstruct physical local score using residue chemistry definitions
                rn=str(rr.get("res_name") or "").upper()
                defs=core.CHEMISTRIES[MAP_TO_INTERNAL[pub]].get("repulsive",{})
                if rn in defs:
                    state=defs[rn][2]
                    rep[k]=max(rep[k],scrsa[k]*core.state_availability(rn,state,ph))

    ngc=(sum(fav.values())-sum(rep.values()))/denom if denom>0 else np.nan
    gfc=sum(fav.values())/denom if denom>0 else np.nan
    grc=sum(rep.values())/denom if denom>0 else np.nan

    amap=analyze_all_maps(pdb_text=raw,chain=None,pH=float(ph),ionic_mM=150.0,temp_K=298.0)
    patches=[]
    for pub in channels:
        pk=PUBLIC_KEY[pub]
        for p in amap["maps"].get(pk,{}).get("patches",[]):
            if int(p.get("pareto_front",999))!=1: continue
            members=[str(x) for x in p.get("members",[])]
            pd=sum(scrsa.get(k,0.0) for k in members)
            if pd<=0: continue
            pf=sum(fav.get(k,0.0) for k in members)
            pr=sum(rep.get(k,0.0) for k in members)
            net=(pf-pr)/pd
            favcomp=pf/pd
            cov=pd/denom if denom>0 else np.nan
            patches.append((net,favcomp,cov,pk,str(p.get("center_key",""))))
    if patches:
        patches.sort(key=lambda x:(-x[0],-x[1],-x[2],x[3],x[4]))
        bpnc=patches[0][0]
        bpsc=patches[0][2]
        bpfc=max(x[1] for x in patches)
        lns=max(max(0.0,x[0])*x[2] for x in patches)
        n_front1=len(patches)
    else:
        bpnc=bpfc=bpsc=lns=0.0
        n_front1=0

    mw=structure_mw_from_prepared(prepared)
    return {
        "GFC":gfc,
        "GRC":grc,
        "NGC":ngc,
        "BPNC":bpnc,
        "BPFC":bpfc,
        "BPSC":bpsc,
        "LNS":lns,
        "n_front1_patches":n_front1,
        "structure_mw_Da":mw,
        "n_surface_residues":len(surface),
        "surface_scrsa_sum":denom,
        "selected_chain":prep.get("selected_chain","ALL"),
    }

def ex(v):
    if v is None:return ""
    return f"{v:.12g}" if isinstance(v,float) else str(v).strip()

def main():
    df,url=fetch_bad()
    prelim=[];excl=[];needed=set()
    for _,r in df.iterrows():
        rid=str(r.get("ID","")).strip();pdb=str(r.get("PDB","")).strip().upper()
        protein=str(r.get("Protein","")).strip();surf=str(r.get("Adsorbing Surface","")).strip()
        gamma=fnum(r.get("Surface Concentration"));csol=fnum(r.get("Solution Concentration"))
        ph=fnum(r.get("pH"));ionic=fnum(r.get("Ionic Strength"));temp=fnum(r.get("Temperature"))
        method=str(r.get("Measurement Method","")).strip();etype=str(r.get("Experiment Type","")).strip()
        val=str(r.get("Data Validation",r.get("Data validation",""))).strip();notes=str(r.get("Explanatory Notes","")).strip()
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
        if reasons:
            excl.append({**base,"exclusion_reasons":";".join(reasons)});continue
        key=(pdb,float(ph),tuple(chans));needed.add(key);prelim.append((base,key))

    scores={}
    with ThreadPoolExecutor(max_workers=8) as exr:
        fut={exr.submit(compute_descriptors,k):k for k in needed}
        for f in as_completed(fut):
            k=fut[f]
            try:scores[k]=f.result()
            except Exception as e:scores[k]={"error":repr(e)}

    rows=[]
    for base,key in prelim:
        sc=scores[key]
        if "error" in sc:
            excl.append({**base,"exclusion_reasons":"v1_or_structure_failure","score_error":sc["error"]})
            continue
        z={**base,**sc}
        mw=sc.get("structure_mw_Da")
        z["molar_surface_density_umol_m2"]=(
            float(base["surface_concentration_mg_m2"])*1e3/float(mw)
            if mw and mw>0 else np.nan
        )
        rows.append(z)

    pd.DataFrame(rows).to_csv(OUT/"scored_rows.csv",index=False)
    pd.DataFrame(excl).to_csv(OUT/"exclusions.csv",index=False)

    # exact strata; collapse repeated same-PDB rows to median response
    groups=defaultdict(list)
    for r in rows:
        k=(ns(r["adsorbing_surface"]),ex(r["solution_concentration_mg_ml"]),ex(r["pH"]),
           ex(r["ionic_strength"]),ex(r["temperature_C"]),r["measurement_method"].lower(),r["experiment_type"].lower())
        groups[k].append(r)

    descriptors=["GFC","GRC","NGC","BPNC","BPFC","BPSC","LNS"]
    targets=["surface_concentration_mg_m2","molar_surface_density_umol_m2"]
    results=[]
    collapsed=[]
    for k,arr in groups.items():
        bypdb=defaultdict(list)
        for x in arr:bypdb[x["pdb_id"]].append(x)
        if len(bypdb)<3:continue
        reps=[]
        for pdb,aa in bypdb.items():
            first=dict(aa[0])
            for target in targets:
                vals=[float(x[target]) for x in aa if x.get(target) is not None and np.isfinite(float(x[target]))]
                first[target]=float(np.median(vals)) if vals else np.nan
            first["n_replicate_rows_collapsed"]=len(aa)
            reps.append(first)
            collapsed.append({**first,"stratum_surface":k[0],"stratum_solution_concentration":k[1]})
        for target in targets:
            for desc in descriptors:
                xs=[];ys=[]
                for x in reps:
                    try: xv=float(x[desc]);yv=float(x[target])
                    except:continue
                    if np.isfinite(xv) and np.isfinite(yv):
                        xs.append(xv);ys.append(yv)
                if len(xs)<3 or len(set(xs))<2 or len(set(ys))<3:continue
                rho,p=spearmanr(xs,ys)
                if not np.isfinite(rho):continue
                results.append({
                    "surface":k[0],"solution_concentration_mg_ml":k[1],"pH":k[2],
                    "ionic_strength":k[3],"temperature_C":k[4],"measurement_method":k[5],
                    "experiment_type":k[6],"target":target,"descriptor":desc,
                    "n_unique_pdb":len(xs),"spearman_rho":float(rho),
                    "spearman_p":float(p) if np.isfinite(p) else "",
                    "protein_pdbs":";".join(sorted(bypdb)),
                })

    pd.DataFrame(collapsed).to_csv(OUT/"collapsed_exact_strata_rows.csv",index=False)
    pd.DataFrame(results).to_csv(OUT/"stratum_correlations.csv",index=False)

    summary=[]
    for target in targets:
        for desc in descriptors:
            vals=np.array([r["spearman_rho"] for r in results if r["target"]==target and r["descriptor"]==desc],float)
            if len(vals)==0:continue
            p=""
            if np.any(np.abs(vals)>1e-15):
                try:p=float(wilcoxon(vals,zero_method="wilcox",alternative="two-sided").pvalue)
                except:pass
            summary.append({
                "target":target,"descriptor":desc,"n_evaluable_strata":len(vals),
                "median_rho":float(np.median(vals)),"q1_rho":float(np.quantile(vals,.25)),
                "q3_rho":float(np.quantile(vals,.75)),"fraction_positive":float(np.mean(vals>0)),
                "wilcoxon_p_vs_zero":p
            })
    pd.DataFrame(summary).to_csv(OUT/"descriptor_summary.csv",index=False)

    audit={
        "source_url":url,"live_rows":len(df),"eligible_scored_rows":len(rows),"excluded_rows":len(excl),
        "unique_scored_pdb":len(set(r["pdb_id"] for r in rows)),"unique_scored_surfaces":len(set(r["adsorbing_surface"] for r in rows)),
        "unique_descriptor_calculations_success":sum("error" not in v for v in scores.values()),
        "unique_descriptor_calculations_failed":sum("error" in v for v in scores.values()),
        "exact_groups_with_3plus_distinct_pdb":sum(len(set(x["pdb_id"] for x in a))>=3 for a in groups.values()),
        "correlation_rows":len(results),"post_gsc_failure_exploratory":True,
        "coefficients_fitted":False,"BAD_values_used_in_descriptor_definition":False,
        "same_pdb_replicates_collapsed_to_median":True,
        "monolayer_claim_made":False
    }
    (OUT/"audit.json").write_text(json.dumps(audit,indent=2))
    print(json.dumps(audit,indent=2))
    print(json.dumps(summary,indent=2))

if __name__=="__main__":main()
