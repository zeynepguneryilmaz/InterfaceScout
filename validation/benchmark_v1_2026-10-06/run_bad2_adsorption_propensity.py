from __future__ import annotations
import csv, io, json, math, re, sys, urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

ROOT=Path(__file__).resolve().parents[2]
HERE=ROOT/"validation"/"benchmark_v1_2026-10-06"
OUT=HERE/"direct_adsorption_bad_v1_2026-10-07"
OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/"backend"))

import core
from prepare import prepare_pdb_text

UA="InterfaceScout-BAD-validation/1.0"

MAP_TO_INTERNAL={
    "anionic":"anionic",
    "cationic":"cationic",
    "hydrophobic":"hydrophobic",
    "pi/aromatic":"pi_carbon",
    "H-bond donor":"hbond_donor",
    "H-bond acceptor":"hbond_acceptor",
    "metal-oxide":"oxide",
    "calcium/phosphate charged sites":"hydroxyapatite",
    "transition-metal coordination":"metal_coord",
    "soft-metal sulfur affinity":"gold",
    "phosphate-rich":"phosphate",
}

def norm_surface(s):
    return re.sub(r"\s+"," ",str(s or "").strip().lower())

def map_surface(surface):
    s=norm_surface(surface)
    # explicit outer functionalization first
    if re.search(r"(carboxyl|\bcooh\b)",s):
        return ["anionic","H-bond acceptor"],"explicit carboxyl outer surface"
    if re.search(r"(amine|amino|\bnh2\b)",s):
        return ["cationic","H-bond donor"],"explicit amine outer surface"
    if re.search(r"(hydroxy|hydroxyl|\boh[- ]?terminated|hydroxy-thiol)",s):
        return ["H-bond donor","H-bond acceptor"],"explicit hydroxyl outer surface"
    if re.search(r"(phosphonate|phosphate[- ]?terminated)",s):
        return ["phosphate-rich","anionic","H-bond acceptor"],"explicit phosphate/phosphonate outer surface"
    if re.search(r"(polyethylene glycol|\bpeg\b)",s):
        return ["H-bond acceptor"],"PEG outer surface"
    if re.search(r"(methyl[- ]?thiol|methyl[- ]?terminated|\bch3\b)",s):
        return ["hydrophobic"],"explicit methyl outer surface"

    # bare/material identities
    if s in {"gold","au","bare gold","unmodified gold"}:
        return ["soft-metal sulfur affinity"],"bare gold"
    if re.search(r"(silica|silicon oxide|quartz|\bglass\b)",s):
        return ["metal-oxide"],"silica/silicon oxide class"
    if re.search(r"(titanium oxide|titania|tio2)",s):
        return ["metal-oxide"],"titania class"
    if re.search(r"(iron oxide|magnetite|fe3o4)",s):
        return ["metal-oxide"],"iron oxide class"
    if re.search(r"(alumina|aluminum oxide|aluminium oxide)",s):
        return ["metal-oxide"],"alumina class"
    if re.search(r"(hydroxyapatite|calcium phosphate)",s):
        return ["calcium/phosphate charged sites"],"calcium phosphate class"
    if re.search(r"(poly\(styrene\)|polystyrene|\bps\b)",s):
        return ["hydrophobic","pi/aromatic"],"polystyrene class"
    if re.search(r"(poly\(tetrafluoroethylene\)|tetrafluorethylene|ptfe)",s):
        return ["hydrophobic"],"PTFE fluorocarbon class"
    if re.search(r"(graphite|graphene|carbon nanotube|\bcnt\b)",s):
        return ["hydrophobic","pi/aromatic"],"graphitic carbon class"
    if re.search(r"(diamond-like carbon|\bdlc\b)",s):
        return ["hydrophobic","pi/aromatic"],"diamond-like carbon class"
    if re.search(r"\bmica\b",s):
        return ["anionic","H-bond acceptor"],"mica"
    return [],"ambiguous/unmapped"

def fetch_tables():
    url="https://molecularsense.com/bad-2/"
    req=urllib.request.Request(url,headers={"User-Agent":UA})
    with urllib.request.urlopen(req,timeout=120) as r:
        html=r.read().decode("utf-8",errors="replace")
    (OUT/"bad2_snapshot.html").write_text(html,encoding="utf-8")
    tabs=pd.read_html(io.StringIO(html))
    # first table with adsorption columns
    for df in tabs:
        cols=[str(c).strip() for c in df.columns]
        if "Surface Concentration" in cols and "Adsorbing Surface" in cols and "PDB" in cols:
            df.columns=cols
            return df, url
    raise RuntimeError("BAD adsorption table not found")

def fnum(x):
    try:
        if pd.isna(x): return None
        s=str(x).strip()
        if not s or s.lower() in {"nan","none","na","n/a"}: return None
        return float(s)
    except:
        return None

def fetch_pdb(pdb):
    url=f"https://files.rcsb.org/download/{pdb}.pdb"
    req=urllib.request.Request(url,headers={"User-Agent":UA})
    with urllib.request.urlopen(req,timeout=60) as r:
        return r.read().decode("utf-8",errors="replace")

def compute_gsc(pdb,pH,channels):
    raw=fetch_pdb(pdb)
    prepared,prep=prepare_pdb_text(raw,chain=None)
    req=core.AnalyzeRequest(
        pdb_text=prepared,
        chain=None,
        env=core.EnvParams(pH=float(pH),ionic=150.0,temp=298.0)
    )
    res=core.analyze(req)
    if hasattr(res,"body"):
        raise RuntimeError("core returned HTTP response")
    surface=res.get("surface_residues",[])
    denom=sum(float(r.get("scrsa") or 0.0) for r in surface)
    bykey={str(r.get("key")):0.0 for r in surface}
    for pub in channels:
        internal=MAP_TO_INTERNAL[pub]
        ch=res.get("chemistries",{}).get(internal,{})
        for rr in ch.get("residues",[]):
            k=str(rr.get("key"))
            if k in bykey:
                bykey[k]=max(bykey[k],float(rr.get("local_score") or 0.0))
    numer=sum(bykey.values())
    compatible=sum(1 for v in bykey.values() if v>0)
    n_surface=len(surface)

    # predefined diagnostics
    # compatible residue fraction and raw compatible scRSA-like support
    comp_frac=compatible/n_surface if n_surface else float("nan")

    # derive patches from public engine-level implementation for diagnostic descriptors
    from interface_engine import analyze_all_maps
    allmaps=analyze_all_maps(pdb_text=raw,chain=None,pH=float(pH),ionic_mM=150.0,temp_K=298.0)
    max_support=0.0; max_coh=0.0; n_front1=0; max_prop=0.0
    for pub in channels:
        public_key={
            "anionic":"anionic","cationic":"cationic","hydrophobic":"hydrophobic",
            "pi/aromatic":"pi_carbon","H-bond donor":"hbond_donor","H-bond acceptor":"hbond_acceptor",
            "metal-oxide":"oxide","calcium/phosphate charged sites":"calcium_phosphate",
            "transition-metal coordination":"metal_coord","soft-metal sulfur affinity":"soft_metal_sulfur",
            "phosphate-rich":"phosphate"
        }[pub]
        mm=allmaps["maps"].get(public_key,{})
        f1=[p for p in mm.get("patches",[]) if int(p.get("pareto_front",999))==1]
        n_front1 += len(f1)
        if f1:
            max_support=max(max_support,max(float(p.get("chemistry_support") or 0) for p in f1))
            max_coh=max(max_coh,max(float(p.get("patch_coherence") or 0) for p in f1))
        if mm.get("residues"):
            max_prop=max(max_prop,max(float(r.get("propensity") or 0) for r in mm["residues"]))
    patch_density=n_front1/n_surface if n_surface else float("nan")
    return {
        "GSC":numer/denom if denom>0 else float("nan"),
        "compatible_residue_fraction":comp_frac,
        "compatible_support_sum":numer,
        "max_front1_chemistry_support":max_support,
        "max_front1_patch_coherence":max_coh,
        "front1_patch_density":patch_density,
        "max_residue_propensity":max_prop,
        "n_surface_residues":n_surface,
        "surface_scrsa_sum":denom,
        "selected_chain":prep.get("selected_chain","ALL"),
    }

def exact_key(v):
    if v is None: return ""
    if isinstance(v,float):
        return f"{v:.12g}"
    return str(v).strip()

def main():
    df,url=fetch_tables()
    df.to_csv(OUT/"bad2_live_snapshot.csv",index=False)

    cache={}
    rows=[]; exclusions=[]
    for _,r in df.iterrows():
        rid=str(r.get("ID","")).strip()
        protein=str(r.get("Protein","")).strip()
        pdb=str(r.get("PDB","")).strip().upper()
        surf=str(r.get("Adsorbing Surface","")).strip()
        gamma=fnum(r.get("Surface Concentration"))
        csol=fnum(r.get("Solution Concentration"))
        ph=fnum(r.get("pH"))
        ionic=fnum(r.get("Ionic Strength"))
        temp=fnum(r.get("Temperature"))
        method=str(r.get("Measurement Method","")).strip()
        etype=str(r.get("Experiment Type","")).strip()
        validation=str(r.get("Data Validation",r.get("Data validation",""))).strip()
        notes=str(r.get("Explanatory Notes","")).strip()
        channels,rule=map_surface(surf)
        reasons=[]
        if not re.fullmatch(r"[0-9A-Za-z]{4}",pdb): reasons.append("invalid_or_missing_pdb")
        if gamma is None or gamma<=0: reasons.append("nonpositive_or_missing_surface_concentration")
        if csol is None or csol<=0: reasons.append("nonpositive_or_missing_solution_concentration")
        if ph is None: reasons.append("missing_pH")
        if not channels: reasons.append("surface_unmapped")
        if re.search(r"(uncertain|invalid)",validation+" "+notes,re.I): reasons.append("source_flag_uncertain_or_invalid")
        base={
            "source_row_id":rid,"protein":protein,"pdb_id":pdb,
            "surface_concentration_mg_m2":gamma if gamma is not None else "",
            "solution_concentration_mg_ml":csol if csol is not None else "",
            "adsorbing_surface":surf,"pH":ph if ph is not None else "",
            "ionic_strength":ionic if ionic is not None else "",
            "temperature_C":temp if temp is not None else "",
            "measurement_method":method,"experiment_type":etype,
            "reference":str(r.get("Reference","")).strip(),"DOI":str(r.get("DOI","")).strip(),
            "year":str(r.get("Year","")).strip(),"data_validation":validation,
            "explanatory_notes":notes,"mapped_channels":";".join(channels),
            "surface_mapping_rule":rule,
        }
        if reasons:
            exclusions.append({**base,"exclusion_reasons":";".join(reasons)})
            continue
        key=(pdb,float(ph),tuple(channels))
        if key not in cache:
            try:
                cache[key]=compute_gsc(pdb,float(ph),channels)
            except Exception as e:
                cache[key]={"score_error":repr(e)}
        sc=cache[key]
        if "score_error" in sc:
            exclusions.append({**base,"exclusion_reasons":"v1_or_structure_failure","score_error":sc["score_error"]})
            continue
        rows.append({**base,**sc})

    pd.DataFrame(rows).to_csv(OUT/"bad2_scored_rows.csv",index=False)
    pd.DataFrame(exclusions).to_csv(OUT/"bad2_exclusions.csv",index=False)

    # Primary exact-comparability strata
    descriptors=[
        "GSC","compatible_residue_fraction","compatible_support_sum",
        "max_front1_chemistry_support","max_front1_patch_coherence",
        "front1_patch_density","max_residue_propensity"
    ]
    groups=defaultdict(list)
    for r in rows:
        k=(
            norm_surface(r["adsorbing_surface"]),
            exact_key(r["solution_concentration_mg_ml"]),
            exact_key(r["pH"]),
            exact_key(r["ionic_strength"]),
            exact_key(r["temperature_C"]),
            r["measurement_method"].strip().lower(),
            r["experiment_type"].strip().lower(),
        )
        groups[k].append(r)

    strata=[]
    for k,arr in groups.items():
        proteins=sorted(set(x["pdb_id"] for x in arr))
        if len(proteins)<3:
            continue
        for desc in descriptors:
            vals=[]; gam=[]
            for x in arr:
                v=x.get(desc)
                try: vf=float(v); gf=float(x["surface_concentration_mg_m2"])
                except: continue
                if np.isfinite(vf) and np.isfinite(gf):
                    vals.append(vf); gam.append(gf)
            if len(set(x["pdb_id"] for x in arr))<3 or len(vals)<3 or len(set(vals))<2 or len(set(gam))<2:
                continue
            rho,p=spearmanr(vals,gam)
            if not np.isfinite(rho): continue
            strata.append({
                "surface":k[0],"solution_concentration_mg_ml":k[1],"pH":k[2],
                "ionic_strength":k[3],"temperature_C":k[4],
                "measurement_method":k[5],"experiment_type":k[6],
                "descriptor":desc,"n_rows":len(vals),"n_unique_pdb":len(set(x["pdb_id"] for x in arr)),
                "spearman_rho":float(rho),"spearman_p":float(p) if np.isfinite(p) else "",
                "protein_pdbs":";".join(proteins),
            })

    pd.DataFrame(strata).to_csv(OUT/"bad2_stratum_correlations.csv",index=False)

    summary=[]
    for desc in descriptors:
        vals=[float(x["spearman_rho"]) for x in strata if x["descriptor"]==desc]
        if vals:
            p=""
            if any(abs(v)>1e-15 for v in vals):
                try:p=float(wilcoxon(np.asarray(vals),zero_method="wilcox",alternative="two-sided").pvalue)
                except: p=""
            summary.append({
                "descriptor":desc,"n_evaluable_strata":len(vals),
                "median_rho":float(np.median(vals)),
                "q1_rho":float(np.quantile(vals,.25)),
                "q3_rho":float(np.quantile(vals,.75)),
                "fraction_rho_positive":float(np.mean(np.asarray(vals)>0)),
                "wilcoxon_p_vs_zero":p
            })
    pd.DataFrame(summary).to_csv(OUT/"bad2_descriptor_summary.csv",index=False)

    audit={
        "source_url":url,
        "live_rows":int(len(df)),
        "eligible_scored_rows":len(rows),
        "excluded_rows":len(exclusions),
        "unique_scored_pdb":len(set(r["pdb_id"] for r in rows)),
        "unique_scored_surfaces":len(set(r["adsorbing_surface"] for r in rows)),
        "unique_score_calculations":len([v for v in cache.values() if "score_error" not in v]),
        "failed_score_calculations":len([v for v in cache.values() if "score_error" in v]),
        "primary_comparable_groups_total":sum(1 for arr in groups.values() if len(set(x["pdb_id"] for x in arr))>=3),
        "stratum_descriptor_rows":len(strata),
        "primary_GSC_strata":sum(x["descriptor"]=="GSC" for x in strata),
        "experimental_adsorption_values_used_in_GSC_definition":False,
        "v1_parameters_refit":False,
    }
    (OUT/"bad2_audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")
    print(json.dumps(audit,indent=2))
    print(json.dumps([x for x in summary if x["descriptor"]=="GSC"],indent=2))

if __name__=="__main__":
    main()
