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
from interface_engine import _prepare_shared_context, _build_map

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
    ctx=_prepare_shared_context(
        pH=float(ph), ionic_mM=150.0, temp_K=298.0,
        pdb_id=pdb, pdb_text=None, chain=None
    )
    cres=ctx["core_result"]
    prepared=ctx["prepared"]
    prep=ctx["prep_report"]
    surface=cres.get("surface_residues",[])
    denom=sum(float(r.get("scrsa") or 0.0) for r in surface)
    scrsa={str(r["key"]):float(r.get("scrsa") or 0.0) for r in surface}
    fav={k:0.0 for k in scrsa}
    rep={k:0.0 for k in scrsa}

    for pub in channels:
        internal=MAP_TO_INTERNAL[pub]
        ch=cres.get("chemistries",{}).get(internal,{})
        for rr in ch.get("residues",[]):
            k=str(rr.get("key"))
            if k in fav:
                fav[k]=max(fav[k],float(rr.get("local_score") or 0.0))
        for rr in ch.get("repulsive_residues",[]):
            k=str(rr.get("key"))
            if k in rep:
                rn=str(rr.get("res_name") or "").upper()
                defs=core.CHEMISTRIES[internal].get("repulsive",{})
                if rn in defs:
                    state=defs[rn][2]
                    rep[k]=max(rep[k],scrsa[k]*core.state_availability(rn,state,ph))

    ngc=(sum(fav.values())-sum(rep.values()))/denom if denom>0 else np.nan
    gfc=sum(fav.values())/denom if denom>0 else np.nan
    grc=sum(rep.values())/denom if denom>0 else np.nan

    patches=[]
    spec_lookup={
        "anionic":("anionic","Anionic surface",""),
        "cationic":("cationic","Cationic surface",""),
        "hydrophobic":("hydrophobic","Hydrophobic surface",""),
        "pi/aromatic":("pi_carbon","π / aromatic surface",""),
        "H-bond donor":("hbond_donor","H-bond donor surface",""),
        "H-bond acceptor":("hbond_acceptor","H-bond acceptor surface",""),
        "metal-oxide":("oxide","Metal-oxide surface",""),
        "calcium/phosphate charged sites":("hydroxyapatite","Calcium/phosphate charged sites",""),
        "transition-metal coordination":("metal_coord","Transition-metal coordination",""),
        "soft-metal sulfur affinity":("gold","Soft-metal sulfur affinity",""),
        "phosphate-rich":("phosphate","Phosphate-rich surface",""),
    }
    for pub in channels:
        internal,label,desc=spec_lookup[pub]
        mm=_build_map(ctx, PUBLIC_KEY[pub], internal, label, desc)
        for p in mm.get("patches",[]):
            if int(p.get("pareto_front",999))!=1: continue
            members=[str(x) for x in p.get("members",[])]
            pd=sum(scrsa.get(k,0.0) for k in members)
            if pd<=0: continue
            pf=sum(fav.get(k,0.0) for k in members)
            pr=sum(rep.get(k,0.0) for k in members)
            net=(pf-pr)/pd
            favcomp=pf/pd
            cov=pd/denom if denom>0 else np.nan
            patches.append((net,favcomp,cov,PUBLIC_KEY[pub],str(p.get("center_key",""))))

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
        "GFC":gfc,"GRC":grc,"NGC":ngc,"BPNC":bpnc,"BPFC":bpfc,"BPSC":bpsc,"LNS":lns,
        "n_front1_patches":n_front1,"structure_mw_Da":mw,
        "n_surface_residues":len(surface),"surface_scrsa_sum":denom,
        "selected_chain":prep.get("selected_chain","ALL"),
    }

