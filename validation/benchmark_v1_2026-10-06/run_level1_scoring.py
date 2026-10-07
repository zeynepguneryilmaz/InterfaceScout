from __future__ import annotations
import csv, io, json, math, os, re, sys, time, unicodedata, urllib.parse, urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, spearmanr, wilcoxon

ROOT = Path(__file__).resolve().parents[2]
BENCH = ROOT / "validation" / "benchmark_v1_2026-10-06"
OUT = BENCH / "scoring_v1_2026-10-07"
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / "backend"))

from interface_engine import analyze_all_maps

PUBLIC_CHANNELS = [
    "anionic","cationic","hydrophobic","pi_carbon","hbond_donor","hbond_acceptor",
    "oxide","calcium_phosphate","metal_coord","soft_metal_sulfur","phosphate"
]
MAP_COL_TO_PUBLIC = {
    "anionic":"anionic",
    "cationic":"cationic",
    "hydrophobic":"hydrophobic",
    "pi_aromatic":"pi_carbon",
    "hbond_donor":"hbond_donor",
    "hbond_acceptor":"hbond_acceptor",
    "metal_oxide":"oxide",
    "calcium_phosphate_charged_sites":"calcium_phosphate",
    "transition_metal_coordination":"metal_coord",
    "soft_metal_sulfur_affinity":"soft_metal_sulfur",
    "phosphate_rich":"phosphate",
}
DESCRIPTORS = [
    "max_front1_chemistry_support",
    "max_front1_patch_coherence",
    "max_front1_mean_accessibility",
    "max_front1_orientation_coherence",
    "n_front1_patches",
    "max_residue_propensity",
]

UA = "InterfaceScout-benchmark/1.0 (academic reproducibility; GitHub Actions)"

def http_get(url, timeout=90, retries=4):
    last=None
    for attempt in range(retries):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":UA,"Accept":"application/json,text/plain,*/*"})
            with urllib.request.urlopen(req,timeout=timeout) as r:
                return r.read()
        except Exception as e:
            last=e
            time.sleep(min(8,1.5**attempt))
    raise last

def norm_name(s):
    s=unicodedata.normalize("NFKD",s or "")
    s=s.replace("β","beta").replace("α","alpha").replace("γ","gamma").replace("κ","kappa").replace("λ","lambda")
    s=s.casefold()
    s=re.sub(r"[^a-z0-9]+"," ",s)
    return " ".join(s.split())

def primary_uniprot_name(s):
    # UniProt TSV "Protein names" places the recommended name first, followed by parenthesized alternatives.
    x=(s or "").strip()
    cut=x.find(" (")
    if cut>0: x=x[:cut]
    return x.strip()

def download_human_uniprot():
    query=urllib.parse.urlencode({
        "query":"(organism_id:9606)",
        "format":"tsv",
        "fields":"accession,reviewed,protein_name,length",
    })
    data=http_get("https://rest.uniprot.org/uniprotkb/stream?"+query,timeout=180).decode("utf-8",errors="replace")
    rows=list(csv.DictReader(io.StringIO(data),delimiter="\t"))
    by=defaultdict(list)
    for r in rows:
        pn=primary_uniprot_name(r.get("Protein names",""))
        k=norm_name(pn)
        if k:
            by[k].append({
                "accession":r.get("Entry",""),
                "reviewed":(r.get("Reviewed","").strip().lower()=="reviewed"),
                "protein_name_uniprot":pn,
                "length":r.get("Length",""),
            })
    return by, len(rows)

def resolve_pcdb_names(protein_rows):
    by,total=download_human_uniprot()
    out=[]
    for r in protein_rows:
        if r.get("source_database")!="Chou_PCDB": continue
        k=norm_name(r.get("protein_name",""))
        cands=by.get(k,[])
        status=""
        chosen=None
        if len(cands)==1:
            chosen=cands[0]; status="unique_exact_normalized"
        elif len(cands)>1:
            rev=[x for x in cands if x["reviewed"]]
            if len(rev)==1:
                chosen=rev[0]; status="unique_reviewed_exact_normalized"
            else:
                status="ambiguous_exact_normalized"
        else:
            status="no_exact_normalized_match"
        out.append({
            "source_database":"Chou_PCDB",
            "source_identifier":r.get("source_identifier",""),
            "protein_name":r.get("protein_name",""),
            "species":r.get("species","Human"),
            "uniprot_accession":chosen["accession"] if chosen else "",
            "identity_resolution_status":status,
            "identity_candidate_count":len(cands),
            "uniprot_name":chosen["protein_name_uniprot"] if chosen else "",
            "uniprot_reviewed":chosen["reviewed"] if chosen else "",
            "uniprot_length":chosen["length"] if chosen else "",
        })
    return out,total

def resolve_payne(protein_rows):
    out=[]
    for r in protein_rows:
        if r.get("source_database")!="Payne_PXD053700": continue
        acc=(r.get("uniprot_accession") or "").strip()
        notes=r.get("notes","")
        multi=("Ambiguous protein group contains" in notes) or not acc
        out.append({
            "source_database":"Payne_PXD053700",
            "source_identifier":r.get("source_identifier",""),
            "protein_name":r.get("protein_name",""),
            "species":r.get("species","Bos taurus"),
            "uniprot_accession":acc if not multi else "",
            "identity_resolution_status":"source_single_accession" if not multi else "source_multi_accession_unresolved",
            "identity_candidate_count":1 if not multi else "",
            "uniprot_name":"",
            "uniprot_reviewed":"",
            "uniprot_length":r.get("length",""),
        })
    return out


def _select_pdbe_candidate(arr):
    if not isinstance(arr,list):
        return None
    cand=[]
    for x in arr:
        try: cov=float(x.get("coverage"))
        except: cov=0.0
        res=x.get("resolution")
        try: resf=float(res) if res not in (None,"","None") else None
        except: resf=None
        if cov < 0.70: continue
        if resf is not None and resf>4.0: continue
        pdb=(x.get("pdb_id") or x.get("pdbId") or "").upper()
        chain=(x.get("chain_id") or x.get("chainId") or "").strip()
        if not pdb or not chain: continue
        cand.append((cov, resf, pdb, chain, x))
    if not cand:
        return None
    cand.sort(key=lambda z:(-z[0], z[1] if z[1] is not None else 999.0, z[2], z[3]))
    cov,res,pdb,chain,x=cand[0]
    return {
        "structure_origin":"experimental_pdb",
        "pdb_id":pdb,"chain":chain,"coverage":cov,"resolution":res,
        "structure_url":f"https://files.rcsb.org/download/{pdb}.pdb",
    }

def batch_pdbe_best(accessions, chunk_size=200):
    accessions=sorted(set(a for a in accessions if a))
    out={}
    endpoint="https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/"
    for i in range(0,len(accessions),chunk_size):
        chunk=accessions[i:i+chunk_size]
        body=",".join(chunk).encode("utf-8")
        req=urllib.request.Request(endpoint,data=body,headers={
            "User-Agent":UA,
            "Accept":"application/json",
            "Content-Type":"text/plain",
        },method="POST")
        last=None
        for attempt in range(4):
            try:
                with urllib.request.urlopen(req,timeout=120) as r:
                    obj=json.loads(r.read().decode("utf-8"))
                for acc in chunk:
                    arr=obj.get(acc) or obj.get(acc.upper()) or obj.get(acc.lower()) or []
                    out[acc]=_select_pdbe_candidate(arr)
                last=None
                break
            except Exception as e:
                last=e
                time.sleep(min(8,1.5**attempt))
        if last is not None:
            # Deterministic fallback to the already-defined single-accession endpoint.
            for acc in chunk:
                exp,_=get_best_pdbe_structure(acc)
                out[acc]=exp
    return out

def enrich_identities_batch(identities):
    accessions=sorted(set((r.get("uniprot_accession") or "").strip() for r in identities if (r.get("uniprot_accession") or "").strip()))
    pdbe=batch_pdbe_best(accessions)
    no_exp=[a for a in accessions if not pdbe.get(a)]
    af_map={}
    with ThreadPoolExecutor(max_workers=16) as ex:
        fut={ex.submit(get_alphafold_structure,a):a for a in no_exp}
        for f in as_completed(fut):
            a=fut[f]
            try:
                af,status=f.result()
            except Exception:
                af,status=None,"alphafold_exception"
            af_map[a]=(af,status)
    enriched=[]
    for r in identities:
        z=dict(r)
        acc=(r.get("uniprot_accession") or "").strip()
        if not acc:
            z.update({"structure_status":"identity_unresolved","structure_origin":"","pdb_id":"","chain":"","coverage":"","resolution":"","structure_url":""})
        elif pdbe.get(acc):
            z.update(pdbe[acc]); z["structure_status"]="experimental_selected"
        else:
            af,status=af_map.get(acc,(None,"alphafold_no_model"))
            if af:
                z.update(af); z["structure_status"]="alphafold_selected"; z["pdbe_status"]="pdbe_no_eligible_structure"
            else:
                z.update({"structure_status":"no_qualifying_structure","structure_origin":"","pdb_id":"","chain":"","coverage":"","resolution":"","structure_url":"","pdbe_status":"pdbe_no_eligible_structure","alphafold_status":status})
        enriched.append(z)
    return enriched

def get_best_pdbe_structure(accession):
    url=f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{urllib.parse.quote(accession)}"
    try:
        obj=json.loads(http_get(url,timeout=60).decode("utf-8"))
    except Exception:
        return None, "pdbe_unavailable_or_no_mapping"
    arr=obj.get(accession) or obj.get(accession.upper()) or obj.get(accession.lower())
    if arr is None and len(obj)==1:
        arr=next(iter(obj.values()))
    if not isinstance(arr,list):
        return None,"pdbe_no_mapping"
    cand=[]
    for x in arr:
        try: cov=float(x.get("coverage"))
        except: cov=0.0
        res=x.get("resolution")
        try: resf=float(res) if res not in (None,"","None") else None
        except: resf=None
        if cov < 0.70: continue
        if resf is not None and resf>4.0: continue
        pdb=(x.get("pdb_id") or x.get("pdbId") or "").upper()
        chain=(x.get("chain_id") or x.get("chainId") or "").strip()
        if not pdb or not chain: continue
        cand.append((cov, resf, pdb, chain, x))
    if not cand:
        return None,"pdbe_no_eligible_structure"
    cand.sort(key=lambda z:(-z[0], z[1] if z[1] is not None else 999.0, z[2], z[3]))
    cov,res,pdb,chain,x=cand[0]
    return {
        "structure_origin":"experimental_pdb",
        "pdb_id":pdb,"chain":chain,"coverage":cov,"resolution":res,
        "structure_url":f"https://files.rcsb.org/download/{pdb}.pdb",
    },"experimental_selected"

def get_alphafold_structure(accession):
    url=f"https://alphafold.ebi.ac.uk/api/prediction/{urllib.parse.quote(accession)}"
    try:
        obj=json.loads(http_get(url,timeout=60).decode("utf-8"))
    except Exception as e:
        return None,"alphafold_no_model"
    if isinstance(obj,list) and obj:
        x=obj[0]
    elif isinstance(obj,dict):
        x=obj
    else:
        return None,"alphafold_no_model"
    pdburl=x.get("pdbUrl") or x.get("pdb_url")
    if not pdburl:
        return None,"alphafold_no_pdb_url"
    return {
        "structure_origin":"alphafold_predicted",
        "pdb_id":"",
        "chain":"A",
        "coverage":"",
        "resolution":"",
        "structure_url":pdburl,
        "alphafold_entry":x.get("entryId") or x.get("entry_id") or "",
        "alphafold_version":x.get("latestVersion") or x.get("version") or "",
    },"alphafold_selected"

def enrich_one_identity(r):
    z=dict(r)
    acc=(r.get("uniprot_accession") or "").strip()
    if not acc:
        z.update({"structure_status":"identity_unresolved","structure_origin":"","pdb_id":"","chain":"","coverage":"","resolution":"","structure_url":""})
        return z
    exp,status=get_best_pdbe_structure(acc)
    if exp:
        z.update(exp); z["structure_status"]=status; return z
    af,afstatus=get_alphafold_structure(acc)
    if af:
        z.update(af); z["structure_status"]=afstatus; z["pdbe_status"]=status; return z
    z.update({"structure_status":"no_qualifying_structure","structure_origin":"","pdb_id":"","chain":"","coverage":"","resolution":"","structure_url":"","pdbe_status":status,"alphafold_status":afstatus})
    return z

def extract_channel_descriptors(result):
    out={}
    for ch in PUBLIC_CHANNELS:
        m=result.get("maps",{}).get(ch,{})
        front1=[p for p in m.get("patches",[]) if int(p.get("pareto_front",999))==1]
        out[ch]={
            "max_front1_chemistry_support":max([float(p.get("chemistry_support",0.0)) for p in front1],default=0.0),
            "max_front1_patch_coherence":max([float(p.get("patch_coherence",0.0)) for p in front1],default=0.0),
            "max_front1_mean_accessibility":max([float(p.get("mean_accessibility",0.0)) for p in front1],default=0.0),
            "max_front1_orientation_coherence":max([float(p.get("orientation_coherence",0.0)) for p in front1],default=0.0),
            "n_front1_patches":len(front1),
            "max_residue_propensity":max([float(x.get("propensity") or 0.0) for x in m.get("residues",[])],default=0.0),
        }
    return out

def score_structure(r):
    z=dict(r)
    if r.get("structure_status") not in ("experimental_selected","alphafold_selected"):
        z["scoring_status"]="not_scored_no_structure"
        return z
    try:
        pdb_text=http_get(r["structure_url"],timeout=120).decode("utf-8",errors="replace")
        result=analyze_all_maps(
            pdb_text=pdb_text,
            chain=r.get("chain") or None,
            pH=7.4, ionic_mM=150.0, temp_K=298.0,
        )
        z["scoring_status"]="scored"
        z["v1_version"]=result.get("version","")
        z["n_surface_residues"]=result.get("diagnostics",{}).get("n_surface_residues","")
        z["channel_descriptors_json"]=json.dumps(extract_channel_descriptors(result),separators=(",",":"))
    except Exception as e:
        z["scoring_status"]="failed_structure_or_v1"
        z["scoring_error"]=repr(e)
    return z

def parse_csv(path):
    with path.open(encoding="utf-8-sig",newline="") as f: return list(csv.DictReader(f))

def write_csv(path, rows):
    fields=[]
    for r in rows:
        for k in r:
            if k not in fields: fields.append(k)
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader(); w.writerows(rows)

def mapping_lookup():
    rows=parse_csv(BENCH/"material_to_interfacescout_map.csv")
    out={}
    for r in rows:
        src=r.get("source_database") or "Chou_PCDB"
        chans=[]
        for col,pub in MAP_COL_TO_PUBLIC.items():
            if str(r.get(col,"")).strip() in ("1","1.0","true","True"):
                chans.append(pub)
        out[(src,r.get("core_material","").strip(),r.get("surface_modification","").strip())]=chans
    return out

def composite_desc(channel_json, channels):
    d=json.loads(channel_json)
    out={}
    for name in DESCRIPTORS:
        vals=[float(d.get(ch,{}).get(name,0.0)) for ch in channels]
        if name=="n_front1_patches":
            out[name]=sum(vals)
        else:
            out[name]=max(vals,default=0.0)
    return out

def auroc(y,score):
    y=np.asarray(y,dtype=int); score=np.asarray(score,dtype=float)
    pos=(y==1); neg=(y==0)
    n1=int(pos.sum()); n0=int(neg.sum())
    if n1==0 or n0==0: return None
    ranks=rankdata(score,method="average")
    u=float(ranks[pos].sum()-n1*(n1+1)/2)
    return u/(n1*n0)

def fdr_bh(pvals):
    idx=[i for i,p in enumerate(pvals) if p is not None and np.isfinite(p)]
    out=[None]*len(pvals)
    if not idx: return out
    order=sorted(idx,key=lambda i:pvals[i])
    m=len(order); prev=1.0
    for rank,i in reversed(list(enumerate(order,start=1))):
        q=min(prev,pvals[i]*m/rank); prev=q; out[i]=q
    return out

def summarize_metric_rows(rows, metric_name, null_value):
    groups=defaultdict(list)
    for r in rows:
        try: v=float(r[metric_name])
        except: continue
        groups[(r["source_database"],r["structure_group"],r["descriptor"])].append(v)
    summary=[]
    for (src,origin,desc),vals in sorted(groups.items()):
        arr=np.asarray(vals,float)
        p=None
        if len(arr)>0 and np.any(arr!=null_value):
            try: p=float(wilcoxon(arr-null_value,alternative="two-sided",zero_method="wilcox").pvalue)
            except: p=None
        summary.append({
            "source_database":src,"structure_group":origin,"descriptor":desc,
            "n_evaluable_experiments":len(arr),
            "median":float(np.median(arr)),
            "q1":float(np.quantile(arr,.25)),
            "q3":float(np.quantile(arr,.75)),
            "fraction_above_null":float(np.mean(arr>null_value)),
            "null_value":null_value,
            "wilcoxon_p":p if p is not None else "",
        })
    # FDR within source/origin across six descriptors
    buckets=defaultdict(list)
    for i,r in enumerate(summary):
        buckets[(r["source_database"],r["structure_group"])].append(i)
    for inds in buckets.values():
        pv=[]
        for i in inds:
            try: pv.append(float(summary[i]["wilcoxon_p"]))
            except: pv.append(None)
        q=fdr_bh(pv)
        for ii,qq in zip(inds,q):
            summary[ii]["wilcoxon_fdr_bh"]=qq if qq is not None else ""
    return summary

def main():
    proteins=parse_csv(BENCH/"proteins.csv")
    pcdb, human_uniprot_count=resolve_pcdb_names(proteins)
    payne=resolve_payne(proteins)
    identities=pcdb+payne
    print(f"Identity rows: {len(identities)}; human UniProt downloaded: {human_uniprot_count}",flush=True)

    enriched=enrich_identities_batch(identities)
    print(f"structure enrichment batch complete: {len(enriched)} identities",flush=True)
    enriched.sort(key=lambda r:(r["source_database"],r["source_identifier"]))
    write_csv(OUT/"level1_structure_enrichment.csv",enriched)

    scoreable=[r for r in enriched if r.get("structure_status") in ("experimental_selected","alphafold_selected")]
    scored=[]
    # Structure scoring is CPU-heavy; process-level parallelism changes execution only, not the frozen scientific protocol.
    with ProcessPoolExecutor(max_workers=4) as ex:
        fut={ex.submit(score_structure,r):r for r in scoreable}
        for i,f in enumerate(as_completed(fut),1):
            scored.append(f.result())
            if i%50==0: print(f"InterfaceScout structures {i}/{len(scoreable)}",flush=True)
    # add unscoreable audit rows
    scored += [dict(r,scoring_status="not_scored_no_structure") for r in enriched if r not in scoreable]
    scored.sort(key=lambda r:(r["source_database"],r["source_identifier"]))
    write_csv(OUT/"level1_structure_scores.csv",scored)

    by_pc_name={r["protein_name"]:r for r in scored if r["source_database"]=="Chou_PCDB"}
    by_pay_id={r["source_identifier"]:r for r in scored if r["source_database"]=="Payne_PXD053700"}
    maps=mapping_lookup()

    occurrence=[]
    abundance=[]

    def evaluate_experiment(source, experiment_id, study_id, channels, obs):
        # obs entries: {protein_key, label, value, score_row}
        if not channels: return
        for structure_group in ("experimental_pdb","alphafold_predicted","pooled"):
            if structure_group=="pooled":
                sub=[x for x in obs if x["score_row"] and x["score_row"].get("scoring_status")=="scored"]
            else:
                sub=[x for x in obs if x["score_row"] and x["score_row"].get("scoring_status")=="scored" and x["score_row"].get("structure_origin")==structure_group]
            if len(sub)<10: continue
            comp=[]
            for x in sub:
                comp.append((x, composite_desc(x["score_row"]["channel_descriptors_json"],channels)))
            y=[x[0]["label"] for x in comp]
            for desc in DESCRIPTORS:
                vals=[x[1][desc] for x in comp]
                a=auroc(y,vals)
                if a is not None:
                    occurrence.append({
                        "source_database":source,"study_id":study_id,"experiment_id":experiment_id,
                        "structure_group":structure_group,"descriptor":desc,
                        "n_scored":len(comp),"n_present":sum(y),"n_nonpositive":len(y)-sum(y),
                        "auroc":a,"channels":";".join(channels),
                    })
                pos=[(x[0]["value"],x[1][desc]) for x in comp if x[0]["label"]==1 and x[0]["value"] is not None and x[0]["value"]>0]
                if len(pos)>=10 and len(set(v for v,_ in pos))>=3:
                    rho,p=spearmanr([v for v,_ in pos],[s for _,s in pos])
                    if np.isfinite(rho):
                        abundance.append({
                            "source_database":source,"study_id":study_id,"experiment_id":experiment_id,
                            "structure_group":structure_group,"descriptor":desc,
                            "n_positive_scored":len(pos),"n_distinct_abundance":len(set(v for v,_ in pos)),
                            "spearman_rho":float(rho),"spearman_p_within_experiment":float(p) if np.isfinite(p) else "",
                            "channels":";".join(channels),
                        })

    # PC-DB
    pcmeta={}
    for r in parse_csv(BENCH/"raw_sources/chou_pcdb/PC-DB_for_Meta-Analysis.csv"):
        pcmeta[r["NP Entry ID"].strip()]=r
    with (BENCH/"raw_sources/chou_pcdb/Dataset-for-protein-freq-analysis-v2.csv").open(encoding="utf-8-sig",newline="") as f:
        rr=csv.reader(f); header=next(rr); pnames=header[2:]
        for row in rr:
            if not row: continue
            eid=row[0].strip(); sid=row[1].strip(); m=pcmeta.get(eid,{})
            channels=maps.get(("Chou_PCDB",m.get("Core NP Material","").strip(),m.get("Surface Modification","").strip()),[])
            if not channels: continue
            obs=[]
            for j,pname in enumerate(pnames,start=2):
                raw=row[j].strip() if j<len(row) else ""
                try: val=float(raw)
                except: val=None
                obs.append({"protein_key":pname,"label":1 if val is not None and val>0 else 0,"value":val,"score_row":by_pc_name.get(pname)})
            evaluate_experiment("Chou_PCDB",eid,sid,channels,obs)

    # Payne
    with (BENCH/"raw_sources/payne2025/NP_Database_BovOnly_v5_parsed.csv").open(encoding="utf-8-sig",newline="") as f:
        npmeta={str(r["Sample_num"]).strip():r for r in csv.DictReader(f)}
    with (BENCH/"raw_sources/payne2025/Bov_SP_06252024_Pers.txt").open(encoding="utf-8-sig",newline="") as f:
        rr=csv.reader(f,delimiter="\t"); h=next(rr); types=next(rr); data=list(rr)
    id_idx=h.index("Majority protein IDs")
    intensity=[(i,n) for i,n in enumerate(h) if n.startswith("Intensity ") and not (n.endswith("_FBS") or "_SPQC" in n)]
    for ci,colname in intensity:
        label=colname.replace("Intensity ","",1); sample=label.split("_",1)[0]
        pm=npmeta.get(sample)
        if pm:
            core="Polystyrene" if pm.get("Core_PS")=="1" else ("Iron oxide/Au" if pm.get("Core_Iron_Oxide_Au")=="1" else ("Iron oxide" if pm.get("Core_Iron_Oxide")=="1" else pm.get("Core Material","")))
            mod=pm.get("Ligands","")
        elif "PS Carb@ PEG" in label:
            core="Polystyrene"; mod="Carboxylate_PEG2k"
        else:
            core="";mod=""
        channels=maps.get(("Payne_PXD053700",core,mod),[])
        if not channels: continue
        obs=[]
        for gi,r in enumerate(data,start=1):
            raw=r[ci].strip() if ci<len(r) else ""
            try: val=float(raw)
            except: val=None
            obs.append({"protein_key":f"protein_group_{gi}","label":1 if val is not None and val>0 else 0,"value":val,"score_row":by_pay_id.get(f"protein_group_{gi}")})
        evaluate_experiment("Payne_PXD053700",f"PAYNE_{sample}","PXD053700",channels,obs)

    write_csv(OUT/"level1_occurrence_experiment_metrics.csv",occurrence)
    write_csv(OUT/"level1_abundance_experiment_metrics.csv",abundance)
    occ_summary=summarize_metric_rows(occurrence,"auroc",0.5)
    abu_summary=summarize_metric_rows(abundance,"spearman_rho",0.0)
    write_csv(OUT/"level1_occurrence_summary.csv",occ_summary)
    write_csv(OUT/"level1_abundance_summary.csv",abu_summary)

    # Audit
    def cnt(rows,key,val): return sum(r.get(key)==val for r in rows)
    audit={
        "human_uniprot_entries_downloaded":human_uniprot_count,
        "pcdb_source_entities":len(pcdb),
        "pcdb_identity_resolved":sum(bool(r.get("uniprot_accession")) for r in pcdb),
        "pcdb_identity_unresolved":sum(not bool(r.get("uniprot_accession")) for r in pcdb),
        "payne_source_entities":len(payne),
        "payne_identity_resolved":sum(bool(r.get("uniprot_accession")) for r in payne),
        "payne_identity_unresolved":sum(not bool(r.get("uniprot_accession")) for r in payne),
        "experimental_pdb_selected":cnt(enriched,"structure_origin","experimental_pdb"),
        "alphafold_selected":cnt(enriched,"structure_origin","alphafold_predicted"),
        "no_qualifying_structure":sum(not r.get("structure_origin") for r in enriched),
        "v1_structures_scored":cnt(scored,"scoring_status","scored"),
        "v1_structures_failed":cnt(scored,"scoring_status","failed_structure_or_v1"),
        "occurrence_experiment_metric_rows":len(occurrence),
        "abundance_experiment_metric_rows":len(abundance),
        "mapped_experiments_expected_from_freeze":464,
        "descriptor_count":len(DESCRIPTORS),
        "interface_scout_results_used_for_identity_or_eligibility":False,
    }
    (OUT/"level1_scoring_audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")
    print(json.dumps(audit,indent=2),flush=True)

if __name__=="__main__":
    main()
