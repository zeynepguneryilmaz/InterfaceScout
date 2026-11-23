"""Join the four independently completed ISV80 data audits.

This is a reproducible structural triage, not a final independent validation
panel. Engine-compatibility does not equal biological/coordinate acceptance.
"""
from __future__ import annotations
import csv,json,sys,collections
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
POOL=ROOT/"validation/ISV_protein_candidate_pool_80_v2.csv"
def loadcsv(p):
    with Path(p).open(newline="",encoding="utf-8-sig") as f:return list(csv.DictReader(f))
def num(x,default=None):
    try:return float(x)
    except (ValueError,TypeError):return default
def main(work):
    work=Path(work)
    candidate=loadcsv(POOL)
    coord={r["proposed_pdb_id"]:r for r in loadcsv(work/"isv_coordinate_80.csv")}
    metadata={r["proposed_pdb_id"]:r for r in loadcsv(work/"isv_rcsb_metadata_68.csv")}
    metadata.update({r["proposed_pdb_id"]:r for r in loadcsv(work/"biofilm_12_rcsb_records.csv")})
    detailed=[]
    for i in range(8):
        fn=work/f"details_{i}.json"
        detailed.extend(json.loads(fn.read_text(encoding="utf-8")))
    engine={r["pdb"]:r for r in detailed}
    if any(len(x)!=80 for x in (candidate,coord,metadata,engine)):
        raise ValueError(f"Expected complete 80 candidate records, got {[len(x) for x in (candidate,coord,metadata,engine)]}")
    final=[]
    for entry in candidate:
        pdb=entry["proposed_pdb_id"]
        e=engine[pdb]; c=coord[pdb]; m=metadata[pdb]
        method=m["experimental_method"]
        resolution=num(m.get("resolution_A"))
        coverage=num(c.get("best_chain_CA_coverage"))
        sc_fraction=e.get("sidechain_incomplete_residues",0)/max(e.get("n_residues",1),1)
        warnings=[]
        if coverage is None or coverage<.90:warnings.append("DEPOSITED_CA_COVERAGE_LT90_OR_MISSING")
        if sc_fraction>.10:warnings.append("SIDECHAIN_MISSING_GT10_PERCENT")
        if int(e.get("backbone_incomplete_residues",0))>0:warnings.append("BACKBONE_ATOM_MISSING")
        if "X-RAY" in method.upper() and resolution is not None and resolution>2.5:
            warnings.append("XRAY_RESOLUTION_GT2P5")
        if ("NMR" in method.upper()):warnings.append("NMR_ENSEMBLE_POLICY")
        if ("ELECTRON" in method.upper()):warnings.append("CRYOEM_RESOLUTION_MODEL_REVIEW")
        if int(e.get("n_chains",0))>12 or int(e.get("n_residues",0))>2500:
            warnings.append("LARGE_ASSEMBLY_REVIEW")
        caveat=entry["structural_caveat"].lower()
        if any(term in caveat for term in ("domain","fragment","truncated","partial construct","repeat only")):
            warnings.append("PARTIAL_CONSTRUCT_SCOPE")
        if any(term in caveat for term in ("heme","glycosylation","ligand","zinc","copper","nickel","calcium","metal","cofactor","fad")):
            warnings.append("LIGAND_CHEMICAL_STATE_POLICY")
        if entry["category"] in ("biofilm_membrane_complex","biofilm_fibril") or pdb in ("1HC1",):
            warnings.append("SPECIAL_STRUCTURAL_REPRESENTATION")
        engine_pass=e["engine_status"]=="PASS"
        if not engine_pass: bucket="ENGINE_NOT_QUALIFIED"
        elif "SPECIAL_STRUCTURAL_REPRESENTATION" in warnings:bucket="SPECIAL_CHALLENGE_CONDITIONAL"
        elif warnings:bucket="CONDITIONAL_QA_REVIEW"
        else:bucket="PRIORITY_STRUCTURAL_QA_CANDIDATE"
        final.append({
          "candidate_id":entry["candidate_id"],"pdb":pdb,"protein":entry["protein_name"],
          "category":entry["category"],"engine_pass":int(engine_pass),
          "n_maps":(e.get("engine") or {}).get("n_maps",""),
          "n_maps_with_patches":(e.get("engine") or {}).get("n_nonempty_patch_maps",""),
          "method":method,"resolution_A":resolution if resolution is not None else "",
          "observed_CA_coverage":coverage if coverage is not None else "",
          "n_residues":e.get("n_residues",""),"n_chains":e.get("n_chains",""),
          "missing_backbone_residues":e.get("backbone_incomplete_residues",""),
          "missing_sidechain_residues":e.get("sidechain_incomplete_residues",""),
          "missing_sidechain_fraction":round(sc_fraction,5),
          "structural_priority":bucket,"reasons":";".join(warnings),
          "source_sha256":e.get("assembly_cif_sha256",""),
          "prepared_sha256":e.get("prepared_pdb_sha256",""),
          "experimental_residue_evidence_accepted":"NO_SOURCE_LABEL_AUDIT_YET",
          "TV48_domain_family_independence":"NOT_ASSESSED",
          "ISV_final_membership":"NOT_YET_ACCEPTED"})
    with (work/"ISV80_TRIAGE.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(final[0]));w.writeheader();w.writerows(final)
    counts=dict(collections.Counter(x["structural_priority"] for x in final))
    summary={"n_total":len(final),"engine_pass":sum(x["engine_pass"] for x in final),
        "n_maps_executed":sum(int(x["n_maps"]) for x in final),
        "n_nonempty_patch_maps":sum(int(x["n_maps_with_patches"]) for x in final),
        "observed_ca_coverage_90pct_or_more":sum(num(x["observed_CA_coverage"],-1)>=.90 for x in final),
        "full_core_priority_count":counts.get("PRIORITY_STRUCTURAL_QA_CANDIDATE",0),
        "tier_counts":counts,
        "all_candidates_retained_for_review":True,
        "interpretation":"Every engine-PASS candidate remains potentially usable. Quality tiers are screening priorities only; true final acceptance awaits source/chain/sidechain/biological context and sequence/family matching to TV48.",
        "origin_runs":{"engine":"45270447667","coordinates":"45270701434","metadata_68":"45259307414","biofilm_12":"45268134181"}}
    (work/"ISV80_TRIAGE_SUMMARY.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False))
    print("ISV80_TRIAGE_SUMMARY "+json.dumps(summary),flush=True)
if __name__=="__main__":
    if len(sys.argv)!=2:raise SystemExit("Usage: isv_triage_80.py INPUT_DIR")
    main(sys.argv[1])
