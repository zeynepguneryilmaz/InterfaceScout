"""Run actual 48 frozen assembly-1 structures with canonical and physics-aware maps.

The 48-PDB identifier set was selected prior to electrostatic scoring.
No results are labeled PASS unless PDB2PQR, PROPKA, APBS and the new map
analysis all completed with actual PQR/DX data and original assembly hashes.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,os,sys,traceback,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"backend"),str(ROOT/"validation")]
from strict_apo_compat_48 import prepare
from interface_engine import analyze_all_maps
from optional_electrostatics import ElectrostaticSettings,compute
from physics_aware_11maps import build_physical_11,comparison_rows
PANEL=ROOT/"validation"/"strict_apo_48_selected_candidates_v2_final.csv"
PDBS=list(csv.DictReader(PANEL.open(encoding="utf-8")))
assert len(PDBS)==48 and len({r["pdb_id"] for r in PDBS})==48

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--shard",type=int,default=0)
    ap.add_argument("--n-shards",type=int,default=6)
    ap.add_argument("--smoke",action="store_true")
    ap.add_argument("--out",default="physics_aware_rerun_48")
    args=ap.parse_args()
    selected=PDBS[args.shard::args.n_shards]
    if args.smoke:selected=[next(x for x in PDBS if x["pdb_id"]=="1O83")]
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    case_rows=[];maps=[]
    for idx,item in enumerate(selected,1):
        pdb=item["pdb_id"];record={"pdb":pdb,"status":"FAIL","shard":args.shard}
        try:
            provenance=prepare(pdb)
            prepared=Path("strict_apo_panel_48")/(pdb+"_prepared.pdb")
            inputtext=prepared.read_text(encoding="utf-8")
            record["prepared_sha256"]=hashlib.sha256(inputtext.encode()).hexdigest()
            baseline=analyze_all_maps(pdb_text=inputtext,pH=7.4,ionic_mM=150,temp_K=298.15)
            if baseline["n_maps"]!=11:raise ValueError("Baseline did not produce 11 maps")
            settings=ElectrostaticSettings(pH=7.4,ionic_mM=150,temperature_K=298.15)
            physical_result,blob=compute(pdb_text=inputtext,chain=None,settings=settings)
            folder=out/pdb;folder.mkdir(exist_ok=True)
            (folder/"external_APBS_PROPKA_PQR.zip").write_bytes(blob)
            with zipfile.ZipFile(io.BytesIO(blob)) as z:
                pqr=z.read("data/prepared.pqr").decode()
                dx=z.read("data/potential.dx").decode()
            informed=build_physical_11(prepared_protein_pdb=inputtext,
                         pqr_text=pqr,dx_text=dx,pH=7.4,temperature_K=298.15,
                         baseline_analysis=baseline)
            permap=comparison_rows(baseline,informed,pdb)
            if len(permap)!=11:raise ValueError("Expected 11 compared channels")
            maps.extend(permap)
            (folder/"11_map_model_comparison.json").write_text(json.dumps({
               "pdb":pdb,"source_assembly_sha256":item["assembly_sha256"],
               "prepared_sha256":record["prepared_sha256"],
               "baseline_model":baseline["version"],
               "physics_model":informed["model"],
               "n_maps":11,"physical_surface":informed["physics"],
               "comparisons":permap},indent=2),encoding="utf-8")
            record.update({"status":"PASS","n_maps":11,
                    "sampled_potential_residues":informed["physics"]["potential_sampled_residues"],
                    "n_surface_residues":informed["physics"]["surface_residues"],
                    "external_zip_bytes":len(blob)})
            print("RERUN_48",idx,pdb,"PASS",flush=True)
        except Exception as exc:
            record["error"]=type(exc).__name__+": "+str(exc)
            (out/f"{pdb}_FAIL.log").write_text(traceback.format_exc(),encoding="utf-8")
            print("RERUN_48",idx,pdb,"FAIL",record["error"][:300],flush=True)
        case_rows.append(record)
        (out/"partial_status.json").write_text(json.dumps(case_rows,indent=2))
    def save_csv(path,records):
        headers=sorted({k for row in records for k in row})
        with path.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=headers);w.writeheader();w.writerows(records)
    save_csv(out/"case_status.csv",case_rows)
    if maps:save_csv(out/"eleven_map_comparisons.csv",maps)
    passed=sum(r["status"]=="PASS" for r in case_rows)
    (out/"manifest.json").write_text(json.dumps({
       "frozen_panel":str(PANEL),"all_48_ids":[r["pdb_id"] for r in PDBS],
       "shard":args.shard,"n_shards":args.n_shards,"smoke":args.smoke,
       "n_selected":len(selected),"n_pass":passed,
       "n_failed":len(selected)-passed,
       "experimental_adsorption_validation":False,
       "apbs_solution":"linearized Poisson-Boltzmann",
       "do_not_pool_protein_48_with_evaluation_40":True,
       "models_compared":"original 2.4.0 vs versioned 2.5.0 prototype"},indent=2))
    print("RERUN_48_SUMMARY",passed,"/",len(selected),flush=True)
    if passed!=len(selected):raise SystemExit(1)

if __name__=="__main__":main()
