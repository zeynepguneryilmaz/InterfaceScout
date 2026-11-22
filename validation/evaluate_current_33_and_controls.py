"""Technical InterfaceScout 11-map rerun for the CURRENT 33-PDB selection and three controls.
This verifies map computation, NOT experimental adsorption, state, or biological validity.
"""
import argparse, csv, json, subprocess, sys
from pathlib import Path
import audit_literature_33_engine_20261122 as previous

CANONICAL = """1A5R 1AKI 1AML 1AO6 1BEB 1F6R 1FNF 1H3X 1HRC 1HZH 1IGT 1JNJ 1MBN 1TGU 1UBQ 2LZM 2CBA 2GB1 2HAV 2HHB 2OED 2PTN 2VH7 4EUL 3ECA 3GHG 3STN 4F5S 5J17 6O5E 7KWI 7RSA 8ONG""".split()
CONTROLS = ["1F6S", "3V83", "1JY3"]
IDS = CANONICAL + CONTROLS
assert len(CANONICAL) == 33 and len(set(IDS)) == 36
OUT = Path("selected_33_results")
OUT.mkdir(exist_ok=True)
previous.OUT = OUT

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--worker")
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--shards", type=int, default=6)
    args=p.parse_args()
    if args.worker:
        try:
            data=previous.work(args.worker)
            assert data["n_maps"] == len(data["n_patches_per_map"]) == 11
            print("RESULT_JSON " + json.dumps(data), flush=True)
        except Exception:
            import traceback
            traceback.print_exc()
            sys.exit(1)
        return
    assert args.shards > 0 and 0 <= args.shard < args.shards
    rows=[]
    for pid in IDS[args.shard::args.shards]:
        row={"pdb":pid,"cohort":"CURRENT_CANONICAL" if pid in CANONICAL else "CONTEXT_COMPARATOR",
             "baseline":"pH7.4,150mM,298K; not necessarily experiment-matched",
             "disclaimer":"Protein-only, first model, altloc collapsed, all assembly chains; technical compatibility ONLY"}
        try:
            row.update(previous.prepare(pid))
            proc=subprocess.run([sys.executable,__file__,"--worker",pid],capture_output=True,text=True,timeout=900)
            if proc.returncode or "RESULT_JSON " not in proc.stdout:
                raise RuntimeError((proc.stderr or proc.stdout or "no result")[-2400:])
            row.update(json.loads(proc.stdout.split("RESULT_JSON ",1)[1].splitlines()[0]))
            row["status"]="PASS_11MAP_TECHNICAL_ONLY"
        except Exception as e:
            row["status"]="FAIL_PREP_OR_ENGINE"
            row["error"]=str(e)
        print(pid,row["status"],row.get("n_residues",""),row.get("error","")[:80],flush=True)
        rows.append(row)
    dst=OUT / f"shard_{args.shard}"
    dst.mkdir(parents=True,exist_ok=True)
    (dst/"results.json").write_text(json.dumps(rows,indent=2),encoding="utf-8")
    cols=["pdb","cohort","status","n_maps","n_residues","n_surface","n_chains",
          "n_patches_per_map","combined_checks","original_sha256","protein_only_sha256",
          "baseline","disclaimer","preparation_impact","error"]
    with (dst/"results.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=cols)
        writer.writeheader()
        for row in rows:
            writer.writerow({k:json.dumps(row.get(k),ensure_ascii=False) if isinstance(row.get(k),(list,dict)) else row.get(k,"") for k in cols})
    if len(rows)!=len(IDS[args.shard::args.shards]):
        raise RuntimeError("incomplete output")

if __name__=="__main__":
    main()
