from __future__ import annotations

import json
import re
import urllib.request
from pathlib import Path

from openpyxl import Workbook

from v2.interface_engine import analyze_all_maps

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"

CASES = [
    {"pdb": "4F5S", "chain": "A", "pH": 8.0, "map": "anionic", "folder": "4F5S"},
    {"pdb": "2VUF", "chain": "A", "pH": 6.2, "map": "anionic", "folder": "2VUF"},
    {"pdb": "2CBA", "chain": "A", "pH": 6.3, "map": "anionic", "folder": "2CBA/pH_6.3"},
    {"pdb": "2CBA", "chain": "A", "pH": 8.3, "map": "anionic", "folder": "2CBA/pH_8.3"},
    {"pdb": "1UBQ", "chain": "A", "pH": 7.7, "map": "anionic", "folder": "1UBQ"},
    {"pdb": "1JNJ", "chain": "A", "pH": 7.0, "map": "anionic", "folder": "1JNJ"},
    {"pdb": "3ECA", "chain": "A,B,C,D", "pH": 7.5, "map": "cationic", "folder": "3ECA"},
    {"pdb": "1SUV", "chain": "C,E", "pH": 7.4, "map": "hydrophobic", "folder": "1SUV"},
    {"pdb": "1TGU", "chain": "A,B,C,D", "pH": 7.4, "map": "hydrophobic", "folder": "1TGU"},
]


def fetch_pdb(pdb_id: str) -> str:
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=60) as r:
        return r.read().decode("utf-8", errors="replace")


def flatten(obj, prefix=""):
    rows = []
    if obj is None:
        return [(prefix, "")]
    if isinstance(obj, list):
        return [(prefix, "; ".join(json.dumps(v, separators=(",", ":")) if isinstance(v, (dict, list)) else str(v) for v in obj))]
    if isinstance(obj, dict):
        for k, v in obj.items():
            rows.extend(flatten(v, f"{prefix}.{k}" if prefix else k))
        return rows
    return [(prefix, obj)]


def residue_rows(data):
    all_rows = {}
    for key in data["map_order"]:
        for r in data["maps"][key].get("residues", []):
            x = all_rows.setdefault(r["key"], {
                "ResidueKey": r["key"], "Chain": r.get("chain"), "ResidueNumber": r.get("res_seq"),
                "InsertionCode": r.get("icode", ""), "ResidueName": r.get("res_name"), "scRSA": r.get("scrsa")
            })
            x[f"{key}_local_score"] = r.get("local_score")
            x[f"{key}_propensity_0_100"] = r.get("propensity")
            x[f"{key}_multiscale_persistence"] = r.get("multiscale_persistence")
    return sorted(all_rows.values(), key=lambda r: (str(r.get("Chain") or ""), int(r.get("ResidueNumber") or 0), str(r.get("InsertionCode") or "")))


def patch_rows(data):
    rows = []
    for key in data["map_order"]:
        m = data["maps"][key]
        for p in m.get("patches", []):
            rows.append({
                "MapKey": key, "MapLabel": m["label"], "DisplayRank": p.get("display_rank"),
                "ParetoFront": p.get("pareto_front"), "Classification": p.get("classification", ""),
                "CenterKey": p.get("center_key"), "MemberCount": p.get("n_members"),
                "ChemistrySeedCount": p.get("n_chemistry_seeds"), "Chains": ";".join(p.get("chains", [])),
                "ResidueSpanMin": p.get("residue_span_min"), "ResidueSpanMax": p.get("residue_span_max"),
                "Diameter_A": p.get("diameter_A"), "ChemistrySupport": p.get("chemistry_support"),
                "MeanAccessibility": p.get("mean_accessibility"), "PatchCoherence": p.get("patch_coherence"),
                "OrientationCoherence": p.get("orientation_coherence"), "Members": ";".join(p.get("members", [])),
                "ChemistrySeeds": ";".join(p.get("seed_members", [])),
            })
    return rows


def patch_member_rows(data):
    rows = []
    for key in data["map_order"]:
        m = data["maps"][key]
        for p in m.get("patches", []):
            seeds = set(p.get("seed_members", []))
            for r in p.get("member_residues", []):
                rows.append({
                    "MapKey": key, "MapLabel": m["label"], "PatchRank": p.get("display_rank"),
                    "ParetoFront": p.get("pareto_front"), "ResidueKey": r.get("key", ""),
                    "Chain": r.get("chain"), "ResidueNumber": r.get("res_seq"), "InsertionCode": r.get("icode", ""),
                    "ResidueName": r.get("res_name"), "IsChemistrySeed": r.get("key") in seeds,
                })
    return rows


def add_sheet(wb, name, rows):
    ws = wb.create_sheet(name)
    if not rows:
        ws.append(["Info"]); ws.append(["No data"]); return
    headers = list(rows[0].keys())
    ws.append(headers)
    for row in rows:
        ws.append([row.get(h) for h in headers])
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = ws.dimensions
    for i, h in enumerate(headers, 1):
        max_len = len(str(h))
        for cell in ws.iter_cols(min_col=i, max_col=i, min_row=2, max_row=min(ws.max_row, 201)):
            for c in cell:
                max_len = max(max_len, len(str(c.value or "")))
        ws.column_dimensions[ws.cell(1, i).column_letter].width = min(28, max_len + 2)


def write_excel(data, path):
    wb = Workbook()
    wb.remove(wb.active)
    run = []
    run.append({"Parameter": "engine", "Value": data.get("engine")})
    run.append({"Parameter": "model_version", "Value": data.get("version")})
    run.append({"Parameter": "public_version", "Value": data.get("public_version", "1.0-publication")})
    for section in ["input", "structure_preparation", "method", "scope", "diagnostics", "composite_policy"]:
        for k, v in flatten(data.get(section), section):
            run.append({"Parameter": k, "Value": v})
    add_sheet(wb, "Run_Method", run)
    defs = []
    for key in data["map_order"]:
        m = data["maps"][key]
        defs.append({"MapKey": key, "MapLabel": m["label"], "Description": m["description"],
                     "InternalChemistryKey": m["internal_chemistry_key"], "PatchCount": m["n_patches"],
                     "ParetoFront1Count": m["n_pareto_primary_patches"]})
    add_sheet(wb, "Map_Definitions", defs)
    add_sheet(wb, "Residue_Values", residue_rows(data))
    add_sheet(wb, "Patches", patch_rows(data))
    add_sheet(wb, "Patch_Members", patch_member_rows(data))
    wb.save(path)


def rewrite_bfactor(raw_pdb, rows, label):
    lookup = {}
    for r in rows:
        key = (str(r.get("chain") or "").strip(), int(r.get("res_seq")), str(r.get("icode") or "").strip())
        v = r.get("propensity")
        lookup[key] = max(0.0, min(100.0, float(v or 0.0)))
    out = []
    for line in raw_pdb.replace("\r", "").split("\n"):
        if line.startswith(("ATOM  ", "HETATM")):
            padded = line.ljust(66)
            chain = padded[21:22].strip()
            try: resi = int(padded[22:26].strip())
            except ValueError: resi = 0
            icode = padded[26:27].strip()
            b = f"{lookup.get((chain, resi, icode), 0.0):6.2f}"
            line = padded[:60] + b + padded[66:]
        out.append(line)
    remarks = [
        "REMARK 900 INTERFACESCOUT B-FACTOR TRACK",
        f"REMARK 900 MAP: {label}",
        "REMARK 900 FIELD: canonical propensity; NORMALIZED RANGE 0-100",
        "REMARK 900 B-FACTOR VALUES ARE INTERFACESCOUT SCORES, NOT EXPERIMENTAL CRYSTALLOGRAPHIC B-FACTORS",
    ]
    return "\n".join(remarks + out)


def main():
    if EXAMPLES.exists():
        import shutil; shutil.rmtree(EXAMPLES)
    for case in CASES:
        print("Generating", case)
        raw = fetch_pdb(case["pdb"])
        data = analyze_all_maps(pdb_text=raw, pdb_id=case["pdb"], chain=case["chain"], pH=case["pH"], ionic_mM=150.0, temp_K=298.0)
        data["public_version"] = "1.0-publication"
        outdir = EXAMPLES / case["folder"]
        outdir.mkdir(parents=True, exist_ok=True)
        suffix = f"_pH{case['pH']}" if case["pdb"] == "2CBA" else ""
        write_excel(data, outdir / f"InterfaceScout_{case['pdb']}{suffix}_all_results.xlsx")
        m = data["maps"][case["map"]]
        pdb_out = rewrite_bfactor(raw, m.get("residues", []), m["label"])
        (outdir / f"InterfaceScout_{case['pdb']}{suffix}_{case['map']}_Bfactor.pdb").write_text(pdb_out, encoding="utf-8")

if __name__ == "__main__":
    main()
