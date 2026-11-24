"""Real browser full 11-map PDF/3D export of all 76 ISV candidates.

Produces genuine 3Dmol screenshots via InterfaceScout's own UI, not placeholders.
SHA-checks frozen assembly-1 mmCIF, uses production API and ZIP export, and
validates all 11 map PDFs, images, PDBs, Excel and real APBS residues.
No unrelated-protein consensus is constructed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from io import BytesIO
from pathlib import Path
import subprocess
import sys
import time
import traceback
from urllib.request import urlopen
import zipfile

from openpyxl import load_workbook
from PIL import Image, ImageChops, ImageStat
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "backend"), str(ROOT / "validation")]
from isv_full_80_qa import prepare, ROOTOUT as PREPARED_ROOT
from export_bundle import MAP_NAMES

ACTIVE = list(csv.DictReader(
    (ROOT / "validation" / "ISV_active_candidates_76_v3.csv").open(
        encoding="utf-8",newline=""
    )
))
assert len(ACTIVE) == 76 and len({x["proposed_pdb_id"] for x in ACTIVE}) == 76
assert not ({x["proposed_pdb_id"] for x in ACTIVE} & {"3GHG","1QJS","4DN0","6LQH"})
MANIFEST = [{"pdb_id": row["proposed_pdb_id"],"row": row} for row in ACTIVE]
OUT = ROOT / "isv76_full_3d_atlas"
OUT.mkdir(exist_ok=True)
URL = "http://127.0.0.1:8765"

def record(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def inspect_zip(path, pid, expected_pdb, expected_residues):
    prefix = f"structures/{pid}"
    with zipfile.ZipFile(path) as archive:
        if archive.testzip() is not None:
            raise AssertionError("ZIP CRC mismatch")
        names = set(archive.namelist())
        def read(relative):
            full = f"{prefix}/{relative}"
            if full not in names:
                raise AssertionError("Missing ZIP member: " + full)
            return archive.read(full)

        original = read("original.pdb")
        if hashlib.sha256(original).hexdigest() != hashlib.sha256(expected_pdb.encode()).hexdigest():
            raise AssertionError("Exported original PDB differs from analyzed prepared PDB")
        analysis = json.loads(read("analysis_complete.json"))
        expect_apbs = pid != "1GCT"
        electro = json.loads(read("electrostatics/APBS_PROPKA_residue_potentials.json")) if expect_apbs else {"status":"not_requested"}

        if analysis.get("n_maps") != 11 or list(analysis["map_order"]) != list(MAP_NAMES):
            raise AssertionError("Missing or reordered canonical maps")
        if analysis["electrostatics"] != electro or electro.get("status") != ("computed" if expect_apbs else "not_requested"):
            raise AssertionError("Unexpected APBS state in full browser export")
        if len(analysis["residue_profiles"]) != expected_residues:
            raise AssertionError("Residue count mismatch")
        if expect_apbs and electro.get("n_matched_residues",0) != expected_residues:
            raise AssertionError("PQR/PDB residue identity mismatch")
        exposed = [r for r in electro.get("residues",[]) if r["surface_exposed"]]
        sas = [r for r in exposed if r["surface_potential_kT_per_e"] is not None]
        if expect_apbs and len(sas) != electro["n_surface_potential_residues"]:
            raise AssertionError("Primary SAS count inconsistent")
        if expect_apbs and len(sas) < 0.90 * len(exposed):
            raise AssertionError("Incomplete primary SAS coverage")
        if expect_apbs and "surface_sampling" not in electro:
            raise AssertionError("Missing SAS method provenance")
        if expect_apbs and not any(r["res_name"] == "GLY" and
                   r["surface_potential_kT_per_e"] is not None and
                   r["sidechain_mean_potential_kT_per_e"] is None
                   for r in electro.get("residues",[])):
            raise AssertionError("No SAS data for exposed glycine")


        screenshots, pdfs, bfactor = [], [], []
        for i, key in enumerate(MAP_NAMES, 1):
            stem = f"{i:02d}_{key}"
            png = read(f"screenshots/{stem}.png")
            pdf = read(f"pdf_maps/{stem}.pdf")
            if not pdf.startswith(b"%PDF") or b"/Subtype /Image" not in pdf:
                raise AssertionError(f"PDF with actual viewer image missing: {key}")
            with Image.open(BytesIO(png)) as im:
                im.verify()
            with Image.open(BytesIO(png)) as im:
                rgb = im.convert("RGB")
                if rgb.width < 200 or rgb.height < 100 or max(ImageStat.Stat(rgb).stddev) < 3:
                    raise AssertionError(f"Blank/undersized actual 3D map screenshot: {key}")
            screenshots.append(png)
            pdfs.append(pdf)
            bfactor.append(read(f"pdb_bfactor/{stem}_Bfactor.pdb"))
            for k in (1, 3, 5):
                bfactor.append(read(f"pdb_bfactor/{stem}_Top{k}_Bfactor.pdb"))
        with Image.open(BytesIO(screenshots[0])) as a, Image.open(BytesIO(screenshots[2])) as b:
            if a.size != b.size or ImageChops.difference(a.convert("RGB"), b.convert("RGB")).getbbox() is None:
                raise AssertionError("11 map screenshots cannot be the same frozen image")
        atlas = read("Complete_Surface_Atlas.pdf")
        if not atlas.startswith(b"%PDF") or atlas != read("full_report.pdf"):
            raise AssertionError("Complete Surface Atlas missing or inconsistent")
        if b"/Subtype /Image" not in atlas:
            raise AssertionError("Surface Atlas lacks real screenshot images")
        for csv_path in ("Protein_Residues", "Map_Residue_Scores", "Map_Favorable_Residues",
                         "All_Patches", "Patch_Members", "Regional_Channels",
                         "Regional_Pair_Distances", "Regional_Method"):
            read(f"csv/{csv_path}.csv")
        if expect_apbs:
            read("csv/Electrostatics_Residues.csv")
            read("csv/Electrostatics_Settings.csv")
        read("combined/combined_maps.json")
        read("combined/selected_combination.png")
        wb = load_workbook(BytesIO(read("complete_results.xlsx")), read_only=True, data_only=True)
        expected_sheets = {"Map_Residue_Scores",
                           "Regional_Channels", "Regional_Pair_Distances"}
        if not expected_sheets.issubset(wb.sheetnames):
            raise AssertionError("Required Excel sheets missing: " + str(expected_sheets - set(wb.sheetnames)))
        if expect_apbs:
            first_line = read("csv/Electrostatics_Residues.csv").decode("utf-8").splitlines()[0]
            if "surface_potential_kT_per_e" not in first_line or "sidechain_mean_potential_kT_per_e" not in first_line:
                raise AssertionError("Both primary SAS and legacy atom values must export as CSV")
        if expect_apbs:
            electro_rows = wb["Electrostatics_Residues"]
            if electro_rows.max_row - 1 != expected_residues:
                raise AssertionError("Electrostatic Excel residue count mismatch")
        wb.close()
        return {
            "status": "PASS", "pdb_id": pid, "APBS_status": "computed" if expect_apbs else "unavailable_structural_gap", "zip_bytes": path.stat().st_size,
            "zip_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "n_maps": 11, "n_screenshots": len(screenshots), "n_map_pdfs": len(pdfs),
            "n_bfactor_pdbs": len(bfactor), "n_residues": expected_residues,
            "n_apbs_residues": electro.get("n_matched_residues",0),
            "n_pqr_atoms": electro.get("n_pqr_atoms",0),
            "n_surface_sas_potentials": len(sas),
            "n_zip_members": len(names), "electrostatics_status": electro["status"],
            "source_assembly_sha256_verified": True,
            "has_regional_colocation": True,
            "has_raw_APBS_DX_or_PQR": False,
        }


def run_one(browser, item):
    pid = item["pdb_id"].upper()
    details, pdbfile = prepare(item["row"])
    src = PREPARED_ROOT / f"{pid}-assembly1.cif"
    sha = hashlib.sha256(src.read_bytes()).hexdigest()
    if sha != details["assembly_cif_sha256"]:
        raise AssertionError(f"{pid}: source biological-assembly SHA mismatch")
    pdb = pdbfile.read_text(encoding="utf-8")
    page = browser.new_page(viewport={"width": 1200, "height": 850}, accept_downloads=True)
    js_errors = []
    page.on("pageerror", lambda error: js_errors.append(str(error)))
    try:
        page.goto(URL, wait_until="domcontentloaded", timeout=45000)
        page.wait_for_function("() => !!document.querySelector('#run') && !document.querySelector('#run').disabled", timeout=60000)
        page.set_input_files("#file", str(pdbfile))
        page.wait_for_function("() => typeof pdbText === 'string' && pdbText.includes('ATOM')")
        # The local prepared filename is *_prepared.pdb, but preserve the
        # frozen PDB accession as the exported structure directory name.
        page.evaluate("(pdbId) => { structureName = pdbId; }", pid)
        if pid != "1GCT":
            page.locator("#electrostatics").check()
        page.locator("#run").click()
        page.wait_for_function("(expected) => data && data.n_maps===11 && data.electrostatics?.status===expected",
                               arg=("computed" if pid != "1GCT" else "not_requested"),timeout=480000)
        if page.evaluate("() => document.querySelector('#error').textContent") and not page.evaluate("() => data"):
            raise AssertionError("Analysis returned an error")
        if pid != "1GCT" and not page.evaluate("() => data.electrostatics.n_matched_residues === data.residue_profiles.length"):
            raise AssertionError("APBS residue count mismatch")
        expected = page.evaluate("() => data.residue_profiles.length")
        # Same user-selectable hydrophobic + aromatic overlay, with full
        # regional co-localization and exact intersections in final exports.
        page.locator('[data-preset="hydrophobic_aromatic"]').click()
        page.wait_for_function(
            "() => decisionCombo && decisionCombo.selected_maps.length === 2 && selected.size===2",
            timeout=90000)
        if page.locator("#isFullExport").is_disabled():
            raise AssertionError("Full exporter disabled after analysis")
        filename = f"InterfaceScout_ISV76_{pid}_Final_Full_Atlas.zip"
        filepath = OUT / "zips" / filename
        filepath.parent.mkdir(parents=True, exist_ok=True)
        with page.expect_download(timeout=420000) as downloaded:
            page.locator("#isFullExport").click()
        download = downloaded.value
        download.save_as(filepath)
        report = inspect_zip(filepath, pid, pdb, expected)
        report["selected_regional_maps"] = ["hydrophobic", "pi_carbon"]
        report["prepared_pdb_sha256"] = hashlib.sha256(pdb.encode()).hexdigest()
        report["n_chains"] = details["n_chains"]
        report["candidate_category"] = item["row"]["category"]
        report["source_assembly_sha256"] = sha
        report["browser_js_errors"] = js_errors
        return report
    finally:
        page.close()


def run_shard(index, shards):
    server = subprocess.Popen(
        [sys.executable, "-m", "uvicorn", "app:app", "--app-dir", str(ROOT / "backend"),
         "--host", "127.0.0.1", "--port", "8765"],
        cwd=str(ROOT), stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT,
        env={**os.environ, "PYTHONUNBUFFERED": "1"}
    )
    records = []
    try:
        for attempt in range(120):
            if server.poll() is not None:
                raise RuntimeError("InterfaceScout backend exited")
            try:
                with urlopen(URL + "/health", timeout=1) as res:
                    if res.status == 200:
                        break
            except Exception:
                time.sleep(0.5)
        else:
            raise RuntimeError("Backend not responding")
        with sync_playwright() as pw:
            browser = pw.chromium.launch(
                headless=True,
                args=["--use-gl=angle", "--use-angle=swiftshader",
                      "--enable-unsafe-swiftshader", "--enable-webgl", "--disable-dev-shm-usage"]
            )
            try:
                for i, item in enumerate(MANIFEST):
                    if i % shards != index:
                        continue
                    pid = item["pdb_id"].upper()
                    try:
                        result = run_one(browser, item)
                    except Exception as exc:
                        result = {"pdb_id": pid, "status": "FAIL",
                                  "error_type": type(exc).__name__, "error": str(exc),
                                  "traceback": traceback.format_exc()[-5000:]}
                    records.append(result)
                    record(OUT / "records" / f"{pid}.json", result)
                    print("ISV76_ATLAS_CASE " + json.dumps({
                        "pdb_id": pid, "status": result["status"],
                        "bytes": result.get("zip_bytes"),
                        "error": result.get("error", "")[:500],
                    }), flush=True)
            finally:
                browser.close()
    finally:
        server.terminate()
        try:
            server.wait(timeout=7)
        except subprocess.TimeoutExpired:
            server.kill()
        record(OUT / f"shard_{index}.json", records)
    print("ISV76_ATLAS_SHARD " + json.dumps({
        "index": index, "attempted": len(records),
        "passed": sum(x["status"] == "PASS" for x in records),
        "failed": sum(x["status"] != "PASS" for x in records)
    }), flush=True)


def aggregate(shards):
    rows = []
    missing = []
    for n in range(shards):
        path = OUT / f"shard_{n}.json"
        if path.exists():
            rows += json.loads(path.read_text(encoding="utf-8"))
        else:
            missing.append(path.name)
    counts = {}
    for x in rows:
        counts[x["pdb_id"]] = counts.get(x["pdb_id"], 0) + 1
    absent = sorted({x["pdb_id"] for x in MANIFEST} - set(counts))
    failures = [x for x in rows if x["status"] != "PASS"]
    successful = [x for x in rows if x["status"] == "PASS"]
    summary = {
        "model": "InterfaceScout frozen post-TV48 2.5.0-dual-accessibility",
        "electrostatics": "real PDB2PQR/PROPKA/APBS; pH 7.4; 150 mM; 298 K",
        "requested": 76, "attempted": len(rows), "passed": len(successful),
        "n_real_APBS_atlases": sum(x.get("APBS_status") == "computed" for x in successful),
        "n_map_only_atlases": sum(x.get("APBS_status") == "unavailable_structural_gap" for x in successful),
        "failed": len(failures), "missing_shards": missing, "missing_ids": absent,
        "duplicate_ids": sorted(k for k,v in counts.items() if v > 1),
        "failed_structures": [{"pdb_id": x["pdb_id"], "error": x.get("error")}
                              for x in failures],
        "full_76_atlas_pass": len(successful) == 76 and not failures and not missing and not absent
                       and len(counts) == 76 and all(v == 1 for v in counts.values()),
        "real_viewer_pngs": 11 * len(successful),
        "real_map_pdf_files": 11 * len(successful),
        "BFactor_PDB_outputs": 44 * len(successful),
        "total_zip_bytes": sum(x["zip_bytes"] for x in successful),
        "scope": "76 candidate biological assemblies; 1GCT atlas maps only because APBS structural preparation fails; 8C50 electrostatic outlier is exploratory/quarantined.",
        "raw_solver_grids": "Not retained by the application's current APBS engine.",
    }
    record(OUT / "SUMMARY.json", summary)
    cols = ["pdb_id", "status", "APBS_status", "zip_bytes", "zip_sha256", "n_residues",
            "n_apbs_residues", "n_pqr_atoms", "n_screenshots", "n_map_pdfs",
            "n_bfactor_pdbs", "error_type", "error"]
    with (OUT / "ISV76_FULL_ATLAS_STATUS.csv").open("w", newline="", encoding="utf-8") as output:
        w = csv.DictWriter(output, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        by_id = {x["pdb_id"]: x for x in rows}
        for item in MANIFEST:
            w.writerow(by_id.get(item["pdb_id"], {"pdb_id": item["pdb_id"], "status": "NOT_ATTEMPTED"}))
    print("ISV76_ATLAS_SUMMARY " + json.dumps(summary), flush=True)
    if not summary["full_76_atlas_pass"]:
        raise SystemExit("Incomplete ISV76 full 3D/PDF atlas; inspect failed structure records.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard", type=int)
    parser.add_argument("--shards", type=int, default=16)
    parser.add_argument("--aggregate", action="store_true")
    args = parser.parse_args()
    if args.aggregate:
        aggregate(args.shards)
    elif args.shard is not None and 0 <= args.shard < args.shards:
        run_shard(args.shard, args.shards)
    else:
        parser.error("Provide --shard or --aggregate")


if __name__ == "__main__":
    main()
