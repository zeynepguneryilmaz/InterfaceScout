"""Corrected biological-unit InterfaceScout batch runner.

Uses the public app unchanged: real RCSB coordinates -> prepared PDB -> real
Chromium/3Dmol UI -> canonical 11-map ZIP export. Experimental outcomes are
not supplied to the engine. A failed/ambiguous case is explicitly recorded.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
from io import BytesIO
import itertools
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import traceback
import urllib.request
import zipfile

import gemmi
from PIL import Image
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))
from export_bundle import MAP_NAMES
from research_workflows import combine_maps, compare_results

GROUPS = [
    ["1A5R", "1AKI", "1AML", "1AO6", "1BEB"],
    ["1F6R", "1F6S", "1FNF", "1H3X", "1HRC"],
    ["1HZH", "1IGT", "1JNJ", "1MBN", "1TGU"],
    ["1UBQ", "2LZM", "2CBA", "2GB1", "2HAV", "3QYT"],
    ["2HHB", "2OED", "2PTN", "2VH7", "4EUL"],
    ["3ECA", "3GHG", "3STN", "4F5S", "5J17"],
    ["6O5E", "7KWI", "7RSA", "8ONG"],
]
assert len([p for group in GROUPS for p in group]) == 35
assert len(set(p for group in GROUPS for p in group)) == 35

# The chain selections are explicitly pre-outcome. An assembly choice does not
# assert that all experimentally unobserved subunits were structurally resolved.
ASSEMBLY = {"1BEB", "1H3X", "1HZH", "1IGT", "1TGU", "2HHB", "3ECA", "3GHG"}
SELECTED_CHAIN = {
    "1F6S": "B",  # separate Ca-bound alpha-lactalbumin state
}
DEFAULT_CHAIN = "A"
SCOPE_NOTES = {
    "1A5R": "NMR model 1 only; not an ensemble population",
    "1AML": "Amyloid-beta fragment; not folded whole precursor",
    "1FNF": "FNIII7-10 fragment; not full-length fibronectin",
    "1H3X": "Fc fragment; not whole IgG",
    "1HZH": "Intact antibody; source clone may not match a literature assay",
    "1IGT": "Mouse mAb231; clone-specific experimental transfer unsupported",
    "1HRC": "Covalent heme is stripped from the protein-only chemistry model",
    "1MBN": "Heme is stripped; no ligand-aware electrostatics claim",
    "1TGU": "Catalase heme tetramer; heme groups stripped",
    "2CBA": "Zn environment omitted by the protein-only chemistry kernel",
    "2HAV": "Apo transferrin; pair with iron-bound 3QYT",
    "3QYT": "Diferric transferrin; pair with apo 2HAV; iron removed from model",
    "1F6R": "Apo alpha-lactalbumin; pair with Ca-containing 1F6S",
    "1F6S": "Ca-bound alpha-lactalbumin; Ca removed; incomplete GLU121 sidechain",
    "2HHB": "Deoxy Hb tetramer; heme removed from chemistry model",
    "3GHG": "Only experimentally resolved fibrinogen material; NOT entire fibrinogen",
    "4EUL": "EGFP engineered protein; engineered chromophore chemistry not modeled",
    "6O5E": "Vitronectin domain; not full-length vitronectin",
    "7KWI": "NMR model 1; R2ab construct only",
    "8ONG": "NMR model 1; EMAP II construct only",
}
CANON = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
SERVER = "http://127.0.0.1:8765"
OUT = Path("curated_biological_units_outputs")

def sha256(data):
    return hashlib.sha256(data).hexdigest()

def fetch_bytes(url, attempts=3):
    last = None
    for i in range(attempts):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "InterfaceScout-curation/1.0 (scientific-reproducibility)"})
            with urllib.request.urlopen(req, timeout=100) as response:
                blob = response.read()
            if len(blob) < 1000:
                raise RuntimeError("RCSB returned suspiciously short content")
            return blob
        except Exception as exc:
            last = exc
            time.sleep(2 ** i)
    raise RuntimeError(f"RCSB retrieval failed {url}: {last}")

def prepare_case(pid, folder):
    folder.mkdir(parents=True, exist_ok=True)
    use_assembly = pid in ASSEMBLY
    source_url = ("https://files.rcsb.org/download/" +
                  (f"{pid}-assembly1.cif" if use_assembly else f"{pid}.cif"))
    source = fetch_bytes(source_url)
    src_path = folder / ("rcsb_assembly1.cif" if use_assembly else "rcsb_entry.cif")
    src_path.write_bytes(source)
    # Store the entry coordinates as well; the assembly coordinates can contain
    # generated symmetry mates and have a different source hash.
    if use_assembly:
        try:
            entry_blob = fetch_bytes(f"https://files.rcsb.org/download/{pid}.cif")
            (folder / "rcsb_entry.cif").write_bytes(entry_blob)
        except Exception as exc:
            (folder / "entry_download_warning.txt").write_text(str(exc))
    st = gemmi.read_structure(str(src_path))
    n_models = len(st)
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    chain_want = SELECTED_CHAIN.get(pid, DEFAULT_CHAIN)
    chosen = list(st[0]) if use_assembly else [chain for chain in st[0] if chain.name == chain_want]
    if not chosen:
        raise ValueError(f"Required chain {chain_want} absent from {pid}; "
                         f"observed {[x.name for x in st[0]]}. No fallback allowed.")
    clean = gemmi.Structure()
    clean.name = pid
    clean.cell = st.cell
    clean.spacegroup_hm = st.spacegroup_hm
    model = gemmi.Model(1)
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    mapping = []
    count_ca_only = 0
    for chain in chosen:
        residues = [r for r in chain.get_polymer() if r.name in CANON]
        if not residues:
            continue
        if len(mapping) >= len(alphabet):
            raise ValueError("PDB format chain-ID limit exceeded; hold for mmCIF-aware export")
        out_chain = gemmi.Chain(alphabet[len(mapping)])
        for residue in residues:
            names = {atom.name for atom in residue}
            if "CA" in names and not {"N", "C", "O"} & names:
                count_ca_only += 1
            out_chain.add_residue(residue.clone())
        model.add_chain(out_chain)
        mapping.append({"auth_chain": chain.name, "export_chain": out_chain.name,
                        "modeled_standard_residues": len(residues)})
    if not mapping:
        raise ValueError("No standard protein residues in chosen biological unit")
    clean.add_model(model)
    pdb = clean.make_pdb_string()
    protein_file = folder / "prepared_protein_only.pdb"
    protein_file.write_text(pdb, encoding="utf-8")
    nres = sum(m["modeled_standard_residues"] for m in mapping)
    if nres >= 20 and count_ca_only / nres >= 0.80:
        raise ValueError("C-alpha-only structure; unsuitable for physicochemical side-chain maps")
    provenance = {
        "pdb_id": pid,
        "source_url": source_url,
        "source_sha256": sha256(source),
        "protein_only_sha256": sha256(pdb.encode()),
        "selection": "wwPDB biological assembly 1" if use_assembly else f"entry model 1 auth chain {chain_want}",
        "source_model_count": n_models,
        "model_used": 1,
        "n_protein_chains": len(mapping),
        "chain_mapping": mapping,
        "n_modelled_residues": nres,
        "ca_only_residues": count_ca_only,
        "reference_condition": {"pH": 7.4, "ionic_mM": 150, "temperature_K": 298},
        "modeling_limitations": [
            "Cofactors, ions, nonstandard residues, hydrogens and glycans excluded from protein-only mapping",
            "pH may affect chemical channels; ionic strength and temperature are documented metadata only",
            "No experimental adsorption-site labels supplied to engine",
            "Structural maps do not predict affinity, adsorption pose or conformer probabilities",
        ],
        "scope_note": SCOPE_NOTES.get(pid, "Check original study construct/state before any experimental interpretation"),
    }
    (folder / "source_and_selection.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False), encoding="utf-8")
    return provenance, protein_file

def audit_zip(path, pid):
    with zipfile.ZipFile(path) as z:
        bad = z.testzip()
        if bad is not None:
            raise AssertionError(f"Invalid archive member: {bad}")
        names = z.namelist()
        roots = sorted({n.split("/")[1] for n in names if n.startswith("structures/") and len(n.split("/")) > 2})
        if len(roots) != 1:
            raise AssertionError("Expected one correctly selected protein structure per ZIP")
        prefix = "structures/" + roots[0] + "/"
        def ext(prefix_sub, extension):
            return sorted(n for n in names if n.startswith(prefix + prefix_sub) and n.endswith(extension))
        pngs = ext("screenshots/", ".png")
        primary_pngs = [x for x in pngs if not x.endswith("_reverse.png")]
        pdfs = ext("pdf_maps/", ".pdf")
        pdbs = ext("pdb_bfactor/", ".pdb")
        if len(primary_pngs) != 11 or len(pdfs) != 11 or len(pdbs) != 44:
            raise AssertionError(f"Insufficient outputs: PNG={len(primary_pngs)}, PDF={len(pdfs)}, PDB={len(pdbs)}")
        for req in ["analysis_complete.json", "complete_results.xlsx", "full_report.pdf", "original.pdb"]:
            if prefix + req not in names:
                raise AssertionError(f"Missing {req}")
        a = json.loads(z.read(prefix + "analysis_complete.json"))
        if a["n_maps"] != 11 or list(a["map_order"]) != MAP_NAMES:
            raise AssertionError("11 frozen canonical map channels not present")
        original = z.read(prefix + "original.pdb")
        if a["input"].get("source_pdb_sha256") != sha256(original):
            raise AssertionError("Input SHA mismatch in engine's primary export")
        for fn in pdfs:
            blob = z.read(fn)
            if not blob.startswith(b"%PDF-") or b"/Subtype /Image" not in blob:
                raise AssertionError(f"Map PDF lacks embedded image: {fn}")
        for fn in primary_pngs:
            with Image.open(BytesIO(z.read(fn))) as im:
                im.verify()
        if not ext("csv/", ".csv"):
            raise AssertionError("No per-map CSV exports")
        return a, {"png": len(primary_pngs), "reverse_png": len(pngs)-len(primary_pngs),
                   "pdf_maps": len(pdfs), "bfactor_pdb": len(pdbs),
                   "csv": len(ext("csv/", ".csv")),
                   "xlsx": 1, "json": 1,
                   "zip_members": len(names), "zip_sha256": sha256(path.read_bytes())}

def save_55_pairwise(analysis, folder):
    rows = []
    for c1,c2 in itertools.combinations(MAP_NAMES,2):
        result = combine_maps(analysis, [c1,c2], max_regions=None)
        rows.append({"map_1": c1, "map_2": c2,
                     "n_regions": result.get("regions_total",0),
                     "n_exact_overlap": result.get("n_exact_overlap",0),
                     "n_colocated_regions": sum(bool(x.get("all_maps_colocated")) for x in result.get("regions",[])),
                     "method": "separate unweighted geometric overlay; no adsorption binding prediction"})
    path = folder/"55_pairwise_map_colocations.csv"
    with path.open("w",newline="",encoding="utf-8-sig") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)

def start_server():
    proc=subprocess.Popen([sys.executable,"-m","uvicorn","app:app","--app-dir",str(ROOT/"backend"),
                           "--host","127.0.0.1","--port","8765"],
                          stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT,cwd=ROOT)
    for _ in range(100):
        if proc.poll() is not None:
            raise RuntimeError("InterfaceScout server terminated")
        try:
            with urllib.request.urlopen(SERVER+"/health",timeout=2) as r:
                if r.status == 200:
                    return proc
        except Exception:
            time.sleep(1)
    raise TimeoutError("InterfaceScout server health check failed")

def browser_case(browser, pid, prepared, folder):
    page=browser.new_page(viewport={"width": 1200,"height": 850},accept_downloads=True)
    page.set_default_timeout(120000)
    errors=[]
    page.on("pageerror",lambda err: errors.append(str(err)))
    try:
        page.goto(SERVER,wait_until="domcontentloaded",timeout=90000)
        page.locator("#run").wait_for(state="visible")
        page.set_input_files("#file",str(prepared.resolve()))
        page.wait_for_function("() => typeof pdbText === 'string' && pdbText.includes('ATOM')",timeout=60000)
        page.locator("#chain").fill("")
        page.locator("#run").click()
        page.wait_for_function("() => typeof data!=='undefined' && data!==null && data.n_maps===11",timeout=480000)
        page.locator('[data-preset="hydrophobic_aromatic"]').click()
        page.wait_for_function("() => selected.has('hydrophobic') && selected.has('pi_carbon') && selected.size === 2",timeout=90000)
        # Invoke the same complete export as a user, including the real WebGL
        # 3Dmol images generated by ensembleCaptureEleven().
        path=folder/(pid+"_Complete_Results.zip")
        with page.expect_download(timeout=600000) as info:
            page.locator("#isFullExport").click()
        info.value.save_as(str(path))
        if not path.exists() or path.stat().st_size < 30000:
            raise RuntimeError("Full browser export not downloaded")
        (folder/"browser_page_errors.json").write_text(json.dumps(errors,indent=2))
        return path
    except Exception:
        (folder/"browser_page_errors.json").write_text(json.dumps(errors,indent=2))
        try: (folder/"ui_error.txt").write_text(page.locator("#isError").inner_text()[:5000])
        except Exception: pass
        raise
    finally:
        page.close()

def compare_pair(left, right, output_dir):
    folder = output_dir / f"comparison_{left}_{right}"
    folder.mkdir(exist_ok=True)
    def analysis(pid):
        root = output_dir/pid
        arc = root/(pid+"_Complete_Results.zip")
        with zipfile.ZipFile(arc) as z:
            key=next(x for x in z.namelist() if x.endswith("/analysis_complete.json"))
            return json.loads(z.read(key))
    try:
        a, b=analysis(left),analysis(right)
        comp=compare_results(a,b)
        (folder/"comparison.json").write_text(json.dumps(comp,indent=2,ensure_ascii=False))
        return {"pair":f"{left}/{right}","status":"PASS_MATCHED_RESIDUE_COMPARISON","output":"comparison.json"}
    except Exception as exc:
        (folder/"comparison_failure.txt").write_text(str(exc))
        return {"pair":f"{left}/{right}","status":"UNMATCHED_OR_FAILED","reason":str(exc)}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--shard",type=int,required=True)
    args=parser.parse_args()
    if not 0<=args.shard<len(GROUPS):
        raise SystemExit("shard must be 0-6")
    group=GROUPS[args.shard]
    OUT.mkdir(exist_ok=True)
    root=OUT/f"shard_{args.shard}"
    root.mkdir(exist_ok=True)
    records=[]
    server=start_server()
    try:
        with sync_playwright() as p:
            browser=p.chromium.launch(headless=True,args=[
                "--use-gl=angle","--use-angle=swiftshader","--enable-unsafe-swiftshader",
                "--enable-webgl","--disable-dev-shm-usage"])
            try:
                for pid in group:
                    path=root/pid
                    row={"pdb_id":pid,"selection_type":"assembly1" if pid in ASSEMBLY else "single_auth_chain",
                         "status":"NOT_STARTED"}
                    try:
                        source, prepared=prepare_case(pid,path)
                        row.update({"selected_chains":source["chain_mapping"],
                                    "n_residues":source["n_modelled_residues"],
                                    "coordinate_sha256":source["protein_only_sha256"]})
                        bundle=browser_case(browser,pid,prepared,path)
                        analysis,counts=audit_zip(bundle,pid)
                        pairs=save_55_pairwise(analysis,path)
                        row.update({"status":"PASS_FULL_REAL_BROWSER_EXPORT",
                                    "n_maps":11,"n_pairwise_colocations":pairs,
                                    "n_surface":analysis["diagnostics"]["n_surface_residues"],
                                    "output_counts":counts})
                    except Exception as exc:
                        row.update({"status":"FAILED_OR_HELD","error":f"{type(exc).__name__}: {exc}"})
                        (path/"traceback.txt").write_text(traceback.format_exc(),encoding="utf-8")
                    records.append(row)
                    (root/"case_status.json").write_text(json.dumps(records,indent=2,ensure_ascii=False),encoding="utf-8")
                    print("CURATED_MAP_RUN",pid,row["status"],row.get("error","")[:300],flush=True)
            finally:
                browser.close()
    finally:
        server.terminate()
        try: server.wait(timeout=6)
        except subprocess.TimeoutExpired: server.kill()
    comparisons=[]
    for left,right in [("1F6R","1F6S"),("2HAV","3QYT")]:
        if left in group and right in group and all(
            next((r for r in records if r["pdb_id"]==pid),{}).get("status")=="PASS_FULL_REAL_BROWSER_EXPORT"
            for pid in (left,right)):
            comparisons.append(compare_pair(left,right,root))
    (root/"comparison_status.json").write_text(json.dumps(comparisons,indent=2))
    with (root/"case_status.csv").open("w",newline="",encoding="utf-8-sig") as f:
        cols=["pdb_id","status","selection_type","n_residues","n_surface","n_maps",
              "n_pairwise_colocations","coordinate_sha256","selected_chains","output_counts","error"]
        writer=csv.DictWriter(f,fieldnames=cols)
        writer.writeheader()
        for row in records:
            writer.writerow({k:json.dumps(row.get(k),ensure_ascii=False) if isinstance(row.get(k),(list,dict)) else row.get(k,"") for k in cols})
    (root/"README.md").write_text(
        "InterfaceScout exact-engine real-browser atlas results. Corrected biological units are selected "
        "without experimental outcome inspection. Every passing PDB contains the original RCSB mmCIF, "
        "prepared protein-only PDB, source and selection JSON, 55 pairwise co-location table and a "
        "real 3Dmol browser-generated complete ZIP (11 map PDFs, PNGs, B-factor PDBs, CSV, XLSX, JSON). "
        "Any failed PDB has an explicit error and is NOT a passing analysis. "
        "Bioassembly annotations, construct/ligand state, omitted components and missing residues "
        "remain scientific limitations. No site/affinity prediction accuracy implied.",
        encoding="utf-8")
    # Deliberately fail the job if any requested case was omitted; artifact
    # upload uses if: always() so partial output remains accessible for audit.
    if len(records)!=len(group) or any(r["status"]!="PASS_FULL_REAL_BROWSER_EXPORT" for r in records):
        raise SystemExit("One or more cases have no fully validated export; inspect case_status.csv")
    print("SHARD SUCCESS",args.shard,len(records),flush=True)

if __name__=="__main__":
    main()
