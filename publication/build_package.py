"""Assemble the v1.1.0 RC archival bundle ONLY from passing CI evidence.

The source snapshot is exactly the tested Git checkout; external/untracked
material is never silently substituted for scientific source data.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import pathlib
import subprocess
import sys
import zipfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
SUPPORTED_SYSTEMS = ("Linux", "Windows", "macOS")
PUBLIC_VERSION = "1.1.0-rc1"
MODEL_VERSION = "2.4.0-selected-parameters"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_json(path: pathlib.Path) -> dict:
    if not path.is_file():
        raise ValueError(f"Missing evidence file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def audit_numerical(evidence: pathlib.Path) -> dict:
    folder = evidence / "numerical"
    result = load_json(folder / "verification.json")
    if result.get("status") != "PASS" or result.get("suite") != "full":
        raise ValueError("Full numerical verification is not PASS")
    cases = result.get("cases")
    if not isinstance(cases, list) or len(cases) != 21:
        raise ValueError("Expected 21 case-level numerical verification entries")
    if any(c.get("status") != "PASS" for c in cases):
        raise ValueError("At least one numerical case was not verified")
    if result.get("model_version") != MODEL_VERSION:
        raise ValueError("Numerical verification used the wrong model version")
    source_root = pathlib.Path(result.get("git_commit") or "")
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    if source_root.as_posix() != commit:
        raise ValueError("Numerical checks were executed on a different source revision")
    for item in cases:
        p = folder / (item["id"] + ".json")
        if not p.is_file():
            raise ValueError("Missing numerical per-condition output: " + item["id"])
        if sha256(p.read_bytes()) != item.get("output_sha256"):
            raise ValueError("Output SHA-256 mismatch: " + item["id"])
        if not item.get("oracle") or not item.get("bfactor"):
            raise ValueError("Numerical evidence omitted independent checks: " + item["id"])
    return {
        "status": "PASS", "n_conditions": len(cases), "n_maps": len(cases) * 11,
        "case_names": [c["id"] for c in cases],
        "verification_sha256": sha256((folder / "verification.json").read_bytes()),
        "source_commit": commit,
    }


def audit_browser(evidence: pathlib.Path) -> dict:
    verified = {}
    for os_name in SUPPORTED_SYSTEMS:
        base = evidence / "browser" / os_name
        rows = []
        for stem in ("single_pdb_browser_export", "multi_pdb_browser_export"):
            meta = load_json(base / (stem + ".json"))
            zipped = base / (stem + ".zip")
            if not zipped.is_file():
                raise ValueError("Missing real browser ZIP: " + str(zipped))
            if meta.get("passed_export_checks") is not True:
                raise ValueError("Browser checks did not pass: " + str(zipped))
            if sha256(zipped.read_bytes()) != meta.get("sha256"):
                raise ValueError("Browser ZIP SHA-256 mismatch: " + str(zipped))
            with zipfile.ZipFile(zipped) as archive:
                if archive.testzip() is not None or len(archive.namelist()) < 80:
                    raise ValueError("Damaged or incomplete browser ZIP: " + str(zipped))
            rows.append({
                "name": stem, "sha256": meta["sha256"], "bytes": zipped.stat().st_size
            })
        verified[os_name] = rows
    return {"status": "PASS", "systems": verified, "n_real_browser_exports": 6}


def tracked_files() -> list[pathlib.Path]:
    result = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT)
    paths = [ROOT / p.decode("utf-8") for p in result.split(b"\x00") if p]
    for p in paths:
        if not p.is_file() or p.is_symlink():
            raise ValueError("Invalid tracked source file: " + str(p))
    if not paths:
        raise ValueError("No Git-tracked source files were found")
    return paths


def build(evidence: pathlib.Path, output: pathlib.Path) -> pathlib.Path:
    numerical = audit_numerical(evidence)
    browser = audit_browser(evidence)
    paths = tracked_files()
    payload: dict[str, bytes] = {}
    for p in paths:
        relative = p.relative_to(ROOT).as_posix()
        payload["source/" + relative] = p.read_bytes()
    # Copy *only* the two verified evidence namespaces, no stray generated files.
    for folder in (evidence / "numerical", evidence / "browser"):
        for p in sorted(folder.rglob("*")):
            if p.is_file():
                relative = p.relative_to(evidence).as_posix()
                payload["validation-evidence/" + relative] = p.read_bytes()
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    manifest = {
        "tool": "InterfaceScout", "software_version": PUBLIC_VERSION,
        "scientific_model_version": MODEL_VERSION,
        "git_commit": revision,
        "built_at_utc": datetime.now(timezone.utc).isoformat(),
        "historical_version_preserved": "v1.0.0",
        "numerical_verification": numerical, "browser_verification": browser,
        "source_file_count": len(paths),
        "source_and_evidence_files": {
            key: {"sha256": sha256(data), "bytes": len(data)}
            for key, data in sorted(payload.items())
        },
        "scientific_limitation": (
            "Software correctness evidence and frozen literature benchmark are "
            "distinct; no adsorption-site predictive accuracy is established."
        ),
    }
    payload["PACKAGE_MANIFEST.json"] = (
        json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True)
        + "\n").encode("utf-8")
    payload["SHA256SUMS.txt"] = (
        "".join(f"{sha256(v)}  {k}\n" for k, v in sorted(payload.items()))
    ).encode("utf-8")
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED,
                         compresslevel=6, allowZip64=True) as archive:
        for key, data in sorted(payload.items()):
            archive.writestr(key, data)
    with zipfile.ZipFile(output) as archive:
        if archive.testzip() is not None:
            raise RuntimeError("Output publication ZIP is damaged")
        for key, record in manifest["source_and_evidence_files"].items():
            if sha256(archive.read(key)) != record["sha256"]:
                raise RuntimeError("Final ZIP payload checksum mismatch: " + key)
    digest = sha256(output.read_bytes())
    output.with_suffix(output.suffix + ".sha256").write_text(
        f"{digest}  {output.name}\n", encoding="utf-8")
    print(json.dumps({
        "status": "PASS", "file": str(output), "sha256": digest,
        "bytes": output.stat().st_size, "tracked_source_files": len(paths),
        "verified_numerical_cases": numerical["n_conditions"],
        "verified_maps": numerical["n_maps"],
        "real_browser_zips": browser["n_real_browser_exports"],
    }, indent=2))
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-dir", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    a = parser.parse_args()
    try:
        build(a.evidence_dir, a.output)
    except (ValueError, OSError, RuntimeError, KeyError,
            subprocess.CalledProcessError) as exc:
        print("PUBLICATION BUNDLE REJECTED: " + str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
