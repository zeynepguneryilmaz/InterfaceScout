#!/usr/bin/env python3
"""Create reproducible OS-specific *source* archives for InterfaceScout releases.

These archives contain the application and launch scripts, not bundled Python
interpreters or compiled .exe/.app binaries. Dependencies are installed on
the destination computer using the per-platform instructions.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import zipfile
import hashlib

ROOT = Path(__file__).resolve().parents[1]
IGNORED_DIRS = {".git", ".venv", "venv", "__pycache__", ".pytest_cache",
                ".mypy_cache", ".idea", ".vscode", ".tox", "dist", "build"}
IGNORED_FILES = {".DS_Store", "startup_log.txt"}
COMMON = ("backend", "frontend", "docs", "examples", "validation",
          "tests", "tools", "README.md", "LICENSE", "CITATION.cff")
LAUNCHERS = {
    "Windows": ("run_local.bat", "start.bat"),
    "macOS": ("run_local.sh", "start.command"),
    "Linux": ("run_local.sh", "start.sh"),
}


def included_files(platform: str) -> list[Path]:
    selected = []
    for relative in COMMON + LAUNCHERS[platform]:
        source = ROOT / relative
        if not source.exists():
            raise FileNotFoundError(f"Required release item missing: {relative}")
        if source.is_file():
            selected.append(source)
        else:
            for path in source.rglob("*"):
                if any(part in IGNORED_DIRS for part in path.relative_to(ROOT).parts):
                    continue
                if path.is_file() and path.name not in IGNORED_FILES and path.suffix not in {".pyc", ".pyo"}:
                    selected.append(path)
    return sorted(set(selected))


def build(out_dir: Path, version: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    all_files = set()
    for platform in LAUNCHERS:
        archive = out_dir / f"InterfaceScout-{version}-{platform}-source.zip"
        with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED,
                             compresslevel=8) as output:
            for path in included_files(platform):
                all_files.add(path)
                name = "InterfaceScout/" + path.relative_to(ROOT).as_posix()
                output.write(path, name)
        print(f"{archive.name}: {archive.stat().st_size:,} bytes")
    bundle = out_dir / f"InterfaceScout_{version.replace('-', '_').replace('.', '_')}_publication_bundle.zip"
    with zipfile.ZipFile(bundle, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=8) as output:
        for path in sorted(all_files):
            output.write(path, "InterfaceScout/" + path.relative_to(ROOT).as_posix())
    print(f"{bundle.name}: {bundle.stat().st_size:,} bytes")
    checksums = []
    for path in sorted(out_dir.glob("*.zip")):
        checksums.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}")
    (out_dir / "SHA256SUMS.txt").write_text("\n".join(checksums) + "\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "dist")
    parser.add_argument("--version", required=True, help="Release tag, e.g. v1.1.0")
    args = parser.parse_args()
    build(args.out, args.version)
