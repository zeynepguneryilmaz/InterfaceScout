"""Check pinned InterfaceScout requirements without contacting PyPI.

Returns 0 when all pinned distributions are installed at the specified versions;
returns 1 with human-readable diagnostics otherwise. Uses only Python stdlib.
"""
from __future__ import annotations
from importlib import metadata
from pathlib import Path
import sys

REQUIREMENTS = Path(__file__).with_name("requirements.txt")


def main() -> int:
    problems = []
    for raw in REQUIREMENTS.read_text(encoding="utf-8").splitlines():
        spec = raw.partition("#")[0].strip()
        if not spec:
            continue
        if "==" not in spec:
            problems.append("Unsupported requirement (expected name==version): " + spec)
            continue
        name, required = (part.strip() for part in spec.split("==", 1))
        try:
            installed = metadata.version(name)
        except metadata.PackageNotFoundError:
            problems.append(f"{name}=={required} is not installed")
            continue
        if installed != required:
            problems.append(f"{name}: installed {installed}, required {required}")
    if problems:
        print("Pinned dependencies are not ready:", flush=True)
        for problem in problems:
            print("  - " + problem, flush=True)
        return 1
    print("All pinned InterfaceScout dependencies are already installed.", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
