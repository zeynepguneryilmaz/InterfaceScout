"""First-run launcher compatibility on Linux, Windows and macOS.

Performs actual venv installation through run_local.{bat,sh}, starts the server
using the shipped launcher and tests its HTTP APIs. This is not merely --check.
"""
from __future__ import annotations

import json
import platform
from pathlib import Path
import subprocess
import sys
import time
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
PDB = ROOT / "examples" / "1CRN" / "InterfaceScout_1CRN_anionic_Bfactor.pdb"


def call(url: str, data: dict | None = None):
    if data is not None:
        raw = json.dumps(data).encode("utf-8")
        req = Request(url, data=raw, headers={"Content-Type": "application/json"})
    else:
        req = Request(url)
    with urlopen(req, timeout=120) as r:
        assert r.status == 200, (url, r.status)
        return r.read()


def main():
    if sys.platform == "win32":
        cmd = ["cmd", "/c", str(ROOT / "run_local.bat")]
    else:
        cmd = ["bash", str(ROOT / "run_local.sh")]
    print("Running full first-time launcher:", cmd, platform.platform(), flush=True)
    t0 = time.perf_counter()
    p = subprocess.run(cmd, cwd=str(ROOT), capture_output=True,
                       text=True, timeout=360, errors="replace")
    print("Launcher stdout (last 2500 chars):", p.stdout[-2500:], flush=True)
    if p.returncode:
        print("Launcher stderr (last 2500 chars):", p.stderr[-2500:], flush=True)
    assert p.returncode == 0, f"First-run launcher exited with code {p.returncode}"

    base = "http://127.0.0.1:8000"
    for _ in range(35):
        try:
            health = json.loads(call(base + "/health"))
            if health.get("status") == "ok":
                break
        except Exception:
            time.sleep(.5)
    else:
        raise AssertionError("Launcher returned but local API is not healthy")
    home = call(base + "/")
    assert b"InterfaceScout" in home, "Frontend not served from launcher"
    channels = json.loads(call(base + "/chemistries"))
    assert len(channels) == 11, f"Wrong channel count: {len(channels)}"
    goals = json.loads(call(base + "/research-goals"))
    assert len(goals) >= 3
    result = json.loads(call(base + "/analyze", {
        "pdb_text": PDB.read_text("utf-8"),
        "chain": "A", "pH": 7.4, "ionic_mM": 150.0, "temp_K": 298.0,
    }))
    assert result["n_maps"] == 11
    print(
        "FIRST-INSTALL AND LAUNCH PASS:",
        platform.platform(),
        "11 maps analyzed from the actual application server,",
        f"elapsed {time.perf_counter() - t0:.1f}s",
        flush=True
    )


if __name__ == "__main__":
    main()
