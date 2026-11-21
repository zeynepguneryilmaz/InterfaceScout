"""Smoke-test a downloaded InterfaceScout release on the current OS."""
from __future__ import annotations
import json
from pathlib import Path
import subprocess
import sys
import time
import urllib.request
import urllib.error

root = Path(sys.argv[1]).resolve()
fixture = Path(sys.argv[2]).resolve()
assert (root / "backend/app.py").is_file()
assert (root / "frontend/index.html").is_file()
assert fixture.is_file()
if sys.platform == "win32":
    checks = [["cmd", "/c", "call", "run_local.bat", "--check"], ["cmd", "/c", "call", "start.bat", "--check"]]
elif sys.platform == "darwin":
    checks = [["bash", "run_local.sh", "--check"], ["bash", "start.command", "--check"]]
else:
    checks = [["bash", "run_local.sh", "--check"], ["bash", "start.sh", "--check"]]
for command in checks:
    subprocess.run(command, cwd=root, check=True, timeout=45)

proc = subprocess.Popen(
    [sys.executable, "-m", "uvicorn", "app:app", "--host", "127.0.0.1", "--port", "18761"],
    cwd=root / "backend", stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
)
base = "http://127.0.0.1:18761"
try:
    health = None
    for i in range(60):
        if proc.poll() is not None:
            raise RuntimeError(f"Backend exited: {proc.returncode}")
        try:
            with urllib.request.urlopen(base + "/health", timeout=2) as response:
                health = json.load(response)
            break
        except (urllib.error.URLError, TimeoutError, ConnectionError):
            time.sleep(1)
    assert health and health.get("status") == "ok", "Health endpoint failed"
    print("Release backend HEALTH PASS:", health, flush=True)

    request = urllib.request.Request(
        base + "/analyze",
        data=json.dumps({"pdb_text": fixture.read_text(encoding="utf-8", errors="replace"),
                         "chain": "A", "pH": 7.4}).encode("utf-8"),
        headers={"Content-Type": "application/json"}, method="POST",
    )
    with urllib.request.urlopen(request, timeout=240) as response:
        result = json.load(response)
    assert result.get("n_maps") == 11
    assert len(result.get("maps", {})) == 11
    assert all("patches" in m and "residues" in m for m in result["maps"].values())
    print("Release ANALYZE PASS: 11 chemistry maps", flush=True)
finally:
    proc.terminate()
    try:
        proc.communicate(timeout=10)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.communicate(timeout=5)
