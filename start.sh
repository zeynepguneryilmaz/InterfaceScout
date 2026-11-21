#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
BACKEND="$HERE/backend"
if [ "${1:-}" = "--check" ]; then bash -n "$0"; echo "InterfaceScout Linux launcher check passed"; exit 0; fi
if [ ! -f "$BACKEND/.venv/bin/activate" ]; then echo "ERROR: First-time setup is required. Run: bash run_local.sh" >&2; exit 1; fi
if [ ! -f "$HERE/frontend/index.html" ]; then echo "ERROR: InterfaceScout frontend is incomplete." >&2; exit 1; fi
cd "$BACKEND"
source .venv/bin/activate

healthcheck() {
python - <<'PY' >/dev/null 2>&1
import urllib.request
try:
    r=urllib.request.urlopen("http://127.0.0.1:8000/health",timeout=1.5)
    raise SystemExit(0 if r.status==200 else 1)
except Exception:
    raise SystemExit(1)
PY
}

if healthcheck; then command -v xdg-open >/dev/null 2>&1 && xdg-open "http://localhost:8000" >/dev/null 2>&1 || true; exit 0; fi
nohup python -m uvicorn app:app --host 127.0.0.1 --port 8000 >"${TMPDIR:-/tmp}/interfacescout.log" 2>&1 &
for _ in $(seq 1 15); do
  if healthcheck; then command -v xdg-open >/dev/null 2>&1 && xdg-open "http://localhost:8000" >/dev/null 2>&1 || true; exit 0; fi
  sleep 1
done
echo "ERROR: InterfaceScout did not start. Check ${TMPDIR:-/tmp}/interfacescout.log" >&2
exit 1
