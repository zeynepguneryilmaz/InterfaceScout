#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
BACKEND="$HERE/backend"
if [ ! -f "$BACKEND/app.py" ] || [ ! -f "$HERE/frontend/index.html" ]; then
  echo "ERROR: InterfaceScout files are incomplete." >&2
  exit 1
fi
if [ "${1:-}" = "--check" ]; then
  bash -n "$0"
  echo "InterfaceScout Linux launcher check passed"
  exit 0
fi
if [ ! -f "$BACKEND/.venv/bin/activate" ]; then
  echo "First-time setup is required. Run: bash run_local.sh" >&2
  exit 1
fi
if command -v curl >/dev/null 2>&1 && curl -fsS "http://127.0.0.1:8000/health" >/dev/null 2>&1; then
  command -v xdg-open >/dev/null 2>&1 && xdg-open "http://localhost:8000" >/dev/null 2>&1 || true
  exit 0
fi
if command -v lsof >/dev/null 2>&1 && lsof -ti tcp:8000 >/dev/null 2>&1; then
  echo "ERROR: Port 8000 is already in use." >&2
  exit 1
fi
cd "$BACKEND"
# shellcheck disable=SC1091
source .venv/bin/activate
nohup python -m uvicorn app:app --host 127.0.0.1 --port 8000 >"${TMPDIR:-/tmp}/interfacescout.log" 2>&1 &
for _ in 1 2 3 4 5 6 7 8 9 10; do
  if command -v curl >/dev/null 2>&1 && curl -fsS "http://127.0.0.1:8000/health" >/dev/null 2>&1; then
    command -v xdg-open >/dev/null 2>&1 && xdg-open "http://localhost:8000" >/dev/null 2>&1 || true
    exit 0
  fi
  sleep 1
done
echo "ERROR: InterfaceScout did not start. Check ${TMPDIR:-/tmp}/interfacescout.log" >&2
exit 1
