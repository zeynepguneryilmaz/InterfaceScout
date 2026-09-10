#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
BACKEND="$HERE/backend"
if [ ! -f "$BACKEND/app.py" ] || [ ! -f "$HERE/frontend/index.html" ]; then
  osascript -e 'display dialog "InterfaceScout files are incomplete." buttons {"OK"} default button "OK" with title "InterfaceScout"' >/dev/null 2>&1 || true
  exit 1
fi
if [ "${1:-}" = "--check" ]; then
  bash -n "$0"
  echo "InterfaceScout macOS launcher check passed"
  exit 0
fi
if [ ! -f "$BACKEND/.venv/bin/activate" ]; then
  osascript -e 'display dialog "First-time setup is required. Run: bash run_local.sh" buttons {"OK"} default button "OK" with title "InterfaceScout"' >/dev/null 2>&1 || true
  exit 1
fi
if curl -fsS "http://127.0.0.1:8000/health" >/dev/null 2>&1; then
  open "http://localhost:8000" >/dev/null 2>&1 || true
  exit 0
fi
if command -v lsof >/dev/null 2>&1 && lsof -ti tcp:8000 >/dev/null 2>&1; then
  osascript -e 'display dialog "Port 8000 is already in use." buttons {"OK"} default button "OK" with title "InterfaceScout"' >/dev/null 2>&1 || true
  exit 1
fi
cd "$BACKEND"
# shellcheck disable=SC1091
source .venv/bin/activate
nohup python -m uvicorn app:app --host 127.0.0.1 --port 8000 >"${TMPDIR:-/tmp}/interfacescout.log" 2>&1 &
for _ in 1 2 3 4 5 6 7 8 9 10; do
  if curl -fsS "http://127.0.0.1:8000/health" >/dev/null 2>&1; then
    open "http://localhost:8000" >/dev/null 2>&1 || true
    exit 0
  fi
  sleep 1
done
osascript -e 'display dialog "InterfaceScout did not start. Check the InterfaceScout log in the system temporary folder." buttons {"OK"} default button "OK" with title "InterfaceScout"' >/dev/null 2>&1 || true
exit 1
