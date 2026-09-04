#!/usr/bin/env bash
# ============================================================
# InterfaceScout - Linux daily launcher
# ============================================================

HERE="$(cd "$(dirname "$0")" && pwd)"
BACKEND="$HERE/backend"
FRONTEND="$HERE/frontend"

if [ ! -f "$BACKEND/v52_run.py" ] || [ ! -f "$BACKEND/v52_app.py" ] || [ ! -f "$BACKEND/main.py" ] || [ ! -f "$FRONTEND/index.html" ]; then
  echo "InterfaceScout files are incomplete. Keep start.sh next to backend/ and frontend/." >&2
  exit 1
fi

if [ ! -f "$BACKEND/.venv/bin/activate" ]; then
  echo "First-time setup is required. Run run_local.sh once." >&2
  exit 1
fi

if command -v curl >/dev/null 2>&1 && curl -fsS "http://127.0.0.1:8000/health" >/dev/null 2>&1; then
  command -v xdg-open >/dev/null 2>&1 && xdg-open "http://localhost:8000" >/dev/null 2>&1 || true
  exit 0
fi

if command -v lsof >/dev/null 2>&1 && lsof -ti tcp:8000 >/dev/null 2>&1; then
  echo "Port 8000 is already in use by another process." >&2
  exit 1
fi

cd "$BACKEND" || exit 1
# shellcheck disable=SC1091
source .venv/bin/activate

nohup python v52_run.py >/tmp/interfacescout.log 2>&1 &
sleep 2

command -v xdg-open >/dev/null 2>&1 && xdg-open "http://localhost:8000" >/dev/null 2>&1 || true
exit 0
