#!/usr/bin/env bash
# ============================================================
# InterfaceScout - Linux daily launcher
# ============================================================

HERE="$(cd "$(dirname "$0")" && pwd)"
BACKEND="$HERE/backend"
FRONTEND="$HERE/frontend"

if [ ! -f "$BACKEND/main.py" ] || [ ! -f "$FRONTEND/index.html" ]; then
  echo "InterfaceScout files are incomplete. Keep start.sh next to backend/ and frontend/." >&2
  exit 1
fi

cd "$BACKEND" || exit 1

if [ ! -f ".venv/bin/activate" ]; then
  if command -v zenity >/dev/null 2>&1; then
    zenity --error --title="InterfaceScout" \
      --text="First-time setup is required. Run run_local.sh once from the InterfaceScout folder." 2>/dev/null
  else
    echo "First-time setup is required. Run run_local.sh once." >&2
  fi
  exit 1
fi

# shellcheck disable=SC1091
source .venv/bin/activate

if command -v curl >/dev/null 2>&1 && curl -fsS "http://localhost:8000/health" >/dev/null 2>&1; then
  command -v xdg-open >/dev/null 2>&1 && xdg-open "http://localhost:8000" >/dev/null 2>&1 || true
  exit 0
fi

if command -v lsof >/dev/null 2>&1 && lsof -ti tcp:8000 >/dev/null 2>&1; then
  if command -v zenity >/dev/null 2>&1; then
    zenity --error --title="InterfaceScout" --text="Port 8000 is already in use by another process." 2>/dev/null
  else
    echo "Port 8000 is already in use by another process." >&2
  fi
  exit 1
fi

nohup python main.py >/tmp/interfacescout.log 2>&1 &
sleep 2

if command -v xdg-open >/dev/null 2>&1; then
  xdg-open "http://localhost:8000" >/dev/null 2>&1 || true
fi

exit 0
