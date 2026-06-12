#!/usr/bin/env bash
# ============================================================
# InterfaceScout - One-Click Local Setup (macOS / Linux)
# ============================================================
# Creates a venv, installs everything (APBS via the apbs-binary pip
# package, which has Linux + macOS wheels), creates a Desktop icon,
# and starts the backend.
# This file sits in the InterfaceScout folder, next to backend/.
# ============================================================
set -e
PROJ_DIR="$(cd "$(dirname "$0")" && pwd)"

if [ ! -f "$PROJ_DIR/backend/main.py" ]; then
  echo "ERROR: could not find backend/main.py"
  echo "Looked in: $PROJ_DIR/backend"
  echo "run_local.sh must sit in the InterfaceScout folder, next to backend/."
  exit 1
fi
cd "$PROJ_DIR/backend"

TRUSTED="--trusted-host pypi.org --trusted-host files.pythonhosted.org --trusted-host pypi.python.org"

# ---- pick a python with ssl, prefer 3.11/3.12 ----
PYEXE=""
for c in python3.12 python3.11 python3.10 python3; do
  if command -v "$c" >/dev/null 2>&1 && "$c" -c "import ssl" >/dev/null 2>&1; then
    PYEXE="$c"; break
  fi
done
if [ -z "$PYEXE" ]; then
  echo "ERROR: no Python 3 with SSL found. Install Python 3.11+."; exit 1
fi
echo "==> Using $($PYEXE --version) ($PYEXE)"

echo "==> Creating virtual environment (.venv)..."
"$PYEXE" -m venv .venv
# shellcheck disable=SC1091
source .venv/bin/activate

echo "==> Upgrading pip..."
pip install --upgrade pip $TRUSTED >/dev/null

echo "==> Installing dependencies (a few minutes the first time)..."
pip install -r requirements.txt $TRUSTED

echo ""
echo "==> Checking computational binaries..."
python - <<'PY'
import shutil
print("  pdb2pqr:", shutil.which("pdb2pqr") or "NOT FOUND")
try:
    import apbs_binary
    print("  apbs:   ", apbs_binary.APBS_BIN_PATH)
except Exception as e:
    print("  apbs:    NOT FOUND -", e)
PY

# ---- Create a desktop launcher with the InterfaceScout icon ----
chmod +x "$PROJ_DIR/start.command" "$PROJ_DIR/start.sh" 2>/dev/null || true
DESKTOP="$HOME/Desktop"
if [ -d "$DESKTOP" ]; then
  OS="$(uname -s)"
  if [ "$OS" = "Darwin" ]; then
    cp "$PROJ_DIR/start.command" "$DESKTOP/InterfaceScout.command"
    chmod +x "$DESKTOP/InterfaceScout.command"
    echo "==> Desktop launcher created: InterfaceScout.command"
    echo "    (To set its icon: right-click > Get Info, drag interfacescout.png onto the icon.)"
  else
    LAUNCHER="$DESKTOP/InterfaceScout.desktop"
    cat > "$LAUNCHER" <<EOF
[Desktop Entry]
Version=1.0
Type=Application
Name=InterfaceScout
Comment=Protein Surface Analysis - residue-level surface chemistry mapping
Exec=bash "$PROJ_DIR/start.sh"
Icon=$PROJ_DIR/interfacescout.png
Terminal=false
Categories=Science;Education;
StartupNotify=true
EOF
    chmod +x "$LAUNCHER"
    gio set "$LAUNCHER" "metadata::trusted" true 2>/dev/null || true
    echo "==> Desktop launcher created: InterfaceScout.desktop (no terminal window)"
  fi
fi

echo ""
echo "============================================================"
echo "  Setup complete. Starting... browser opens at http://localhost:8000"
echo "  If it doesn't, open that address manually. Press Ctrl+C to stop."
echo "  Next time, use the Desktop icon."
echo "============================================================"
echo ""
python main.py
