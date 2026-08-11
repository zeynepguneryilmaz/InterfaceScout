#!/usr/bin/env bash
# ============================================================
# InterfaceScout - One-Click Local Setup (macOS / Linux)
# ============================================================
# First-time setup:
#   - finds Python with SSL
#   - creates backend/.venv
#   - installs core Python dependencies
#   - attempts optional apbs-binary installation
#   - creates a Desktop launcher
#   - starts InterfaceScout locally
#
# This file must sit in the InterfaceScout folder, next to backend/
# and frontend/.
# ============================================================
set -e

PROJ_DIR="$(cd "$(dirname "$0")" && pwd)"
BACKEND="$PROJ_DIR/backend"
FRONTEND="$PROJ_DIR/frontend"

if [ ! -f "$BACKEND/main.py" ]; then
  echo "ERROR: could not find backend/main.py"
  echo "Looked in: $BACKEND"
  echo "run_local.sh must sit in the InterfaceScout folder, next to backend/ and frontend/."
  exit 1
fi

if [ ! -f "$FRONTEND/index.html" ]; then
  echo "ERROR: could not find frontend/index.html"
  echo "Looked in: $FRONTEND"
  exit 1
fi

cd "$BACKEND"

TRUSTED="--trusted-host pypi.org --trusted-host files.pythonhosted.org --trusted-host pypi.python.org"

PYEXE=""
for c in python3.12 python3.11 python3.10 python3; do
  if command -v "$c" >/dev/null 2>&1 && "$c" -c "import ssl" >/dev/null 2>&1; then
    PYEXE="$c"
    break
  fi
done

if [ -z "$PYEXE" ]; then
  echo "ERROR: no Python 3 with SSL found. Install Python 3.11 or 3.12."
  exit 1
fi

echo "==> Using $($PYEXE --version) ($PYEXE)"

if [ ! -d ".venv" ]; then
  echo "==> Creating virtual environment (.venv)..."
  "$PYEXE" -m venv .venv
fi

# shellcheck disable=SC1091
source .venv/bin/activate

echo "==> Upgrading pip..."
python -m pip install --upgrade pip $TRUSTED >/dev/null

echo "==> Installing core dependencies..."
python -m pip install -r requirements.txt $TRUSTED

echo ""
echo "==> Attempting optional APBS binary installation..."
if python -m pip install apbs-binary $TRUSTED; then
  echo "==> Optional APBS binary package installed."
else
  echo "WARNING: apbs-binary could not be installed on this platform."
  echo "         Canonical InterfaceScout compatibility analysis will still run."
  echo "         Only optional APBS electrostatic descriptors may be unavailable."
fi

echo ""
echo "==> Checking optional computational binaries..."
python - <<'PY'
import shutil
print("  pdb2pqr:", shutil.which("pdb2pqr") or "NOT FOUND")
try:
    import apbs_binary
    print("  apbs:   ", getattr(apbs_binary, "APBS_BIN_PATH", "package found"))
except Exception:
    print("  apbs:    NOT FOUND (optional)")
PY

chmod +x "$PROJ_DIR/start.command" "$PROJ_DIR/start.sh" 2>/dev/null || true
DESKTOP="$HOME/Desktop"

if [ -d "$DESKTOP" ]; then
  OS="$(uname -s)"

  if [ "$OS" = "Darwin" ]; then
    # Do not copy the real start.command to Desktop because it resolves
    # backend/ relative to itself. Create a tiny wrapper pointing back
    # to the project launcher instead.
    LAUNCHER="$DESKTOP/InterfaceScout.command"
    {
      echo '#!/usr/bin/env bash'
      printf 'exec "%s/start.command"\n' "$PROJ_DIR"
    } > "$LAUNCHER"
    chmod +x "$LAUNCHER"
    echo "==> Desktop launcher created: InterfaceScout.command"
    echo "    It points back to the InterfaceScout project folder."
    echo "    To customize its icon: Finder > Get Info and use interfacescout.png."

  else
    LAUNCHER="$DESKTOP/InterfaceScout.desktop"
    cat > "$LAUNCHER" <<EOF
[Desktop Entry]
Version=1.0
Type=Application
Name=InterfaceScout
Comment=Protein-surface chemistry and multiscale compatibility mapping
Exec=bash "$PROJ_DIR/start.sh"
Icon=$PROJ_DIR/interfacescout.png
Terminal=false
Categories=Science;Education;
StartupNotify=true
EOF
    chmod +x "$LAUNCHER"
    command -v gio >/dev/null 2>&1 && gio set "$LAUNCHER" "metadata::trusted" true 2>/dev/null || true
    echo "==> Desktop launcher created: InterfaceScout.desktop"
  fi
fi

echo ""
echo "============================================================"
echo "  Setup complete."
echo "  InterfaceScout will start at http://localhost:8000"
echo "  On later runs, use the Desktop launcher."
echo "============================================================"
echo ""

nohup python main.py >/tmp/interfacescout.log 2>&1 &
sleep 2

OS="$(uname -s)"
if [ "$OS" = "Darwin" ]; then
  open "http://localhost:8000" 2>/dev/null || true
elif command -v xdg-open >/dev/null 2>&1; then
  xdg-open "http://localhost:8000" >/dev/null 2>&1 || true
fi

echo "InterfaceScout started. Log: /tmp/interfacescout.log"
exit 0
