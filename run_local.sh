#!/usr/bin/env bash
# ============================================================
# InterfaceScout - first-time setup (macOS / Linux)
# ============================================================
set -e

PROJ_DIR="$(cd "$(dirname "$0")" && pwd)"
BACKEND="$PROJ_DIR/backend"
FRONTEND="$PROJ_DIR/frontend"

if [ ! -f "$BACKEND/main.py" ] || [ ! -f "$BACKEND/v52_app.py" ] || [ ! -f "$FRONTEND/index.html" ]; then
  echo "ERROR: InterfaceScout files are incomplete."
  echo "Keep run_local.sh next to backend/ and frontend/."
  exit 1
fi

PYEXE=""
for c in python3.12 python3.11 python3.10 python3; do
  if command -v "$c" >/dev/null 2>&1 && "$c" -c "import ssl" >/dev/null 2>&1; then
    PYEXE="$c"
    break
  fi
done

if [ -z "$PYEXE" ]; then
  echo "ERROR: no supported Python with SSL found. Install Python 3.11 or 3.12."
  exit 1
fi

cd "$BACKEND"

if [ ! -d ".venv" ]; then
  echo "==> Creating virtual environment..."
  "$PYEXE" -m venv .venv
fi

# shellcheck disable=SC1091
source .venv/bin/activate

echo "==> Installing core dependencies..."
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

echo "==> Checking optional APBS support..."
if python -c "import apbs_binary" >/dev/null 2>&1; then
  echo "==> apbs-binary already installed."
else
  if python -m pip install "apbs-binary>=3.4.1.2"; then
    echo "==> Optional apbs-binary installed."
  else
    echo "WARNING: optional APBS binary could not be installed."
    echo "         Canonical InterfaceScout compatibility analysis will still run."
  fi
fi

chmod +x "$PROJ_DIR/start.command" "$PROJ_DIR/start.sh" 2>/dev/null || true

DESKTOP="$HOME/Desktop"
OS="$(uname -s)"

if [ -d "$DESKTOP" ]; then
  if [ "$OS" = "Darwin" ]; then
    LAUNCHER="$DESKTOP/InterfaceScout.command"
    {
      echo '#!/usr/bin/env bash'
      printf 'exec "%s/start.command"\n' "$PROJ_DIR"
    } > "$LAUNCHER"
    chmod +x "$LAUNCHER"
    echo "==> Desktop launcher created: InterfaceScout.command"
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

cd "$PROJ_DIR"
if [ "$OS" = "Darwin" ]; then
  exec "$PROJ_DIR/start.command"
else
  exec "$PROJ_DIR/start.sh"
fi
