@echo off
setlocal EnableExtensions
cd /d "%~dp0"
set "PY=backend\.venv\Scripts\python.exe"
if not exist "%PY%" (
  echo InterfaceScout base environment missing. Run run_local.bat once first.
  exit /b 1
)
echo Installing OPTIONAL PDB2PQR / PROPKA and SciPy...
"%PY%" -m pip install --disable-pip-version-check "pdb2pqr==3.7.1" "scipy>=1.14,<1.18"
if errorlevel 1 (
  echo Optional Python dependencies failed. Base InterfaceScout still works.
  exit /b 1
)
echo.
where apbs >nul 2>&1
if errorlevel 1 (
  echo APBS executable was not found in PATH.
  echo Download/install the APBS CLI from its official distribution,
  echo then reopen Command Prompt and check: apbs --version
  echo The optional panel will remain disabled for APBS until it is installed.
) else (
  echo APBS executable found:
  where apbs
)
"%PY%" -c "import pdb2pqr,scipy;print('PDB2PQR + SciPy Python imports OK')"
echo.
echo Reopen InterfaceScout with run_local.bat.
endlocal
