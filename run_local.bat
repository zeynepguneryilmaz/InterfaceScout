@echo off
setlocal EnableExtensions EnableDelayedExpansion
cd /d "%~dp0"
set "ROOT=%~dp0"
set "BACKEND=%ROOT%backend"

if not exist "%BACKEND%\app.py" (echo ERROR: InterfaceScout backend is incomplete.& exit /b 1)
if not exist "%ROOT%frontend\index.html" (echo ERROR: InterfaceScout frontend is incomplete.& exit /b 1)

set "PYEXE="
set "PYARGS="
for %%V in (3.12 3.11 3.10) do (
  if not defined PYEXE (
    py -%%V -c "import sys,ssl; raise SystemExit(0 if sys.version_info[:2]==tuple(map(int,'%%V'.split('.'))) else 1)" >nul 2>&1
    if not errorlevel 1 (set "PYEXE=py"& set "PYARGS=-%%V")
  )
)
if not defined PYEXE (
  python -c "import sys,ssl; raise SystemExit(0 if (3,10) <= sys.version_info[:2] <= (3,12) else 1)" >nul 2>&1
  if not errorlevel 1 set "PYEXE=python"
)
if not defined PYEXE (echo ERROR: Python 3.10, 3.11, or 3.12 with SSL support is required.& exit /b 1)

if /I "%~1"=="--check" (
  %PYEXE% %PYARGS% -c "import sys; print('InterfaceScout Windows setup check:', sys.version.split()[0])"
  exit /b 0
)

if not exist "%BACKEND%\.venv\Scripts\python.exe" (
  echo Creating InterfaceScout environment...
  %PYEXE% %PYARGS% -m venv "%BACKEND%\.venv" || exit /b 1
)
call "%BACKEND%\.venv\Scripts\activate.bat"
python -m pip install --disable-pip-version-check -r "%BACKEND%\requirements.txt" || exit /b 1
call "%ROOT%start.bat"
