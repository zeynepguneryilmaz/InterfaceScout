@echo off
setlocal EnableExtensions
cd /d "%~dp0"
set "ROOT=%~dp0"
set "BACKEND=%ROOT%backend"
if not exist "%BACKEND%\app.py" (
  echo ERROR: InterfaceScout backend is incomplete.
  exit /b 1
)
if not exist "%ROOT%frontend\index.html" (
  echo ERROR: InterfaceScout frontend is incomplete.
  exit /b 1
)
set "PYEXE="
for %%P in (python py) do (
  %%P -c "import sys,ssl; raise SystemExit(0 if sys.version_info[:2] in [(3,10),(3,11),(3,12)] else 1)" >nul 2>&1
  if not errorlevel 1 if not defined PYEXE set "PYEXE=%%P"
)
if not defined PYEXE (
  echo ERROR: Python 3.10, 3.11, or 3.12 with SSL support is required.
  exit /b 1
)
if /I "%~1"=="--check" (
  %PYEXE% -c "import sys; print('InterfaceScout Windows launcher check:', sys.version.split()[0])"
  exit /b 0
)
if not exist "%BACKEND%\.venv\Scripts\activate.bat" (
  echo Creating InterfaceScout environment...
  %PYEXE% -m venv "%BACKEND%\.venv" || exit /b 1
)
call "%BACKEND%\.venv\Scripts\activate.bat"
python -m pip install --upgrade pip
python -m pip install -r "%BACKEND%\requirements.txt" || exit /b 1
powershell -NoProfile -Command "try { $r=Invoke-WebRequest -UseBasicParsing 'http://127.0.0.1:8000/health' -TimeoutSec 2; if($r.StatusCode -eq 200){exit 0}else{exit 1} } catch { exit 1 }" >nul 2>&1
if not errorlevel 1 (
  start "" "http://localhost:8000"
  exit /b 0
)
cd /d "%BACKEND%"
start "InterfaceScout" /min cmd /c "call .venv\Scripts\activate.bat && python -m uvicorn app:app --host 127.0.0.1 --port 8000 > startup_log.txt 2>&1"
for /L %%I in (1,1,10) do (
  timeout /t 1 /nobreak >nul
  powershell -NoProfile -Command "try { $r=Invoke-WebRequest -UseBasicParsing 'http://127.0.0.1:8000/health' -TimeoutSec 1; if($r.StatusCode -eq 200){exit 0}else{exit 1} } catch { exit 1 }" >nul 2>&1
  if not errorlevel 1 (
    start "" "http://localhost:8000"
    exit /b 0
  )
)
echo ERROR: InterfaceScout did not start. Check backend\startup_log.txt.
exit /b 1
