@echo off
setlocal EnableExtensions
cd /d "%~dp0"
set "ROOT=%~dp0"
set "BACKEND=%ROOT%backend"
if not exist "%BACKEND%\v2\api.py" (
  echo ERROR: InterfaceScout backend is incomplete.
  pause
  exit /b 1
)
if not exist "%ROOT%frontend\index.html" (
  echo ERROR: InterfaceScout frontend is incomplete.
  pause
  exit /b 1
)

set "PYEXE="
for %%P in (python py) do (
  %%P -c "import ssl" >nul 2>&1
  if not errorlevel 1 if not defined PYEXE set "PYEXE=%%P"
)
if not defined PYEXE (
  echo ERROR: Python with SSL support is required.
  pause
  exit /b 1
)

if not exist "%BACKEND%\.venv\Scripts\activate.bat" (
  echo Creating InterfaceScout environment...
  %PYEXE% -m venv "%BACKEND%\.venv" || exit /b 1
  call "%BACKEND%\.venv\Scripts\activate.bat"
  python -m pip install --upgrade pip
  python -m pip install -r "%BACKEND%\requirements.txt" || exit /b 1
) else (
  call "%BACKEND%\.venv\Scripts\activate.bat"
)

powershell -NoProfile -Command "try { $r=Invoke-WebRequest -UseBasicParsing 'http://127.0.0.1:8000/health' -TimeoutSec 2; if($r.StatusCode -eq 200){exit 0}else{exit 1} } catch { exit 1 }" >nul 2>&1
if not errorlevel 1 (
  start "" "http://localhost:8000"
  exit /b 0
)

cd /d "%BACKEND%"
start "InterfaceScout" /min cmd /c "call .venv\Scripts\activate.bat && python -m uvicorn v2.api:app --host 127.0.0.1 --port 8000 > startup_log.txt 2>&1"
timeout /t 2 /nobreak >nul
start "" "http://localhost:8000"
exit /b 0
