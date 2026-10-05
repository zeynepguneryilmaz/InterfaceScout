@echo off
setlocal EnableExtensions
cd /d "%~dp0"
set "ROOT=%~dp0"
set "BACKEND=%ROOT%backend"

if /I "%~1"=="--check" (echo InterfaceScout Windows launcher check passed& exit /b 0)
if not exist "%BACKEND%\.venv\Scripts\python.exe" (echo ERROR: First-time setup is required. Run run_local.bat.& exit /b 1)
if not exist "%ROOT%frontend\index.html" (echo ERROR: InterfaceScout frontend is incomplete.& exit /b 1)

call "%BACKEND%\.venv\Scripts\activate.bat"
python -c "import urllib.request; r=urllib.request.urlopen('http://127.0.0.1:8000/health',timeout=2); raise SystemExit(0 if r.status==200 else 1)" >nul 2>&1
if not errorlevel 1 (start "" "http://localhost:8000"& exit /b 0)

cd /d "%BACKEND%"
start "InterfaceScout" /min cmd /c "call .venv\Scripts\activate.bat && python -m uvicorn app:app --host 127.0.0.1 --port 8000 > startup_log.txt 2>&1"
for /L %%I in (1,1,15) do (
  timeout /t 1 /nobreak >nul
  python -c "import urllib.request; r=urllib.request.urlopen('http://127.0.0.1:8000/health',timeout=1); raise SystemExit(0 if r.status==200 else 1)" >nul 2>&1
  if not errorlevel 1 (start "" "http://localhost:8000"& exit /b 0)
)
echo ERROR: InterfaceScout did not start. Check backend\startup_log.txt.
exit /b 1
