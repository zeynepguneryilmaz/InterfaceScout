@echo off
setlocal enabledelayedexpansion
REM ============================================================
REM InterfaceScout - One-Click Local Setup (Windows)
REM ============================================================
REM Handles everything automatically:
REM   - finds a working Python (prefers 3.11/3.12 with SSL)
REM   - creates a virtual environment
REM   - installs Python deps (with corporate-proxy SSL workaround)
REM   - downloads + extracts APBS for Windows (no manual step)
REM   - creates a Desktop icon, then launches the app
REM ------------------------------------------------------------
REM This file sits in the InterfaceScout folder, NEXT TO backend\.
REM ============================================================

cd /d "%~dp0"
echo.
echo InterfaceScout launcher starting...
echo Folder: %~dp0
echo.
set "ROOT=%~dp0"
set "BACKEND=%ROOT%backend"
set "APBSDIR=%ROOT%apbs-win"
set "TRUSTED=--trusted-host pypi.org --trusted-host files.pythonhosted.org --trusted-host pypi.python.org"

if not exist "%BACKEND%\main.py" (
  echo ERROR: could not find backend\main.py
  echo Looked in: %BACKEND%
  echo.
  echo run_local.bat must sit in the InterfaceScout folder, next to the
  echo backend\ and frontend\ folders. This file is at:
  echo   %~f0
  echo.
  pause
  exit /b 1
)

REM ---- FAST PATH: if already set up, skip installation and just launch ----
REM This lets the Desktop shortcut point at run_local.bat: the first run does
REM the full setup, every later run detects the existing environment and starts
REM immediately (no re-install, no waiting).
if exist "%BACKEND%\.venv\Scripts\activate.bat" goto :launch_existing

echo ============================================================
echo   InterfaceScout - setup starting
echo ============================================================
echo.

REM ---- 1. Find a usable Python (need SSL + version 3.10-3.12) ----
set "PYEXE="
for %%P in (
  "%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
  "%LOCALAPPDATA%\Programs\Python\Python311\python.exe"
  "%LOCALAPPDATA%\Programs\Python\Python310\python.exe"
) do (
  if exist %%P (
    %%P -c "import ssl" >nul 2>&1
    if !errorlevel! equ 0 (
      set "PYEXE=%%~P"
      goto :found_py
    )
  )
)
REM fall back to PATH python if it has ssl
python -c "import ssl" >nul 2>&1
if %errorlevel% equ 0 set "PYEXE=python"

:found_py
if "%PYEXE%"=="" (
  echo ERROR: No Python with working SSL found.
  echo Please install Python 3.11 from https://www.python.org/downloads/
  echo and check "Add python.exe to PATH" during install.
  pause
  exit /b 1
)
echo ==^> Using Python: %PYEXE%
"%PYEXE%" --version

REM ---- 2. Create virtual environment ----
if not exist "%BACKEND%\.venv" (
  echo ==^> Creating virtual environment...
  "%PYEXE%" -m venv "%BACKEND%\.venv"
)
call "%BACKEND%\.venv\Scripts\activate.bat"

REM ---- 3. Install Python dependencies (skip apbs-binary on Windows) ----
echo ==^> Upgrading pip...
python -m pip install --upgrade pip %TRUSTED% >nul 2>&1

echo ==^> Installing dependencies (a few minutes the first time)...
pip install fastapi "uvicorn[standard]" pydantic numpy biopython pdb2pqr %TRUSTED%
if %errorlevel% neq 0 (
  echo ERROR: dependency install failed. See messages above.
  pause
  exit /b 1
)

REM ---- 4. Ensure APBS is available ----
REM Best path on Windows: install via conda. The conda-forge APBS package
REM brings ALL its dependencies (runtime DLLs, python3x.dll, suitesparse, ...)
REM so it "just works", unlike the standalone zip which fails with missing
REM python39.dll / 0xc000007b. If conda is present we use it; otherwise we
REM fall back to downloading the standalone build.
set "APBSEXE="

REM 4a. already on PATH and runnable (e.g. previously conda-installed)?
where apbs >nul 2>&1
if !errorlevel! equ 0 (
  for /f "delims=" %%F in ('where apbs') do if "!APBSEXE!"=="" set "APBSEXE=%%F"
)

REM 4b. if not found, try conda
if "%APBSEXE%"=="" (
  where conda >nul 2>&1
  if !errorlevel! equ 0 (
    echo ==^> Installing APBS via conda ^(conda-forge^). This may take a few minutes...
    call conda install -y -c conda-forge apbs
    where apbs >nul 2>&1
    if !errorlevel! equ 0 (
      for /f "delims=" %%F in ('where apbs') do if "!APBSEXE!"=="" set "APBSEXE=%%F"
    )
  )
)

REM 4c. last resort: standalone download (note: may need python39.dll)
if "%APBSEXE%"=="" (
  echo ==^> conda not available - downloading standalone APBS 3.4.1...
  set "APBSZIP=%ROOT%APBS-3.4.1.Windows.zip"
  set "APBSURL=https://github.com/Electrostatics/apbs/releases/download/v3.4.1/APBS-3.4.1.Windows.zip"
  powershell -Command "[Net.ServicePointManager]::SecurityProtocol=[Net.SecurityProtocolType]::Tls12; try { Invoke-WebRequest -Uri '!APBSURL!' -OutFile '!APBSZIP!' } catch { Write-Host 'Download failed:' $_.Exception.Message }"
  if exist "!APBSZIP!" (
    powershell -Command "Expand-Archive -Path '!APBSZIP!' -DestinationPath '%APBSDIR%' -Force"
    del "!APBSZIP!" >nul 2>&1
    for /r "%APBSDIR%" %%F in (apbs.exe) do (echo %%F | findstr /I "\bin\" >nul && set "APBSEXE=%%F")
  )
)

if "%APBSEXE%"=="" (
  echo WARNING: APBS not available. Backend will use Debye-Huckel fallback.
) else (
  echo ==^> APBS: %APBSEXE%
  set "APBS_PATH=%APBSEXE%"
  "%APBSEXE%" --version >nul 2>&1
  if !errorlevel! GEQ 3221225472 (
    echo WARNING: APBS found but failed to start. Using Debye-Huckel fallback.
    echo If you have conda, run:  conda install -c conda-forge apbs
  ) else (
    echo ==^> APBS runs correctly.
  )
)

REM ---- 4b. Create a desktop shortcut with the InterfaceScout icon ----
echo.
echo Creating desktop shortcut...
REM The shortcut points at run_local.bat (in THIS folder) and runs it minimized.
REM We create it on BOTH the OneDrive Desktop and the classic Desktop (whichever
REM exist), because Windows may redirect the Desktop into OneDrive. PowerShell
REM reports back whether it succeeded.
set "MKLNK=%TEMP%\_iscout_mklnk.ps1"
>  "%MKLNK%" echo $ErrorActionPreference = 'SilentlyContinue'
>> "%MKLNK%" echo $target = '%ROOT%run_local.bat'
>> "%MKLNK%" echo $workdir = '%ROOT%'
>> "%MKLNK%" echo $icon = '%ROOT%interfacescout.ico'
>> "%MKLNK%" echo $targets = @()
>> "%MKLNK%" echo $d1 = [Environment]::GetFolderPath('Desktop')
>> "%MKLNK%" echo if ($d1 -and (Test-Path $d1)) { $targets += $d1 }
>> "%MKLNK%" echo $d2 = Join-Path $env:USERPROFILE 'Desktop'
>> "%MKLNK%" echo if ((Test-Path $d2) -and ($targets -notcontains $d2)) { $targets += $d2 }
>> "%MKLNK%" echo $d3 = Join-Path $env:USERPROFILE 'OneDrive\Desktop'
>> "%MKLNK%" echo if ((Test-Path $d3) -and ($targets -notcontains $d3)) { $targets += $d3 }
>> "%MKLNK%" echo $made = 0
>> "%MKLNK%" echo $ws = New-Object -ComObject WScript.Shell
>> "%MKLNK%" echo foreach ($dt in $targets) {
>> "%MKLNK%" echo   try {
>> "%MKLNK%" echo     $lnk = $ws.CreateShortcut((Join-Path $dt 'InterfaceScout.lnk'))
>> "%MKLNK%" echo     $lnk.TargetPath = $target
>> "%MKLNK%" echo     $lnk.WorkingDirectory = $workdir
>> "%MKLNK%" echo     $lnk.WindowStyle = 7
>> "%MKLNK%" echo     if (Test-Path $icon) { $lnk.IconLocation = $icon }
>> "%MKLNK%" echo     $lnk.Description = 'InterfaceScout - Protein Surface Analysis'
>> "%MKLNK%" echo     $lnk.Save()
>> "%MKLNK%" echo     if (Test-Path (Join-Path $dt 'InterfaceScout.lnk')) { $made++ }
>> "%MKLNK%" echo   } catch {}
>> "%MKLNK%" echo }
>> "%MKLNK%" echo if ($made -gt 0) { Write-Output 'SHORTCUT_OK' } else { Write-Output 'SHORTCUT_FAIL' }
for /f "delims=" %%R in ('powershell -NoProfile -ExecutionPolicy Bypass -File "%MKLNK%" 2^>nul') do set "LNKRESULT=%%R"
del "%MKLNK%" 2>nul
if "%LNKRESULT%"=="SHORTCUT_OK" (
  echo ==^> Desktop shortcut created: InterfaceScout
  echo     Double-click it next time to start the app.
) else (
  echo Note: could not create a Desktop shortcut automatically on this PC.
  echo That is fine - you can always start the app by double-clicking
  echo   run_local.bat
  echo in this folder. ^(You may right-click it ^> Send to ^> Desktop to make
  echo your own shortcut.^)
)

REM ---- 5. Launch backend (it serves the UI and opens the browser) ----
echo.
echo ============================================================
echo   Setup complete. Starting the app...
echo   Your browser will open automatically at http://localhost:8000
echo   If it doesn't, open that address manually.
echo   Press Ctrl+C here to stop. Next time, use the Desktop icon
echo   to start the app again.
echo ============================================================
echo.
cd /d "%BACKEND%"
python main.py
set "RC=%errorlevel%"
echo.
if not "%RC%"=="0" (
  echo ============================================================
  echo   The app stopped with an error ^(code %RC%^).
  echo   Writing a diagnostic log to startup_log.txt ...
  python main.py > "%BACKEND%\startup_log.txt" 2>&1
  echo   Log saved: %BACKEND%\startup_log.txt
  echo   Please send that file if you need help.
  echo ============================================================
)
echo.
echo This window can be closed.
pause
exit /b 0

REM ============================================================
REM :launch_existing - fast path when .venv already exists
REM ============================================================
:launch_existing
echo ==^> Existing installation found. Starting the app...
echo     ^(To force a clean reinstall, delete the backend\.venv folder.^)
call "%BACKEND%\.venv\Scripts\activate.bat"

REM pick APBS from PATH or local folder if present
set "APBSEXE="
where apbs >nul 2>&1
if !errorlevel! equ 0 (
  for /f "delims=" %%F in ('where apbs') do if "!APBSEXE!"=="" set "APBSEXE=%%F"
)
if "!APBSEXE!"=="" if exist "%APBSDIR%\*" (
  for /r "%APBSDIR%" %%F in (apbs.exe) do (echo %%F | findstr /I "\bin\" >nul && set "APBSEXE=%%F")
)
if not "!APBSEXE!"=="" set "APBS_PATH=!APBSEXE!"

REM free port 8000 if still bound from a previous run
for /f "tokens=5" %%P in ('netstat -ano ^| findstr ":8000 " ^| findstr LISTENING') do taskkill /F /PID %%P >nul 2>&1

echo.
echo ============================================================
echo   Starting... your browser will open at http://localhost:8000
echo   If it doesn't, open that address manually.
echo   Press Ctrl+C here to stop.
echo ============================================================
echo.
cd /d "%BACKEND%"
python main.py
set "RC=%errorlevel%"
echo.
if not "%RC%"=="0" (
  echo The app stopped with an error ^(code %RC%^). Writing startup_log.txt ...
  python main.py > "%BACKEND%\startup_log.txt" 2>&1
  echo Log saved: %BACKEND%\startup_log.txt
)
echo.
echo This window can be closed.
pause
exit /b 0
