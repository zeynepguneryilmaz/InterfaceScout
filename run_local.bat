@echo off
setlocal EnableExtensions EnableDelayedExpansion
REM ============================================================
REM InterfaceScout - Windows setup + daily launcher
REM ============================================================
REM Keep this file in the InterfaceScout folder, next to:
REM   backend\   frontend\   interfacescout.ico
REM
REM First run:
REM   - finds Python 3.10-3.12 with SSL
REM   - creates backend\.venv
REM   - installs core dependencies
REM   - attempts optional APBS setup
REM   - creates a Desktop shortcut
REM
REM Later runs:
REM   - reuses the existing environment
REM   - starts InterfaceScout immediately
REM ============================================================

cd /d "%~dp0"
set "ROOT=%~dp0"
set "BACKEND=%ROOT%backend"
set "FRONTEND=%ROOT%frontend"
set "APBSDIR=%ROOT%apbs-win"
set "TRUSTED=--trusted-host pypi.org --trusted-host files.pythonhosted.org --trusted-host pypi.python.org"

if not exist "%BACKEND%\main.py" (
  echo ERROR: could not find backend\main.py
  pause
  exit /b 1
)
if not exist "%FRONTEND%\index.html" (
  echo ERROR: could not find frontend\index.html
  pause
  exit /b 1
)

if exist "%BACKEND%\.venv\Scripts\activate.bat" goto :launch_existing

echo.
echo ============================================================
echo   InterfaceScout - first-time setup
echo ============================================================
echo.

REM ---- Find Python with SSL; prefer 3.12 / 3.11 ----
set "PYEXE="
for %%P in (
  "%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
  "%LOCALAPPDATA%\Programs\Python\Python311\python.exe"
  "%LOCALAPPDATA%\Programs\Python\Python310\python.exe"
) do (
  if exist "%%~P" (
    "%%~P" -c "import ssl" >nul 2>&1
    if !errorlevel! equ 0 (
      set "PYEXE=%%~P"
      goto :found_py
    )
  )
)

python -c "import ssl" >nul 2>&1
if !errorlevel! equ 0 set "PYEXE=python"

:found_py
if "%PYEXE%"=="" (
  echo ERROR: No supported Python with working SSL was found.
  echo Install Python 3.11 or 3.12 and enable "Add python.exe to PATH".
  pause
  exit /b 1
)

echo ==^> Using:
"%PYEXE%" --version

echo ==^> Creating virtual environment...
"%PYEXE%" -m venv "%BACKEND%\.venv"
if errorlevel 1 (
  echo ERROR: virtual environment creation failed.
  pause
  exit /b 1
)

call "%BACKEND%\.venv\Scripts\activate.bat"

echo ==^> Upgrading pip...
python -m pip install --upgrade pip %TRUSTED%
if errorlevel 1 (
  echo ERROR: pip upgrade failed.
  pause
  exit /b 1
)

echo ==^> Installing core dependencies...
python -m pip install -r "%BACKEND%\requirements.txt" %TRUSTED%
if errorlevel 1 (
  echo ERROR: dependency installation failed.
  pause
  exit /b 1
)

REM ---- Optional APBS ----
set "APBSEXE="
where apbs >nul 2>&1
if !errorlevel! equ 0 (
  for /f "delims=" %%F in ('where apbs') do if "!APBSEXE!"=="" set "APBSEXE=%%F"
)

if "!APBSEXE!"=="" (
  where conda >nul 2>&1
  if !errorlevel! equ 0 (
    echo ==^> Optional: attempting APBS installation via conda-forge...
    call conda install -y -c conda-forge apbs
    where apbs >nul 2>&1
    if !errorlevel! equ 0 (
      for /f "delims=" %%F in ('where apbs') do if "!APBSEXE!"=="" set "APBSEXE=%%F"
    )
  )
)

if "!APBSEXE!"=="" (
  echo WARNING: APBS is not available.
  echo          Canonical InterfaceScout compatibility analysis will still run.
  echo          Only optional APBS electrostatic descriptors will be unavailable.
) else (
  echo ==^> Optional APBS found: !APBSEXE!
  set "APBS_PATH=!APBSEXE!"
)

REM ---- Desktop shortcut ----
echo ==^> Creating Desktop shortcut...
set "MKLNK=%TEMP%\_iscout_mklnk.ps1"
>  "%MKLNK%" echo $ErrorActionPreference = 'SilentlyContinue'
>> "%MKLNK%" echo $target = '%ROOT%run_local.bat'
>> "%MKLNK%" echo $workdir = '%ROOT%'
>> "%MKLNK%" echo $icon = '%ROOT%interfacescout.ico'
>> "%MKLNK%" echo $desktop = [Environment]::GetFolderPath('Desktop')
>> "%MKLNK%" echo if ($desktop) {
>> "%MKLNK%" echo   $ws = New-Object -ComObject WScript.Shell
>> "%MKLNK%" echo   $lnk = $ws.CreateShortcut((Join-Path $desktop 'InterfaceScout.lnk'))
>> "%MKLNK%" echo   $lnk.TargetPath = $target
>> "%MKLNK%" echo   $lnk.WorkingDirectory = $workdir
>> "%MKLNK%" echo   $lnk.WindowStyle = 7
>> "%MKLNK%" echo   if (Test-Path $icon) { $lnk.IconLocation = $icon }
>> "%MKLNK%" echo   $lnk.Description = 'InterfaceScout'
>> "%MKLNK%" echo   $lnk.Save()
>> "%MKLNK%" echo }
powershell -NoProfile -ExecutionPolicy Bypass -File "%MKLNK%" >nul 2>&1
del "%MKLNK%" >nul 2>&1

goto :launch

:launch_existing
echo ==^> Existing InterfaceScout environment found.
call "%BACKEND%\.venv\Scripts\activate.bat"

REM If our server is already running, only open the browser.
powershell -NoProfile -Command "try { $r=Invoke-WebRequest -UseBasicParsing 'http://127.0.0.1:8000/health' -TimeoutSec 2; if($r.StatusCode -eq 200){exit 0}else{exit 1} } catch { exit 1 }" >nul 2>&1
if !errorlevel! equ 0 (
  start "" "http://localhost:8000"
  exit /b 0
)

REM Refuse to kill an unrelated process using port 8000.
set "PORTPID="
for /f "tokens=5" %%P in ('netstat -ano ^| findstr ":8000 " ^| findstr LISTENING') do set "PORTPID=%%P"
if defined PORTPID (
  echo ERROR: Port 8000 is already in use by another process ^(PID !PORTPID!^).
  echo Close that process, then run InterfaceScout again.
  pause
  exit /b 1
)

REM Restore optional APBS path if available.
where apbs >nul 2>&1
if !errorlevel! equ 0 (
  for /f "delims=" %%F in ('where apbs') do (
    set "APBS_PATH=%%F"
    goto :apbs_done
  )
)
:apbs_done

:launch
echo.
echo ============================================================
echo   Starting InterfaceScout
echo   http://localhost:8000
echo ============================================================
echo.

cd /d "%BACKEND%"

REM Open the browser after the server has had time to start.
start "" /b powershell -NoProfile -WindowStyle Hidden -Command "Start-Sleep -Seconds 2; Start-Process 'http://localhost:8000'"

python main.py
set "RC=%errorlevel%"

if not "%RC%"=="0" (
  echo.
  echo InterfaceScout stopped with error code %RC%.
  echo Re-running once to write backend\startup_log.txt ...
  python main.py > "%BACKEND%\startup_log.txt" 2>&1
  echo Log saved: %BACKEND%\startup_log.txt
  pause
)

exit /b %RC%
