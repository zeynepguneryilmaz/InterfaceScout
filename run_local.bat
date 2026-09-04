@echo off
setlocal EnableExtensions EnableDelayedExpansion
REM InterfaceScout 2.0 - Windows setup + launcher
cd /d "%~dp0"
set "ROOT=%~dp0"
set "BACKEND=%ROOT%backend"
set "FRONTEND=%ROOT%frontend"
set "TRUSTED=--trusted-host pypi.org --trusted-host files.pythonhosted.org --trusted-host pypi.python.org"

if not exist "%BACKEND%\app.py" (
  echo ERROR: could not find backend\app.py
  pause
  exit /b 1
)
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
  pause
  exit /b 1
)

"%PYEXE%" -m venv "%BACKEND%\.venv"
if errorlevel 1 exit /b 1
call "%BACKEND%\.venv\Scripts\activate.bat"
python -m pip install --upgrade pip %TRUSTED%
python -m pip install -r "%BACKEND%\requirements.txt" %TRUSTED%
if errorlevel 1 (
  echo ERROR: dependency installation failed.
  pause
  exit /b 1
)

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
>> "%MKLNK%" echo   if (Test-Path $icon) { $lnk.IconLocation = $icon }
>> "%MKLNK%" echo   $lnk.Description = 'InterfaceScout 2.0'
>> "%MKLNK%" echo   $lnk.Save()
>> "%MKLNK%" echo }
powershell -NoProfile -ExecutionPolicy Bypass -File "%MKLNK%" >nul 2>&1
del "%MKLNK%" >nul 2>&1

goto :launch

:launch_existing
call "%BACKEND%\.venv\Scripts\activate.bat"
powershell -NoProfile -Command "try { $r=Invoke-WebRequest -UseBasicParsing 'http://127.0.0.1:8000/health' -TimeoutSec 2; if($r.StatusCode -eq 200){exit 0}else{exit 1} } catch { exit 1 }" >nul 2>&1
if !errorlevel! equ 0 (
  start "" "http://localhost:8000"
  exit /b 0
)

:launch
echo Starting InterfaceScout 2.0 at http://localhost:8000
cd /d "%BACKEND%"
start "" /b powershell -NoProfile -WindowStyle Hidden -Command "Start-Sleep -Seconds 2; Start-Process 'http://localhost:8000'"
python app.py
set "RC=%errorlevel%"
if not "%RC%"=="0" (
  python app.py > "%BACKEND%\startup_log.txt" 2>&1
  echo Log saved: %BACKEND%\startup_log.txt
  pause
)
exit /b %RC%
