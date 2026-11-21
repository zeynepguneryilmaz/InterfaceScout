@echo off
setlocal EnableExtensions EnableDelayedExpansion
cd /d "%~dp0"
set "ROOT=%~dp0"
set "BACKEND=%ROOT%backend"
set "NO_PAUSE=0"
if /I "%~1"=="--no-pause" set "NO_PAUSE=1"

if not exist "%BACKEND%\app.py" (
  echo ERROR: InterfaceScout backend is incomplete.
  goto :failed
)
if not exist "%ROOT%frontend\index.html" (
  echo ERROR: InterfaceScout frontend is incomplete.
  goto :failed
)

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
if not defined PYEXE (
  echo ERROR: Python 3.10, 3.11, or 3.12 with SSL support is required.
  goto :failed
)

if /I "%~1"=="--check" (
  %PYEXE% %PYARGS% -c "import sys; print('InterfaceScout Windows setup check:', sys.version.split()[0])"
  exit /b %errorlevel%
)

if not exist "%BACKEND%\.venv\Scripts\python.exe" (
  echo Creating InterfaceScout environment...
  %PYEXE% %PYARGS% -m venv "%BACKEND%\.venv"
  if errorlevel 1 goto :failed
)
set "VENV_PY=%BACKEND%\.venv\Scripts\python.exe"
if not exist "%VENV_PY%" (
  echo ERROR: Virtual environment Python is missing.
  goto :failed
)

rem Do not contact PyPI on subsequent launches if the pinned versions are present.
"%VENV_PY%" "%BACKEND%\dependency_check.py" >nul 2>&1
if errorlevel 1 (
  echo Installing missing or mismatched InterfaceScout dependencies...
  "%VENV_PY%" -m pip install --disable-pip-version-check -r "%BACKEND%\requirements.txt"
  if errorlevel 1 goto :dependency_failed
  "%VENV_PY%" "%BACKEND%\dependency_check.py"
  if errorlevel 1 goto :dependency_failed
) else (
  echo Pinned dependencies are already installed; skipping pip download.
)

call "%ROOT%start.bat"
if errorlevel 1 goto :failed
exit /b 0

:dependency_failed
echo.
echo ERROR: Python dependency installation failed.
echo If the output contains CERTIFICATE_VERIFY_FAILED / self signed certificate:
echo   The Python package download certificate chain is not trusted.
echo   This often occurs behind corporate or university HTTPS inspection.
echo   Obtain a verified CA certificate bundle from your institution's IT.
echo   In a new Command Prompt, run:
echo.
echo     set "PIP_CERT=C:\path\to\trusted-ca-bundle.pem"
echo     run_local.bat
echo.
echo   Replace the example with the REAL path to the verified PEM bundle.
echo   If you already have an IT-approved wheelhouse, install dependencies
echo   locally using the instructions in docs\INSTALL_WINDOWS.md.
echo   Do NOT disable SSL verification or use pip --trusted-host.
echo.
echo   Details and offline-install options: docs\INSTALL_WINDOWS.md
goto :failed

:failed
echo.
echo InterfaceScout setup did not complete. Review the messages above.
if "%NO_PAUSE%"=="0" (
  echo Press any key to close this window...
  pause >nul
)
exit /b 1
