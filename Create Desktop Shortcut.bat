@echo off
setlocal
cd /d "%~dp0"
if /I "%~1"=="--check" (
  if exist "assets\InterfaceScout.ico" if exist "launchers\create_windows_shortcut.ps1" exit /b 0
  exit /b 1
)
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0launchers\create_windows_shortcut.ps1"
if errorlevel 1 (
  echo Could not create the shortcut. Make sure the ZIP was fully extracted.
  pause
  exit /b 1
)
echo.
echo Press any key to finish.
pause >nul
