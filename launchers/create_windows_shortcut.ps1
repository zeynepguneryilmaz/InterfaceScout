$ErrorActionPreference = "Stop"
$root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$runner = Join-Path $root "run_local.bat"
$icon = Join-Path $root "assets\InterfaceScout.ico"
if (!(Test-Path $runner) -or !(Test-Path $icon)) {
  throw "InterfaceScout launcher/icon is missing. Extract the Windows ZIP completely."
}
$desktop = [Environment]::GetFolderPath("Desktop")
if (!(Test-Path $desktop)) { throw "Windows Desktop directory not found." }
$shortcut = (New-Object -ComObject WScript.Shell).CreateShortcut((Join-Path $desktop "InterfaceScout.lnk"))
$shortcut.TargetPath = $runner
$shortcut.WorkingDirectory = $root
$shortcut.IconLocation = "$icon,0"
$shortcut.Description = "Launch InterfaceScout locally"
$shortcut.Save()
Write-Host "InterfaceScout desktop shortcut created:" (Join-Path $desktop "InterfaceScout.lnk")
Write-Host "Keep the extracted InterfaceScout directory at this location. Moving it breaks the shortcut."
