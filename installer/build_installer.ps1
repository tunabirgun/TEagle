# Build the native TEagle Windows installer end to end: freeze the PySide6 app (PyInstaller) then
# compile the Inno Setup installer. Version is derived from teagle_core (single source of truth).
# Requires: PyInstaller (pip) + Inno Setup 6. Output: dist\TEagle-Setup-<ver>.exe (one file to share).
$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
$init = Join-Path $root "app\backend\teagle_core\__init__.py"
$m = Select-String -Path $init -Pattern '__version__\s*=\s*["'']([^"'']+)["'']'
if (-not $m) { throw "could not parse __version__ from $init" }
$ver = $m.Matches[0].Groups[1].Value
Write-Output "TEagle version: $ver"

# 1) freeze the native app into dist\TEagle (skip with -SkipFreeze if already built)
$dist = Join-Path $root "dist\TEagle\TEagle.exe"
if (($args -notcontains "-SkipFreeze") -or (-not (Test-Path $dist))) {
  Write-Output "Freezing native app (PyInstaller)..."
  Push-Location $root
  & python -m PyInstaller "installer\teagle_native.spec" --noconfirm
  Pop-Location
  if ($LASTEXITCODE -ne 0) { throw "PyInstaller failed with exit code $LASTEXITCODE" }
  # sanity: the frozen bundle must import the scientific stack + render QtSvg before we ship it
  $env:TEAGLE_SELFTEST = "1"; $env:QT_QPA_PLATFORM = "offscreen"
  & $dist
  $st = $LASTEXITCODE
  Remove-Item Env:\TEAGLE_SELFTEST; Remove-Item Env:\QT_QPA_PLATFORM
  if ($st -ne 0) { throw "frozen-bundle self-test failed (exit $st) - not shipping a broken build" }
  Write-Output "Frozen-bundle self-test passed."
}
if (-not (Test-Path $dist)) { throw "dist\TEagle\TEagle.exe not found - PyInstaller step did not produce it" }

$candidates = @(
  "$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe",
  "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe",
  "$env:ProgramFiles\Inno Setup 6\ISCC.exe"
)
$iscc = $candidates | Where-Object { Test-Path $_ } | Select-Object -First 1
if (-not $iscc) { throw "ISCC.exe not found. Install Inno Setup 6 (winget install JRSoftware.InnoSetup)." }
Write-Output "ISCC: $iscc"

& $iscc "/DMyAppVersion=$ver" (Join-Path $PSScriptRoot "teagle.iss")
if ($LASTEXITCODE -ne 0) { throw "ISCC failed with exit code $LASTEXITCODE" }
Write-Output "Installer written to: $(Join-Path $root ("dist\TEagle-Setup-$ver.exe"))"
