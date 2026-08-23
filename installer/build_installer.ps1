# Build the native TEagle Windows installer end to end: freeze the PySide6 app (PyInstaller) then
# compile the Inno Setup installer. Version is derived from teagle_core (single source of truth).
# Requires: PyInstaller (pip) + Inno Setup 6. Output: dist\TEagle-Setup-<ver>.exe (one file to share).
param(
  [string]$PythonExe = "",
  [switch]$PreflightOnly
)

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
$init = Join-Path $root "app\backend\teagle_core\__init__.py"
$m = Select-String -Path $init -Pattern '__version__\s*=\s*["'']([^"'']+)["'']'
if (-not $m) { throw "could not parse __version__ from $init" }
$ver = $m.Matches[0].Groups[1].Value
Write-Output "TEagle version: $ver"

function Resolve-BuildPython([string]$Requested) {
  $candidates = @()
  if ($Requested) {
    $candidates += $Requested
  } else {
    # A Microsoft Store alias may precede a real interpreter on PATH. Prefer per-user Python installs,
    # then consider executable applications resolved from PATH.
    $localPythonRoot = Join-Path $env:LOCALAPPDATA "Programs\Python"
    if (Test-Path -LiteralPath $localPythonRoot) {
      $candidates += Get-ChildItem -Path (Join-Path $localPythonRoot "Python*\python.exe") -File -ErrorAction SilentlyContinue |
                     Sort-Object FullName -Descending | Select-Object -ExpandProperty FullName
    }
    foreach ($commandName in @("python.exe", "python3.exe")) {
      $resolvedCommands = Get-Command $commandName -CommandType Application -ErrorAction SilentlyContinue
      if ($resolvedCommands) { $candidates += $resolvedCommands | Select-Object -ExpandProperty Source }
    }
  }
  foreach ($candidate in ($candidates | Select-Object -Unique)) {
    try {
      $resolved = (Resolve-Path -LiteralPath $candidate -ErrorAction Stop).Path
      & $resolved -c "import sys; raise SystemExit(0 if sys.version_info[:2] == (3, 12) else 12)" 2>$null
      if ($LASTEXITCODE -eq 0) { return $resolved }
    } catch {
      continue
    }
  }
  if ($Requested) { throw "-PythonExe '$Requested' is not a runnable Python 3.12 executable" }
  throw "Python 3.12 was not found. Pass -PythonExe C:\path\to\python.exe or install Python 3.12."
}

$python = Resolve-BuildPython $PythonExe
$pythonVersion = (& $python -c "import platform; print(platform.python_version())").Trim()
Write-Output "Build Python: $python ($pythonVersion)"

$runtimeImportGate = Join-Path $PSScriptRoot "check_runtime_imports.py"
if (-not (Test-Path -LiteralPath $runtimeImportGate)) { throw "Runtime import parity gate is missing: $runtimeImportGate" }
& $python $runtimeImportGate --root $root
if ($LASTEXITCODE -ne 0) { throw "Runtime import parity preflight failed" }

# Canonical runtime gate: envcheck reads app/backend/requirements.txt, requires exact versions, imports
# every direct dependency, and enforces its minimum Python version. It never installs during a build.
$envcheck = Join-Path $root "app\backend\envcheck.py"
& $python $envcheck --no-install
if ($LASTEXITCODE -ne 0) {
  throw "Runtime dependency preflight failed. Install the exact pins with: '$python' -m pip install -r app\backend\requirements.txt"
}

$buildLock = Join-Path $PSScriptRoot "requirements-build.txt"
$pyiLine = Get-Content -LiteralPath $buildLock | Where-Object { $_ -match '^\s*PyInstaller==([^\s#]+)\s*$' } | Select-Object -First 1
if (-not $pyiLine) { throw "requirements-build.txt must contain exactly one PyInstaller==<version> pin" }
$pyiRequired = ([regex]::Match($pyiLine, '^\s*PyInstaller==([^\s#]+)\s*$')).Groups[1].Value
$pyiInstalled = & $python -c "import importlib.metadata as m; print(m.version('PyInstaller'))" 2>$null
if ($LASTEXITCODE -ne 0) { throw "PyInstaller $pyiRequired is required but is not installed in $python" }
$pyiInstalled = ($pyiInstalled | Select-Object -Last 1).Trim()
if ($pyiInstalled -ne $pyiRequired) {
  throw "PyInstaller version mismatch: required $pyiRequired from requirements-build.txt, found $pyiInstalled"
}
Write-Output "PyInstaller: $pyiInstalled (locked)"

$candidates = @(
  "$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe",
  "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe",
  "$env:ProgramFiles\Inno Setup 6\ISCC.exe"
)
$iscc = $candidates | Where-Object { Test-Path $_ } | Select-Object -First 1
if (-not $iscc) { throw "ISCC.exe not found. Install Inno Setup 6 (winget install JRSoftware.InnoSetup)." }
$innoLines = @(Get-Content -LiteralPath $buildLock | Where-Object { $_ -match '^\s*InnoSetup==([^\s#]+)\s*$' })
if ($innoLines.Count -ne 1) { throw "requirements-build.txt must contain exactly one InnoSetup==<version> pin" }
$innoRequired = ([regex]::Match($innoLines[0], '^\s*InnoSetup==([^\s#]+)\s*$')).Groups[1].Value
$innoInstall = Split-Path -Parent $iscc
$innoRegistry = @(
  Get-ItemProperty 'HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\Inno Setup 6_is1' -ErrorAction SilentlyContinue
  Get-ItemProperty 'HKLM:\Software\Microsoft\Windows\CurrentVersion\Uninstall\Inno Setup 6_is1' -ErrorAction SilentlyContinue
  Get-ItemProperty 'HKLM:\Software\WOW6432Node\Microsoft\Windows\CurrentVersion\Uninstall\Inno Setup 6_is1' -ErrorAction SilentlyContinue
) | Where-Object {
  $_.DisplayVersion -and $_.InstallLocation -and
  ([IO.Path]::GetFullPath($_.InstallLocation).TrimEnd('\') -eq [IO.Path]::GetFullPath($innoInstall).TrimEnd('\'))
} | Select-Object -First 1
if (-not $innoRegistry) { throw "Cannot verify the installed Inno Setup version for $iscc" }
if ($innoRegistry.DisplayVersion -ne $innoRequired) {
  throw "Inno Setup version mismatch: required $innoRequired from requirements-build.txt, found $($innoRegistry.DisplayVersion)"
}
Write-Output "ISCC: $iscc (Inno Setup $innoRequired, locked)"

if ($PreflightOnly) {
  Write-Output "Preflight passed; no build requested."
  exit 0
}

# 1) always freeze the native app into dist\TEagle. A controlled release never packages a cached exe.
$bundle = Join-Path $root "dist\TEagle"
# PyInstaller preserves read-only attributes from package metadata. A later --noconfirm rebuild must be
# able to replace its own generated bundle, so clear that attribute inside this exact build output only.
if (Test-Path -LiteralPath $bundle) {
  @(Get-Item -LiteralPath $bundle) + @(Get-ChildItem -LiteralPath $bundle -Recurse -Force) | ForEach-Object {
    if ($_.Attributes -band [IO.FileAttributes]::ReadOnly) {
      $_.Attributes = $_.Attributes -band (-bnot [IO.FileAttributes]::ReadOnly)
    }
  }
}
$dist = Join-Path $root "dist\TEagle\TEagle.exe"
Write-Output "Freezing native app (PyInstaller)..."
Push-Location $root
& $python -m PyInstaller "installer\teagle_native.spec" --noconfirm
Pop-Location
if ($LASTEXITCODE -ne 0) { throw "PyInstaller failed with exit code $LASTEXITCODE" }
$metadataHelper = Join-Path $PSScriptRoot "sanitize_inno_metadata.py"
if (-not (Test-Path -LiteralPath $metadataHelper)) { throw "Installer metadata helper is missing: $metadataHelper" }
& $python $metadataHelper sanitize-rich-header $dist
if ($LASTEXITCODE -ne 0) { throw "Frozen executable metadata minimisation failed with exit code $LASTEXITCODE" }
# sanity: the frozen bundle must import the scientific stack + render QtSvg before we ship it
$env:TEAGLE_SELFTEST = "1"; $env:QT_QPA_PLATFORM = "offscreen"
# Start-Process -Wait, NOT '& $dist': the exe is GUI-subsystem (console=False), so the call operator returns
  # immediately and leaves $LASTEXITCODE stale — the gate below silently never fired.
$st = (Start-Process -FilePath $dist -Wait -PassThru -NoNewWindow).ExitCode
Remove-Item Env:\TEAGLE_SELFTEST; Remove-Item Env:\QT_QPA_PLATFORM
if ($st -ne 0) { throw "frozen-bundle self-test failed (exit $st) - not shipping a broken build" }
Write-Output "Frozen-bundle self-test passed."
if (-not (Test-Path $dist)) { throw "dist\TEagle\TEagle.exe not found - PyInstaller step did not produce it" }

$innoSource = Split-Path -Parent $iscc
$innoTemp = Join-Path ([IO.Path]::GetTempPath()) ("TEagle-Inno-" + [guid]::NewGuid().ToString("N"))
$setupOutput = Join-Path $root ("dist\TEagle-Setup-$ver.exe")
New-Item -ItemType Directory -Path $innoTemp | Out-Null
try {
  # Inno's stock Setup/uninstaller templates carry a generator Comments field and generator file version.
  # Compile from an isolated copy in which those exact spans are replaced by blank metadata and the TEagle
  # release version without changing PE layout. The helper verifies every byte delta before compilation.
  Get-ChildItem -LiteralPath $innoSource -Force | Copy-Item -Destination $innoTemp -Recurse -Force
  & $python $metadataHelper prepare-templates $innoSource $innoTemp $ver
  if ($LASTEXITCODE -ne 0) { throw "Inno template metadata minimisation failed with exit code $LASTEXITCODE" }
  $sanitizedIscc = Join-Path $innoTemp "ISCC.exe"
  if (-not (Test-Path -LiteralPath $sanitizedIscc)) { throw "Copied ISCC.exe is missing: $sanitizedIscc" }
  & $sanitizedIscc "/DMyAppVersion=$ver" (Join-Path $PSScriptRoot "teagle.iss")
  if ($LASTEXITCODE -ne 0) { throw "ISCC failed with exit code $LASTEXITCODE" }
  & $python $metadataHelper verify-setup $setupOutput $ver
  if ($LASTEXITCODE -ne 0) { throw "Installer metadata verification failed with exit code $LASTEXITCODE" }
} finally {
  if (Test-Path -LiteralPath $innoTemp) {
    $resolvedTemp = (Resolve-Path -LiteralPath $innoTemp).Path
    $tempRoot = [IO.Path]::GetFullPath([IO.Path]::GetTempPath()).TrimEnd('\') + '\'
    if ($resolvedTemp.StartsWith($tempRoot, [StringComparison]::OrdinalIgnoreCase) -and
        (Split-Path -Leaf $resolvedTemp).StartsWith("TEagle-Inno-")) {
      Remove-Item -LiteralPath $resolvedTemp -Recurse -Force
    } else {
      Write-Warning "Refusing to remove unexpected Inno work directory: $resolvedTemp"
    }
  }
}
Write-Output "Installer written to: $setupOutput"
