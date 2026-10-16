param(
    [switch]$Clean,
    [string]$PythonExecutable
)

$ErrorActionPreference = "Stop"
$visualizeDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repositoryRoot = $visualizeDir
$interfaceMasterRoot = Join-Path (Split-Path -Parent $repositoryRoot) "interface_master"
$specPath = Join-Path $visualizeDir "gb_viewer.spec"
$distPath = Join-Path $visualizeDir "dist"
$workPath = Join-Path $visualizeDir "build"

if (Test-Path -LiteralPath $interfaceMasterRoot) {
    $env:PYTHONPATH = "$interfaceMasterRoot;$env:PYTHONPATH"
}

$pythonCommand = $PythonExecutable
if ([string]::IsNullOrWhiteSpace($pythonCommand)) {
    $python = Get-Command python -ErrorAction SilentlyContinue
    if ($null -eq $python) {
        throw "Python was not found. Activate the build environment first."
    }
    $pythonCommand = $python.Source
} elseif (-not (Test-Path -LiteralPath $pythonCommand)) {
    throw "Python executable was not found: $pythonCommand"
}

# Missing hidden imports can otherwise produce an apparently successful exe.
# Fail before packaging if any required runtime is absent from the build host.
& $pythonCommand -c "import ase, numpy, pymatgen, pyxtal, PySide6, pyvista, pyvistaqt, PyInstaller, interfacemaster"
if ($LASTEXITCODE -ne 0) {
    throw "A viewer dependency is missing. Install requirements-gb-viewer.txt in the build environment."
}

$arguments = @("--noconfirm")
if ($Clean) {
    $arguments += "--clean"
}
$arguments += @("--distpath", $distPath, "--workpath", $workPath)
$arguments += $specPath

Push-Location $repositoryRoot
try {
    & $pythonCommand -m PyInstaller @arguments
    if ($LASTEXITCODE -ne 0) {
        throw "PyInstaller failed with exit code $LASTEXITCODE"
    }
} finally {
    Pop-Location
}

Write-Host "Standalone viewer: $visualizeDir\dist\CCkitGBViewer\CCkitGBViewer.exe"
