@echo off
setlocal
set "EXE=%~dp0dist\CCkitGBViewer\CCkitGBViewer.exe"
if not exist "%EXE%" (
    echo Standalone viewer not found. Run build_gb_viewer.ps1 first.
    exit /b 1
)
"%EXE%" %*
endlocal
