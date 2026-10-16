@echo off
setlocal
if exist "%~dp0..\interface_master\interfacemaster" (
    set "PYTHONPATH=%~dp0..\interface_master;%PYTHONPATH%"
)
python "%~dp0gb_viewer.py" %*
endlocal
