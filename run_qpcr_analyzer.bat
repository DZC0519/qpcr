@echo off
setlocal
cd /d "%~dp0"
set "PYTHONPATH=%~dp0src"
if exist "%~dp0.venv\Scripts\pythonw.exe" (
    "%~dp0.venv\Scripts\pythonw.exe" -m qpcr_analyzer
) else (
    pythonw -m qpcr_analyzer
)
endlocal
