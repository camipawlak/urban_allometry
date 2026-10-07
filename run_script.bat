@echo off
REM run_script.bat - Easily run Python scripts in the allometry conda environment
REM Usage: run_script.bat <script_name> [script_args...]
REM Example: run_script.bat 0_clean_names.py
REM Example: run_script.bat 0_clean_names.py --data-dir data/_raw

REM Get the directory of this script
set SCRIPT_DIR=%~dp0

REM Check if user provided a script name
if "%1"=="" (
    echo Usage: run_script.bat ^<script_name^> [script_args...]
    echo Example: run_script.bat 0_clean_names.py
    echo.
    echo Available scripts in scripts/py/:
    dir "%SCRIPT_DIR%scripts\py\*.py" /b
    exit /b 1
)

REM Activate the allometry conda environment and run the script
call "%USERPROFILE%\miniconda3\Scripts\activate.bat" allometry 2>nul
if errorlevel 1 (
    echo.
    echo ERROR: Could not activate 'allometry' conda environment.
    echo Please create it first by running:
    echo   conda env create -f environment.yml
    exit /b 1
)

REM Collect arguments after the script name (%* would pass the script name again)
set "SCRIPT=%~1"
set "ARGS="
shift
:collect
if "%~1"=="" goto run
set ARGS=%ARGS% %1
shift
goto collect

:run
REM Run the script from the repo root
cd /d "%SCRIPT_DIR%"
python "scripts\py\%SCRIPT%" %ARGS%
