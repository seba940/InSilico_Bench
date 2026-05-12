@echo off
set PYTHONPATH=%PYTHONPATH%;%CD%

echo Checking for required Python libraries...
python -c "import Bio; import pandas; import primer3; import openpyxl" 2>nul
if %errorlevel% neq 0 (
    echo Installing missing libraries. This may take a few minutes...
    pip install biopython pandas primer3-py openpyxl
    if %errorlevel% neq 0 (
        echo.
        echo Error: Failed to install libraries. Please check your internet connection or Python installation.
        pause
        exit /b
    )
    echo Installation complete.
)

echo Starting InSilico_Bench...
start /b pythonw main.py
exit
