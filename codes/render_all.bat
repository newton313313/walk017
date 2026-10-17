@echo off
REM ---------------------------------------------------------------------------
REM render_all.bat
REM One-click rendering of the videos into codes\renders\ (wraps render_all.py).
REM
REM   Double-click it            -> a small menu (scenes+demo / demo / scenes / 480p)
REM   render_all.bat <args>      -> passes the args to render_all.py, no menu,
REM                                 e.g.  render_all.bat W105 W110 -l
REM
REM Existing videos are skipped (finished videos are never overwritten);
REM pass --force to render them again. If manimgl is not on PATH (double-click
REM from Explorer), this script tries to activate Anaconda / Miniconda by itself.
REM Keep this file ASCII-only (cmd reads .bat files as cp950).
REM ---------------------------------------------------------------------------
setlocal
cd /d "%~dp0"

where manimgl >nul 2>nul || call :activate_conda
where manimgl >nul 2>nul || goto :nomanim

if not "%~1"=="" (
    python render_all.py %*
    goto :done
)

echo.
echo ============================================================
echo   One-click render  (output: codes\renders\)
echo ============================================================
echo   1. Everything: scenes W101-W112, W201-W212, W301-W308, W401-W408 + Ch1Demo..Ch4Demo, 1080p
echo   2. Demos only: Ch1Demo (W101-W112), Ch2Demo (W201-W212), Ch3Demo (W301-W308), Ch4Demo (W401-W408), 1080p
echo   3. Scenes only: the 40 single-scene videos, 1080p
echo   4. Quick check: everything at 480p (fast, low quality)
echo   Q. Quit
echo.
echo   Videos that already exist are skipped (add --force to redo).
echo.
set "CHOICE=none"
set /p "CHOICE=Choose 1-4 or Q: "
set "CHOICE=%CHOICE: =%"

set DEMOS=--ch1-demo --ch2-demo --ch3-demo --ch4-demo

if /i "%CHOICE%"=="Q" goto :eof
if "%CHOICE%"=="1" python render_all.py --hd %DEMOS% & goto :done
if "%CHOICE%"=="2" python render_all.py --hd %DEMOS% --demos-only & goto :done
if "%CHOICE%"=="3" python render_all.py --hd & goto :done
if "%CHOICE%"=="4" python render_all.py -l %DEMOS% & goto :done
echo Unknown choice "%CHOICE%". Nothing was rendered.
pause
exit /b 1

:activate_conda
for %%D in ("%USERPROFILE%\anaconda3" "%USERPROFILE%\miniconda3" "%ProgramData%\anaconda3" "%ProgramData%\miniconda3" "%LOCALAPPDATA%\anaconda3" "%LOCALAPPDATA%\miniconda3") do (
    if exist "%%~D\Scripts\activate.bat" (
        call "%%~D\Scripts\activate.bat" "%%~D"
        goto :eof
    )
)
goto :eof

:nomanim
echo.
echo manimgl was not found.
echo Start this file from an "Anaconda Prompt", or install ManimGL first.
pause
exit /b 1

:done
echo.
echo Finished. Videos are in %~dp0renders
pause
