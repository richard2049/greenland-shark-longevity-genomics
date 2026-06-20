@echo off
setlocal EnableDelayedExpansion

for %%I in ("%~dp0..") do set "REPO=%%~fI"
set "SALMON_ARGS="

:arg_loop
if "%~1"=="" goto run_salmon
set "ARG=%~1"
set "ARG=!ARG:\=/!"
set SALMON_ARGS=!SALMON_ARGS! "!ARG!"
shift
goto arg_loop

:run_salmon
docker --context desktop-linux run --rm -v "%REPO%:/work" -w /work combinelab/salmon:latest salmon %SALMON_ARGS%
exit /b %ERRORLEVEL%
