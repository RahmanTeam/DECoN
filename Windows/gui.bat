@echo off

set /p path=<DECoNPath.log

set /p "File=Enter full path to .RData file: "
Echo %File%

echo Launching gui in browser

R CMD BATCH --no-save "--args %File%" scripts/runShiny.R

pause
