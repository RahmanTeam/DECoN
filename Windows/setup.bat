@echo off

setlocal enabledelayedexpansion

R --version > Rversion.log 2>&1
echo Checking your default version of R

echo %path% > path.log 2>&1

set Rversion=incorrect

findstr 3.1.2 Rversion.log
if ERRORLEVEL 0 set Rversion=default

echo %Rversion%

set defaultLocation="temp"
set actualLocation="doubletemp"

if not %Rversion%==default (
	echo "R version 3.1.2 is not your default R version"
	set /p "defaultLocation=Is R-3.1.2 installed in it's default location - C:\Program Files\R\R-3.1.2? Please enter y or n: "
)


if not %Rversion%==default (
	if %defaultLocation%==n (
		set /p "actualLocation=Enter path to R-3.1.2: "
		set path=!actualLocation!;!path!
		echo !path! > DECoNpath.log 2>&1
	) else (
		set actualLocation="C:\Program Files\R\R-3.1.2\bin"
		set path=!actualLocation!;!path!
		echo !path! > DECoNpath.log 2>&1
	)
)


if %Rversion%==default (
	echo !path! > DECoNpath.log 2>&1
)

if exist .Rprofile (del .Rprofile)

echo Running setup scripts...

R CMD BATCH scripts\sessionInfo.R

pause
