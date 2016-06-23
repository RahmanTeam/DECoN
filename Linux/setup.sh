#!/bin/bash

[ -w ".Rprofile" ] && rm .Rprofile

Rscript sessionInfo.R --bootstrap-packrat > setup.log 2>&1

cp packrat/packrat_source/.Rprofile ./
