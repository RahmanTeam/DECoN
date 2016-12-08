@echo off

set /p path=<DECoNPath.log

set /p "Rdata=Enter summary file: "
Echo %Rdata%

set /p mincorr=Enter correlation threshold (optional, default =.98): || set mincorr=0.98 
Echo %mincorr%

set /p mincov=Enter coverage threshold (optional, default=100): || set mincov=100 
Echo %mincov%

set /p exons=Enter custom exon numbering file (optional): || set exons=NULL 
Echo %exons%

set /p custom=Generate output for custom exons (TRUE or FALSE): || set custom=FALSE 
Echo %custom%

set brca=FALSE

set /p out=Enter output prefix (optional): || set out=DECoN
Echo %out%

Echo Running quality checks...


R CMD BATCH --no-save "--args %Rdata% %mincorr% %mincov% %exons% %custom% %brca% %out%" scripts\IdentifyFailures.R

pause
