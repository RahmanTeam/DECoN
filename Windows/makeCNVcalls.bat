@echo off

set /p path=<DECoNPath.log

set /p "Rdata=Enter summary file: "

echo %Rdata%

set /p tp=Enter transition probability (optional, default =.01): || set tp=0.01

echo %tp%

set /p exons=Enter custom exon numbering file (optional): || set exons=NULL 

echo %exons%

set /p custom=Generate output for custom exons (TRUE or FALSE)? || set custom=FALSE 

echo %custom%

set brca=FALSE

set /p plot=Generate plots of variants (None, Custom, All)? || set plot=All

echo %plot%

set /p plotFolder=Enter folder to contain plots: || set plotFolder=DECoNPlots

echo %plotFolder%

set /p out=Enter output prefix (optional): || set out=DECoN

echo %out%

echo Identifying exon deletions/duplications...

R CMD BATCH "--args %Rdata% %tp% %exons% %custom% %brca% %plot% %plotFolder% %out%" scripts\makeCNVcalls.R 


pause
