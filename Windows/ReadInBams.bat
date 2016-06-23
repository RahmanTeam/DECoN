@echo off

set /p path=<DECoNPath.log

set /p "Bams=Enter directory containing bams, or list of bam files (required): "
echo %Bams%

set /p "Bed=Enter bed file (required): "
echo %Bed%

set /p "Fasta=Enter fasta file (required): "
echo %Fasta%

set /p out=Enter output prefix (optional): || set out=DECoN
echo %out%

echo Calculating read depth... 

R CMD BATCH --no-save "--args %Bams% %Bed% %Fasta% %out%" scripts\ReadInBams.R

pause
