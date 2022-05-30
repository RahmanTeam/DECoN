renv::restore()
print("BEGIN ReadInBams.R")

library(R.utils)


########## Process inputs #######


args=commandArgs(TRUE)

bam_file=args[1]                                                                        #location of bam files; can be a directory containing only bam files to be processed or the name of a file containing a list of bam files to be processed.
bedfile=args[2]                                                                         #name of bed file
fasta=args[3]                                                                        #name of fasta file
output=args[4]                                                                        #location and name of file to save the output to; will be saved as 'output_counts.RData'


if(bam_file=="NULL"){bam_file=NULL}
if(is.null(bam_file)){
print("ERROR bam files must be provided. Execution halted")
quit()
}


if(bedfile=="NULL"){bedfile=NULL}
if(is.null(bedfile)){
print("ERROR bed file must be provided. Execution halted")
quit()
}

if(fasta=="NULL"){fasta=NULL}
print(fasta)
print(paste(fasta))

if(output=="NULL"){output=NULL}
if(is.null(output)){output="DECoNBams"}


#################################

library(ExomeDepth)

multi_strsplit<-function(x,splits,y){                                                   #function which recursively splits x by an element of 'splits' then extracts the y element of the split vector
	X<-x
	for(i in 1:length(splits)){X=strsplit(X,splits[i])[[1]][y[i]]}
	return(X)
}


print(bam_file)

if(file.info(bam_file)$isdir){                                                          #if the input bam_file is a directory;
    bams<-list.files(bam_file,pattern=".bam",full.names=T)                            #puts all .bam files in that directory on the 'bams' list
    bais<-grep("bai",bams)
    if(length(bais)>0){
        bams<-bams[-bais]
    }
}else{                                                                                  #else expects bam_file to be a file containing a list of all the bam files to be used
    bams<-apply(read.table(paste(bam_file)),1,toString)                            #and reads in that list
}

if(length(bams)==0){print("ERROR NO BAM FILES DETECTED")}

a<-length(strsplit(bams[1],"/")[[1]])                                                   #works out where to split the paths of the bam files to get the sample names
sample.names<-sapply(bams,multi_strsplit,c("/",".bam"),c(a,1))                          #pulls out the sample names from the bam file paths
names(sample.names)<-NULL

bed.file<-read.table(paste(bedfile))                                                    #reads in the bedfile and gives each column a name - expects 4 columns: chr, start, stop, name/gene.
colnames(bed.file)<-c("chromosone","start","end","name")
 
counts <- getBamCounts(bed.frame = bed.file, bam.files = bams, include.chr = FALSE, referenceFasta = fasta)
                                                                                        #reads in coverage info from each bam file in 'bams'; expects chromosomes to be given as numbers, eg 1, 2 etc, not chr1, chr2 etc.
save(counts,bams,bed.file,sample.names,fasta,file=paste(output,".RData",sep=""))                                   #saves workspace as 'output_counts.RData'
   

print("END ReadInBams.R") 
