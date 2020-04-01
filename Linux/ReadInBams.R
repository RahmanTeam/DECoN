packrat::on()
print("BEGIN ReadInBams.R")

library(R.utils)

########## Process inputs #######


args=commandArgs(asValue=TRUE)

bam_file=args$bams                                                                        #location of bam files; can be a directory containing only bam files to be processed or the name of a file containing a list of bam files to be processed.
bedfile=args$bed                                                                         #name of bed file
fasta=args$fasta                                                                        #name of fasta file
output=args$out                                                                        #location and name of file to save the output to; will be saved as 'output_counts.RData'


if(is.null(bam_file)){
print("ERROR bam files must be provided. Execution halted")
quit()
}


if(is.null(bedfile)){
print("ERROR bed file must be provided. Execution halted")
quit()
}



if(is.null(output)){output="DECoNBams"}


#################################
library(ExomeDepth)
library(parallel)
library(foreach)




myBamCounts <- function(bed.frame = NULL, bed.file = NULL, bam.files, index.files = bam.files,
                         min.mapq = 20, read.width = 300, include.chr = FALSE, referenceFasta = NULL) {

  if (is.null(bed.frame)) {
    if (is.null(bed.file)) {
      stop("If no bed data frame is provided there must be a link to a bed file")
    }
    bed.frame <- read.delim(file = bed.file, header =  FALSE, stringsAsFactors = FALSE)
  }

  names(bed.frame)[1] <- 'seqnames'
  names(bed.frame)[2] <- 'start'
  names(bed.frame)[3] <- 'end'

  if (include.chr) {
    if (sum(grepl(pattern = '^chr', bed.frame$seqnames) > 0)) {
      warning('The option include.chr == TRUE adds the chr prefix to the chromosome name but it looks like the chromosome names already have a chr prefix. The argument to getBamCounts is probably an error.')
    }
    bed.frame$seqnames <- paste('chr', bed.frame$seqnames, sep = '')
  }
  
  chr.names.used <- unique(as.character(bed.frame$seqnames))
  chr.levels <- c(as.character(seq(1, 22)), subset( chr.names.used, ! chr.names.used %in% as.character(seq(1, 22))))

  bed.frame$seqnames <- factor(bed.frame$seqnames, levels = chr.levels)  ####specifying the levels is important here to not mess up the order
  bed.frame <- bed.frame[ order(bed.frame$seqnames, bed.frame$start + bed.frame$end), ]  ##order the data frame by position

  target <- GenomicRanges::GRanges(seqnames = bed.frame$seqnames,  
                    IRanges::IRanges(start=bed.frame$start+1,end=bed.frame$end))
  
  rdata <- IRanges::RangedData(space= GenomicRanges::seqnames(target),
                               ranges=GenomicRanges::ranges(target))
  
  if  ((ncol(bed.frame) >= 4) && (class(bed.frame[,4]) %in% c('character', 'factor'))) {    
    row.names(rdata) <- make.unique(as.character(bed.frame[,4]))  ##add exon names if available
  }
  
############################################################################# add GC content
if (!is.null(referenceFasta)) {
  message('Reference fasta file provided so ExomeDepth will compute the GC content in each window')
    target.dnastringset <- Rsamtools::scanFa(referenceFasta, target)
  
    getGCcontent <- function(x) {
      GC.count <- Biostrings::letterFrequency(x,"GC")
      all.count <- Biostrings::letterFrequency(x,"ATGC")
      as.vector(ifelse(all.count==0,NA,GC.count/all.count))
    }
    rdata[["GC"]] <- getGCcontent(target.dnastringset)
  }

############################################################################# Parse BAM files
  nfiles <- length(bam.files)
  message('Parse ', nfiles, ' BAM files')
  print(bam.files)


  my.param <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isDuplicate = FALSE, isPaired = TRUE, isProperPair = TRUE, isNotPrimaryRead = FALSE),
                                      what = c("mapq", "pos", "isize"), )

  # for (i in 1:nfiles) {
  #   bam <- bam.files[ i ]
  #   index <- index.files[ i ]
  #   rdata[[ basename(bam) ]] <- countBamInGRanges.exomeDepth ( bam.file = bam, index = index, granges = target, min.mapq = min.mapq, read.width = read.width)
  #   message("Number of counted fragments : ", sum(rdata[[ basename(bam) ]]))
  # }
  numCores <- detectCores()
  cl <- parallel::makeForkCluster(numCores)
  doParallel::registerDoParallel(cl)
  r <- foreach(i=1:nfiles, .combine='cbind') %dopar% {
    bam <- bam.files[ i ]
    index <- index.files[ i ]
    ExomeDepth::countBamInGRanges.exomeDepth ( bam.file = bam, index = index, granges = target, min.mapq = min.mapq, read.width = read.width)
  }
  parallel::stopCluster(cl)
  colnames(r) <- basename(bam.files)
  rdata <- data.frame(rdata)
  rdata <- cbind(rdata, r)

  return(data.frame(rdata))
}


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
 
counts <- myBamCounts(bed.frame = bed.file, bam.files = bams, include.chr = FALSE, referenceFasta = paste(fasta))
                                                                                        #reads in coverage info from each bam file in 'bams'; expects chromosomes to be given as numbers, eg 1, 2 etc, not chr1, chr2 etc.
save(counts,bams,bed.file,sample.names,fasta,file=paste(output,".RData",sep=""))                                   #saves workspace as 'output_counts.RData'
   

print("END ReadInBams.R") 
