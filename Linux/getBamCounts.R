library(parallel)
library(foreach)
numCores <- detectCores()

getBamCounts <- function(bed.frame = NULL, bed.file = NULL, bam.files, index.files = bam.files,
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
  r <- foreach(i=1:nfiles, .combine='cbind') %dopar% {
    bam <- bam.files[ i ]
    index <- index.files[ i ]
    countBamInGRanges.exomeDepth ( bam.file = bam, index = index, granges = target, min.mapq = min.mapq, read.width = read.width)
  }
  colnames(r) <- basename(bam.files)
  rdata <- cbind(rdata, r)

  return(data.frame(rdata))
}







getBamCounts <- function(bed.frame = NULL, bed.file = NULL, bam.files, index.files = bam.files,
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

  if  ((ncol(bed.frame) >= 4) && (class(bed.frame[,4]) %in% c('character', 'factor'))) {
      GenomicRanges::values(target) <- cbind(GenomicRanges::values(target), data.frame(exon = as.character(bed.frame[,4]),stringsAsFactors = FALSE))
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
  GenomicRanges::values(target) <- cbind(GenomicRanges::values(target), data.frame(GC = getGCcontent(target.dnastringset)))
}

############################################################################# Parse BAM files
  nfiles <- length(bam.files)
  message('Parse ', nfiles, ' BAM files')
  print(bam.files)


  my.param <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isDuplicate = FALSE, isPaired = TRUE, isProperPair = TRUE, isSecondaryAlignment = FALSE),
                                      what = c("mapq", "pos", "isize"), )


  exon_count_frame <- dplyr::tibble(chromosome = as(GenomicRanges::seqnames(target), 'character'),
                                    start = as(GenomicRanges::start(target), 'numeric'),
                                    end = as(GenomicRanges::end(target), 'numeric'))
  exon_count_frame <- dplyr::bind_cols(exon_count_frame, as(GenomicRanges::values(target), 'data.frame'))

  # for (i in 1:nfiles) {
  #   bam <- bam.files[ i ]
  #   index <- index.files[ i ]
  #   exon_count_frame[[ basename(bam) ]] <- countBamInGRanges.exomeDepth ( bam.file = bam, index = index, granges = target, min.mapq = min.mapq, read.width = read.width)
  #   message("Number of counted fragments : ", sum(exon_count_frame[[ basename(bam) ]]))
  # }

 }

getBamCounts <- function(bed.frame = NULL, bed.file = NULL, bam.files, index.files = bam.files,
                         min.mapq = 20, read.width = 300, include.chr = FALSE, referenceFasta = NULL) {
  if (is.null(bed.frame)) {
          if (is.null(bed.file)) {
              stop("If no bed data frame is provided there must be a link to a bed file")
          }
          bed.frame <- read.delim(file = bed.file, header = FALSE,
              stringsAsFactors = FALSE)
      }
      names(bed.frame)[1] <- "seqnames"
      names(bed.frame)[2] <- "start"
      names(bed.frame)[3] <- "end"
      if (include.chr) {
          if (sum(grepl(pattern = "^chr", bed.frame$seqnames) >
              0)) {
              warning("The option include.chr == TRUE adds the chr prefix to the chromosome name but it looks like the chromosome names already have a chr prefix. The argument to getBamCounts is probably an error.")
          }
          bed.frame$seqnames <- paste("chr", bed.frame$seqnames,
              sep = "")
      }
      chr.names.used <- unique(as.character(bed.frame$seqnames))
      chr.levels <- c(as.character(seq(1, 22)), subset(chr.names.used,
          !chr.names.used %in% as.character(seq(1, 22))))
      bed.frame$seqnames <- factor(bed.frame$seqnames, levels = chr.levels)
      bed.frame <- bed.frame[order(bed.frame$seqnames, bed.frame$start +
          bed.frame$end), ]
      target <- GenomicRanges::GRanges(seqnames = bed.frame$seqnames,
          IRanges::IRanges(start = bed.frame$start + 1, end = bed.frame$end))
      rdata <- IRanges::RangedData(space = GenomicRanges::seqnames(target),
          ranges = GenomicRanges::ranges(target))
      if ((ncol(bed.frame) >= 4) && (class(bed.frame[, 4]) %in%
          c("character", "factor"))) {
          row.names(rdata) <- make.unique(as.character(bed.frame[,
              4]))
      }
      if (!is.null(referenceFasta)) {
          message("Reference fasta file provided so ExomeDepth will compute the GC content in each window")
          target.dnastringset <- Rsamtools::scanFa(referenceFasta,
              target)
          getGCcontent <- function(x) {
              GC.count <- Biostrings::letterFrequency(x, "GC")
              all.count <- Biostrings::letterFrequency(x, "ATGC")
              as.vector(ifelse(all.count == 0, NA, GC.count/all.count))
          }
          rdata[["GC"]] <- getGCcontent(target.dnastringset)
      }
      nfiles <- length(bam.files)
      message("Parse ", nfiles, " BAM files")
      print(bam.files)
      my.param <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isDuplicate = FALSE,
          isPaired = TRUE, isProperPair = TRUE, isNotPrimaryRead = FALSE),
          what = c("mapq", "pos", "isize"), )
      for (i in 1:nfiles) {
          bam <- bam.files[i]
          index <- index.files[i]
          rdata[[basename(bam)]] <- countBamInGRanges.exomeDepth(bam.file = bam,
              index = index, granges = target, min.mapq = min.mapq,
              read.width = read.width)
          message("Number of counted fragments : ", sum(rdata[[basename(bam)]]))
      }
      return(rdata)
}





#' Count the number of everted reads for a set of BAM files.
#'
#' This is the ExomeDepth high level function that takes a GenomicRanges
#' object, a list of indexed/sorted BAM files, and compute the number of
#' everted reads in each of the defined bins.
#'
#' Everted reads are characteristic of the presence of duplications in a BAM
#' files. This routine will parse a BAM files and the suggested use is to
#' provide relatively large bins (for example gene based, and ExomeDepth has a
#' genes.hg19 object that is appropriate for this) to flag the genes that
#' contain such reads suggestive of a duplication. A manual check of the data
#' using IGV is recommended to confirm that these reads are all located in the
#' same DNA region, which would confirm the presence of a copy number variant.
#'
#' @param bed.frame \code{data.frame} containing the definition of the regions.
#' The first three columns must be chromosome, start, end.
#' @param bed.file \code{character} file name. Target BED file with the
#' definition of the regions. This file will only be used if no bed.frame
#' argument is provided. No headers are assumed so remove them if they exist.
#' Either a bed.file or a bed.frame must be provided for this function to run.
#' @param bam.files \code{character}, list of BAM files to extract read count
#' data from.
#' @param index.files Optional \code{character} argument with the list of
#' indexes for the BAM files, without the '.bai' suffix. If the indexes are
#' simply obtained by adding .bai to the BAM files, this argument does not need
#' to be specified.
#' @param min.mapq \code{numeric}, minimum mapping quality to include a read.
#' @param include.chr \code{logical}, if set to TRUE, this function will add
#' the string 'chr' to the chromosome names of the target BED file.
#' @return A data frame that contains the region and the number of identified
#' reads in each bin.
#' @note This function calls a lower level function called XXX that works on
#' each single BAM file.
#' @seealso getBAMCounts
#' @references Computational methods for discovering structural variation with
#' next-generation sequencing, Medvedev P, Stanciu M, Brudno M., Nature Methods
#' 2009
#' @examples
#'
#' \dontrun{  test <- count.everted.reads (bed.frame = genes.hg19,
#'   bed.file = NULL,
#'   bam.files = bam.files,
#'   min.mapq = 20,
#'   include.chr = FALSE)
#' }
#'

count.everted.reads <- function(bed.frame = NULL,
                                bed.file = NULL,
                                bam.files,
                                index.files = bam.files,
                                min.mapq = 20,
                                include.chr = FALSE) {

  if (is.null(bed.frame)) {
    if (is.null(bed.file)) {
      stop("If no bed data frame is provided there must be a link to a bed file")
    }
    bed.frame <- read.delim(file = bed.file, header =  FALSE, stringsAsFactors = FALSE)
  }

  names(bed.frame)[1] <- 'seqnames'
  names(bed.frame)[2] <- 'start'
  names(bed.frame)[3] <- 'end'


  if (include.chr) bed.frame$seqnames <- paste('chr', bed.frame$seqnames, sep = '')

  target <- GenomicRanges::GRanges(seqnames = bed.frame$seqnames,
                    IRanges::IRanges(start=bed.frame$start+1,end=bed.frame$end))

  if  ((ncol(bed.frame) >= 4) && (class(bed.frame[,4]) %in% c('character', 'factor'))) {
    GenomicRanges::values(target) <- cbind(GenomicRanges::values(target), data.frame(exon = as.character(bed.frame[,4]),stringsAsFactors = FALSE))
  }

  exon_count_frame <- dplyr::tibble(chromosome = as(GenomicRanges::seqnames(target), 'character'),
                                    start = as(GenomicRanges::start(target), 'numeric'),
                                    end = as(GenomicRanges::end(target), 'numeric'))
  exon_count_frame <- dplyr::bind_cols(exon_count_frame, as(GenomicRanges::values(target), 'data.frame'))

  nfiles <- length(bam.files)
  for (i in 1:nfiles) {
    bam <- bam.files[ i ]
    index <- index.files[ i ]
    exon_count_frame[[ basename(bam) ]] <- countBam.everted (bam.file = bam,  ## replace old RangedData with GRanges
                                                  index = index,
                                                  granges = target,
                                                  min.mapq = min.mapq)
   }

  return(data.frame(exon_count_frame))
}

