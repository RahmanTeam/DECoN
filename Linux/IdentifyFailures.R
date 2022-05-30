renv::restore()

print("BEGIN IdentifyFailures.R")


library(R.utils)
library(optparse)

option_list<-list(
    make_option("--RData",help="Summary RData file",dest='Rdata'),
    make_option("--mincorr",help='minimum correlation, default .98',default=.98,dest='mincorr'),
    make_option("--mincov",help='minimum coverage, default=100',default=100,dest='mincov'),
    make_option("--exons",default=NULL,help="File containing custom exon annotation",dest='exons'),
    make_option("--custom",default=FALSE,help='Use custom reporting, default=FALSE',dest='custom'),
    make_option("--out",default='DECoN',help='output prefix,default=DECoN',dest='out')
)

opt<-parse_args(OptionParser(option_list=option_list))


count_data=opt$Rdata

if(count_data=="NULL"){count_data=NULL}

if(is.null(count_data)){
print("ERROR: no Rdata summary file provided -- Execution halted")
quit()
}																#R workspace with the coverage data saved in - this workspace already has the bedfile, fasta file saved in it.
load(count_data)



corr_thresh=as.numeric(opt$mincorr)

if(length(corr_thresh)==0){corr_thresh=0.98}

cov_thresh=as.numeric(opt$mincov)

if(length(cov_thresh)==0){cov_thresh=100}


exon_numbers=opt$exons

Custom=as.logical(opt$custom)
if(length(Custom)==0){Custom=FALSE}

BRCA=FALSE

OutFails=opt$out

library(ExomeDepth)

ExomeCount<-as(counts, 'data.frame')											#converts counts, a ranged data object, to a data frame

ExomeCount$chromosome <- gsub(as.character(ExomeCount$chromosome),pattern = 'chr',replacement = '') 
																				#remove any chr letters, and coerce to a string.
colnames(ExomeCount)[1:length(sample.names)+4]=sample.names						#assigns the sample names to each column 

ReadDepths<-ExomeCount[,sample.names]											#extracts just the read depths


Sample<-vector()
Exon<-vector()
Details<-vector()
Types<-vector()
Gene<-vector()


Corr<-cor(ReadDepths)															#calculates correlation matrix
MaxCorr<-apply(Corr,1,function(x)max(x[x!=1]))								#finds the maximum correlation for each sample

for(i in 1:length(MaxCorr)){													#tests correlation for each sample; if below .98, adds to list of fails
	if(MaxCorr[i]<corr_thresh){
		
			Sample<-c(Sample,sample.names[i])
			Exon<-c(Exon,"All")
			Types<-c(Types,"Whole sample")
			Details<-c(Details,paste("Low correlation: ", MaxCorr[i],sep=""))
			Gene<-c(Gene,"All")
			
	}

}



SampleMedian<-apply(ReadDepths,2,median)											#calculates median coverage per sample

for(i in 1:length(SampleMedian)){													#tests median coverage for each sample; if below 100x, adds to list of fails
	if(SampleMedian[i]<cov_thresh){
		
		if(sample.names[i]%in%Sample){
			
			k=which(Sample==sample.names[i])
			Details[k] = paste(Details[k],", Low median read depth (FPKM): ", SampleMedian[i],sep="")
			
		}else{
			Sample<-c(Sample,sample.names[i])
			Exon<-c(Exon,"All")
			Types<-c(Types,"Whole sample")
			Details<-c(Details,paste("Low median read depth (FPKM): ", SampleMedian[i],sep=""))
			Gene<-c(Gene,"All")
		}
	}
	
}


	
ExonMedian<-apply(ReadDepths,1,median)												#calculates median coverage per exon

for(i in 1:length(ExonMedian)){													#tests median coverage for each exon; if below 100x, adds to list of fails
	if(ExonMedian[i]<cov_thresh){
		
		Exon<-c(Exon,i)
		Sample<-c(Sample,"All")
		Types<-c(Types,"Whole exon")
		Details<-c(Details,paste("Low median read depth (FPKM): ",ExonMedian[i],sep=""))
		Gene<-c(Gene,paste(bed.file[i,4]))
		
	}
	
}

#N_all<-length(Types)

#for(i in 1:nrow(ReadDepths)){
	
#	if(!i%in%Exon[1:N_all]){
		
#		for(j in 1:ncol(ReadDepths)){
		
#		if(!sample.names[j]%in%Sample[1:N_all]){
			
#			if(ReadDepths[i,j]<100){
					
#					Sample<-c(Sample,sample.names[j])
#					Exon<-c(Exon,i)
#					Types<-c(Types,"Single")
#					Details<-c(Details,paste("Low read depth (FPKM): ",ReadDepths[i,j],sep=""))
#					Gene<-c(Gene,paste(bed.file[i,4]))
					
#			}
			
#			}
		
#		}
	
#	}
	
#}


if(!is.null(exon_numbers)&any(Exon!="All")){

    exons<-read.table(exon_numbers,sep="\t",header=T)
    
    failed.calls=bed.file[Exon[Exon!="All"],]

    brca<-sapply(failed.calls$chr, '==',exons$Chr) & sapply(failed.calls$start, '>=',exons$Start) & sapply(failed.calls$end, '<=',exons$End) | sapply(failed.calls$chr, '==',exons$Chr) & sapply(failed.calls$start, '<=',exons$Start) & sapply(failed.calls$end, '>=',exons$Start) | sapply(failed.calls$chr, '==',exons$Chr) & sapply(failed.calls$start, '<=',exons$End) & sapply(failed.calls$end, '>=',exons$End)

    Brca<-which(colSums(brca)!=0)

    a<-apply(brca[,Brca,drop=F],2,which)

    Clinical_Numbering=rep("NA",length(Exon))

    Clinical_Numbering[paste(Exon)%in%row.names(failed.calls[Brca,])] = exons$Custom[a]


Failures<-data.frame(Sample,Exon,Types,Gene,Clinical_Numbering,Details)

names(Failures)=c("Sample","Exon","Type","Gene","Custom.numbering","Info")

if(Custom){
    
    Failures_custom=Failures[Failures$Clinical_Numbering!="NA" | Failures$Types=="Whole Sample",]
    if(nrow(Failures_custom)>0){
	write.table(Failures_custom,file=paste(OutFails,"_custom_Failures.txt",sep=""),quote=F,row.names=F,sep="\t")
    }
}

}else{
Failures<-data.frame(Sample,Exon,Types,Gene,Details)
names(Failures)=c("Sample","Exon","Type","Gene","Info")
}

if(nrow(Failures)>0){
	write.table(Failures,file=paste(OutFails,"_Failures.txt",sep=""),quote=F,row.names=F,sep="\t")
}




if(BRCA){

brca12=which(bed.file[,4]=="BRCA1" | bed.file[,4]=="BRCA2")

Failures_b1b2<-Failures[Failures$Type=="Whole Sample" | Failures$Exon%in%brca12,]

if(nrow(Failures_b1b2)>0){
        write.table(Failures_b1b2,file=paste(OutFails,"_b1b2_Failures.txt",sep=""),quote=F,row.names=F,sep="\t")
}


}

print("END IdentifyFailures.R")
