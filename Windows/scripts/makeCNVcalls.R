renv::restore()

library(R.utils)

print("BEGIN makeCNVCalls.R")



######Parsing input options and setting defaults########

args=commandArgs(TRUE)

count_data=args[1]
if(count_data=="NULL"){count_data=NULL}
if(is.null(count_data)){
    print("ERROR: no RData summary file provided -- Execution halted")
    quit()
}
load(count_data)

															#R workspace with the coverage data saved in - this workspace already has the bedfile, fasta file saved in it.
trans_prob=as.numeric(args[2])													#transition probability for the HMM
if(length(trans_prob)==0){trans_prob=0.01}

exon_numbers=args[3]
if(exon_numbers=="NULL"){exon_numbers=NULL}

Custom=as.logical(args[4])
if(length(Custom)==0){Custom=FALSE}

b1b2Calls=as.logical(args[5])
if(length(b1b2Calls)==0){b1b2Calls=FALSE}

plotOutput=args[6]
if(plotOutput=="NULL"){plotOutput=NULL}
if(is.null(plotOutput)){plotOutput="All"}
if(!plotOutput%in%c("None","Custom","All")){print("WARNING: plot argument should be one of: None, Custom or All")}

plotFolder=args[7]
if(plotFolder=="NULL"){plotFolder=NULL}
if(is.null(plotFolder)){plotFolder="DECoNPlots"}
if(!file.exists(plotFolder)){dir.create(plotFolder)}

out=args[8]
if(out=="NULL"){out=NULL}
if(is.null(out)){out="DECoNResults"}


#################################################################
library(grid)
library(ExomeDepth)
library(reshape)
library(ggplot2)
#######Make sure bed file is in chromosome order ################

temp<-gsub('chr','',bed.file[,1])
temp1<-order(as.numeric(temp))
bed.file=bed.file[temp1,]


#################CNV CALLING#######################################

ExomeCount<-as(counts, 'data.frame')											#converts counts, a ranged data object, to a data frame

ExomeCount$chromosome <- gsub(as.character(ExomeCount$chromosome),pattern = 'chr',replacement = '') 
																				#remove any chr letters, and coerce to a string.
colnames(ExomeCount)[1:length(sample.names)+4]=sample.names						#assigns the sample names to each column 

cnv.calls = NULL
refs<-list()
models<-list()

for(i in 1:length(sample.names)){												#for each sample:
    print(paste("Processing sample: ",sample.names[i]," ", i,"/",length(sample.names),sep=""))

    my.test <- ExomeCount[,sample.names[i]]										#extracts the sample to be tested
    my.ref.samples <- sample.names[-i]											#uses all other samples as the reference set
    my.reference.set <- as.matrix(ExomeCount[,my.ref.samples])					#places coverage info for all samples in the reference set into a matrix
    my.choice <- select.reference.set (test.counts = my.test,reference.counts = my.reference.set,bin.length = (ExomeCount$end - ExomeCount$start)/1000,n.bins.reduced = 10000)	#from the reference set, selects correlated samples to be used.
    my.matrix <- as.matrix( ExomeCount[, my.choice$reference.choice, drop = FALSE])             #places selected, correlated samples into a matrix
    my.reference.selected <- apply(X = my.matrix,MAR = 1,FUN = sum)				#sums the selected samples across each exon

    all.exons <- new('ExomeDepth', test = my.test, reference = my.reference.selected, formula = 'cbind(test, reference) ~ 1')       #creates ExomeDepth object containing test data, reference data, and linear relationship between them. Automatically calculates likelihoods
    all.exons <- CallCNVs(x = all.exons, transition.probability = trans_prob, chromosome = ExomeCount$chromosome, start = ExomeCount$start, end = ExomeCount$end, name = ExomeCount$exon)	#fits a HMM with 3 states to read depth data; transition.probability - transition probability for HMM from normal to del/dup. Returns ExomeDepth object with CNVcalls

    my.ref.counts <- apply(my.matrix, MAR = 1, FUN = sum)
    
    if(nrow(all.exons@CNV.calls)>0){
        cnvs= cbind(sample.names[i],cor(my.test, my.ref.counts),length( my.choice[[1]] ),all.exons@CNV.calls)
        cnv.calls = rbind(cnv.calls,cnvs)
    }
    
    refs[[i]] = my.choice$reference.choice	
    models[[i]] = c(all.exons@expected[1],all.exons@phi[1])
}


names(refs) = sample.names
names(models) = sample.names
if(!is.null(cnv.calls)){
names(cnv.calls)[1] = "sample"
cnv.calls$sample = paste(cnv.calls$sample)
names(cnv.calls)[2] = "correlation"
names(cnv.calls)[3] = "N.comp"

#####################CONFIDENCE CALCULATION##############################

Flag<-rep("HIGH",nrow(cnv.calls))

a<-which(cnv.calls$correlation<.985)
if(length(a)>0){Flag[a]="LOW"}

a<-which(cnv.calls$reads.ratio<1.25 & cnv.calls$reads.ratio>.75)
if(length(a)>0){Flag[a]="LOW"}

a<-which(cnv.calls$N.comp<=3)
if(length(a)>0){Flag[a]="LOW"}

if(ncol(bed.file)>=4){
    genes<-apply(cnv.calls,1,function(x)paste(unique(bed.file[x[4]:x[5],4]),collapse=", "))
    a<-which(genes=="PMS2")
    if(length(a)>0){Flag[a]="LOW"}
}

cnv.calls<-cbind(cnv.calls,genes)
#names(cnv.calls)[ncol(cnv.calls)]="Confidence"
names(cnv.calls)[ncol(cnv.calls)]="Gene"

###################FORMAT CALLS######################################

cnv.calls_ids=cbind(1:nrow(cnv.calls),cnv.calls)
names(cnv.calls_ids)[1]="ID"

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

for(i in 1:nrow(cnv.calls_ids)){                        #replaces single calls involving multiple genes with multiple calls with single call ID
    genes=strsplit(paste(cnv.calls_ids[i,]$Gene),",")[[1]]
    genes=trim(genes)
    whole.index=cnv.calls_ids[i,]$start.p:cnv.calls_ids[i,]$end.p
    if(length(genes)>1){
        temp=cnv.calls_ids[rep(i,length(genes)),]
        temp$Gene=genes
        for(j in 1:length(genes)){
            gene.index=which(bed.file[,4]==genes[j])
            overlap=gene.index[gene.index%in%whole.index]
            temp[j,]$start.p=min(overlap)
            temp[j,]$end.p=max(overlap)
        }
        if(i==1){
            cnv.calls_ids=rbind(temp,cnv.calls_ids[(i+1):nrow(cnv.calls_ids),])
        }else if(i==nrow(cnv.calls_ids)){
            cnv.calls_ids=rbind(cnv.calls_ids[1:(i-1),],temp)
        }else{
            cnv.calls_ids=rbind(cnv.calls_ids[1:(i-1),],temp,cnv.calls_ids[(i+1):nrow(cnv.calls_ids),])
        }
    }
}
cnv.calls_ids$Gene=trim(cnv.calls_ids$Gene)


Gene.index=vector()
genes_unique<-unique(bed.file[,4])
for(i in 1:length(genes_unique)){

    Gene.index=c(Gene.index,1:sum(bed.file[,4]==genes_unique[i]))

}

start.b<-Gene.index[cnv.calls_ids$start.p]
end.b<-Gene.index[cnv.calls_ids$end.p]

cnv.calls_ids<-cbind(cnv.calls_ids,start.b,end.b)


#################ADDS CLINICAL EXON NUMBERS#####################

if(!is.null(exon_numbers)){
    Clinical.first = rep(NA,nrow(cnv.calls_ids))
    Clinical.last=rep(NA,nrow(cnv.calls_ids))
    exons<-read.table(exon_numbers,sep="\t",header=T)
    
    brca<-sapply(cnv.calls_ids$Gene, '==',exons$Gene)
    for(i in 1:nrow(brca)){
        for(j in 1:ncol(brca)){
            temp=cnv.calls_ids$start[j]<=exons$End[i] & cnv.calls_ids$end[j]>=exons$Start[i]
            brca[i,j]=brca[i,j] & temp
        }
    }
    Brca<-which(colSums(brca)!=0)
    if(length(Brca)>0){
        a<-list(length=length(Brca))
        for(i in 1:length(Brca)){a[[i]]=which(brca[,Brca[i]])}
        first_exon<-unlist(lapply(a,function(a,b)min(b[a,]$Custom),exons))		#identifies the  first and last clinical exon in the deletion/duplication.
        last_exon<-unlist(lapply(a,function(a,b)max(b[a,]$Custom),exons))
        Clinical.first[Brca] = first_exon
        Clinical.last[Brca] = last_exon
    }
    cnv.calls_ids=cbind(cnv.calls_ids,Clinical.first,Clinical.last)
}

###############################################################

if(!is.null(exon_numbers)){
cnv.calls_ids_out<-data.frame(cnv.calls_ids$ID,cnv.calls_ids$sample,cnv.calls_ids$correlation,cnv.calls_ids$N.comp,cnv.calls_ids$start.p,cnv.calls_ids$end.p,cnv.calls_ids$type,cnv.calls_ids$nexons,cnv.calls_ids$start,cnv.calls_ids$end,cnv.calls_ids$chromosome,cnv.calls_ids$id,cnv.calls_ids$BF,cnv.calls_ids$reads.expected,cnv.calls_ids$reads.observed,cnv.calls_ids$reads.ratio,cnv.calls_ids$Gene,cnv.calls_ids$Clinical.first,cnv.calls_ids$Clinical.last)
names(cnv.calls_ids_out)<-c("CNV.ID","Sample","Correlation","N.comp","Start.b","End.b","CNV.type","N.exons","Start","End","Chromosome","Genomic.ID","BF","Reads.expected","Reads.observed","Reads.ratio","Gene","Custom.first","Custom.last")
}else{
cnv.calls_ids_out<-data.frame(cnv.calls_ids$ID,cnv.calls_ids$sample,cnv.calls_ids$correlation,cnv.calls_ids$N.comp,cnv.calls_ids$start.p,cnv.calls_ids$end.p,cnv.calls_ids$type,cnv.calls_ids$nexons,cnv.calls_ids$start,cnv.calls_ids$end,cnv.calls_ids$chromosome,cnv.calls_ids$id,cnv.calls_ids$BF,cnv.calls_ids$reads.expected,cnv.calls_ids$reads.observed,cnv.calls_ids$reads.ratio,cnv.calls_ids$Gene)
names(cnv.calls_ids_out)<-c("CNV.ID","Sample","Correlation","N.comp","Start.b","End.b","CNV.type","N.exons","Start","End","Chromosome","Genomic.ID","BF","Reads.expected","Reads.observed","Reads.ratio","Gene")

}

write.table(cnv.calls_ids_out,file=paste(out,"_all.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=T)


################BRCA OUTPUT###################################

if(b1b2Calls){
    BRCA = bed.file[bed.file[,4]=="BRCA1"|bed.file[,4]=="BRCA2" ,]
    names(BRCA)=c("Chr","Start","End","Gene") 
    brca<-sapply(cnv.calls_ids$chr, '==',BRCA$Chr) & sapply(cnv.calls_ids$start, '>=',BRCA$Start) & sapply(cnv.calls_ids$end, '<=',BRCA$End) | sapply(cnv.calls_ids$chr, '==',BRCA$Chr) & sapply(cnv.calls_ids$start, '<=',BRCA$Start) & sapply(cnv.calls_ids$end, '>=',BRCA$Start) | sapply(cnv.calls_ids$chr, '==',BRCA$Chr) & sapply(cnv.calls_ids$start, '<=',BRCA$End) & sapply(cnv.calls_ids$end, '>=',BRCA$End)
    Brca<-which(colSums(brca)!=0)
    a<-apply(brca[,Brca,drop=F],2,which)		
    b1b2_calls<-cnv.calls_ids_out[Brca,]
    write.table(b1b2_calls,file=paste(out,"_b1b2.txt",sep=""),sep="\t",quote=F,col.names=TRUE,row.names=FALSE)
}


####################CUSTOM OUTPUT#############################


if(Custom){

    custom_calls=cnv.calls_ids_out[!is.na(cnv.calls_ids$Clinical.first),]
    write.table(custom_calls,file=paste(out,"_custom.txt",sep=""),sep="\t",quote=F,col.names=TRUE,row.names=FALSE)

}



####################Plotting####################################

if(plotOutput!="None"){

    if(plotOutput=="Clinical"){
        cnv.calls_plot=cnv.calls_ids[!is.na(cnv.calls_ids$Clinical.first),]
    }else{
        cnv.calls_plot=cnv.calls_ids
    }

	cnv.calls_plot$chr=paste('chr',cnv.calls_plot$chromosome,sep='')

    Index=vector(length=nrow(bed.file))
    Index[1]=1
    for(i in 2:nrow(bed.file)){
        if(bed.file[i,4]==bed.file[i-1,4]){
            Index[i]=Index[i-1]+1
        }else{
            Index[i]=1
        }   
    }   

    if(!is.null(exon_numbers)){

        for(i in 1:nrow(exons)){
            x=which(paste(bed.file[,4])==paste(exons[i,4]) & bed.file[,2]<=exons[i,3] & bed.file[,3]>=exons[i,2])
            Index[x]=exons[i,5]
        }
    }


    
    for(call_index in 1:nrow(cnv.calls_plot)){
        
        Sample<-cnv.calls_plot[call_index,]$sample
        Gene<-unlist(strsplit(cnv.calls_plot[call_index,]$Gene,split=", "))
        exonRange<-which(bed.file[,4]%in%Gene)

        if((cnv.calls_plot[call_index,]$start.p-5)<min(exonRange) & ((cnv.calls_plot[call_index,]$start.p-5)>=1)){exonRange=(cnv.calls_plot[call_index,]$start.p-5):max(exonRange)}
        if((cnv.calls_plot[call_index,]$start.p-5)<min(exonRange) & ((cnv.calls_plot[call_index,]$start.p-5)<=0)){exonRange=1:max(exonRange)}
        if((cnv.calls_plot[call_index,]$end.p+5)>max(exonRange) & ((cnv.calls_plot[call_index,]$end.p+5)<=nrow(bed.file))){exonRange=min(exonRange):(cnv.calls_plot[call_index,]$end.p+5)}
        if((cnv.calls_plot[call_index,]$end.p+5)>max(exonRange) & ((cnv.calls_plot[call_index,]$end.p+5)>nrow(bed.file))){exonRange=min(exonRange):nrow(bed.file)}

        singlechr=length(unique(bed.file[exonRange,1]))==1

        if(!singlechr){
            if(bed.file[exonRange[1],1]!=cnv.calls_plot[call_index,]$chr){
                prev=TRUE
                newchr=bed.file[exonRange[1],1]
            }else{
                prev=FALSE
                newchr=bed.file[exonRange[length(exonRange)],1]
            }
            exonRange=exonRange[bed.file[exonRange,1]==cnv.calls_plot[call_index,]$chr]
        }

        ######First Plot###############
        VariantExon<- unlist(mapply(function(x,y)x:y,cnv.calls[cnv.calls$sample==Sample,]$start.p,cnv.calls[cnv.calls$sample==Sample,]$end.p))
        refs_sample<-refs[[Sample]]
        
        Data<-cbind(ExomeCount[exonRange,c(Sample,refs_sample)],exonRange)
        Data[,-ncol(Data)]=log(Data[,-ncol(Data)])
        Data1<-melt(Data,id=c("exonRange"))

        testref<-rep("gray",nrow(Data1))
        testref[Data1$variable==Sample]="blue"
        Data1<-data.frame(Data1,testref)
        levels(Data1$variable)=c(levels(Data1$variable),"VAR")
	Data1$testref=as.factor(Data1$testref)
        levels(Data1$testref)=c(levels(Data1$testref),"red")

        data_temp<-Data1[Data1$variable==Sample & Data1$exonRange%in%VariantExon,]

        if(nrow(data_temp)>0){
            data_temp$variable="VAR"
            data_temp$testref="red"
            Data1<-rbind(Data1,data_temp)
        }   
        levels(Data1$testref)=c("Test Sample","Reference Sample","Affected exon")

        new_cols=c("blue","gray","red")

        A1<-ggplot(data=Data1,aes(x=exonRange,y=value,group=variable,colour=testref))
        A1<-A1 + geom_point(cex=2.5,lwd=1.5)                        #Have to set up points with scale to get correct legend, then re-plot in correct order etc.
        A1<-A1 + scale_colour_manual(values=new_cols)  
        A1<-A1 + geom_line(data=subset(Data1,testref=="Reference Sample"),lty="dashed",lwd=1.5,col="grey") 
        A1<-A1 + geom_point(data=subset(Data1,testref=="Reference Sample"),cex=2.5,col="grey")   
        A1<-A1+ geom_line(data=subset(Data1,testref=="Test Sample"),lty="dashed",lwd=1.5,col="blue")  
        A1<-A1 + geom_point(data=subset(Data1,testref=="Test Sample"),cex=2.5,col="blue") 
        A1<-A1 + geom_point(data=subset(Data1,testref=="Affected exon"),cex=3.5,col="red") 
        A1<-A1 + ylab("Log (Coverage)")  + theme_bw() + theme(legend.position="none")+ xlab(" ")

        Data2<-Data1[Data1$testref=="Affected exon",]
        if(nrow(Data2)>1){
            for(i in 1:(nrow(Data2)-1)){
                if((Data2$exonRange[i]+1)==Data2$exonRange[i+1]){A1<-A1 + geom_line(data=Data2[i:(i+1),],aes(x=exonRange,y=value,group=1),lwd=1.5,col="red")}
            }
        }   


        if(!singlechr){
            if(prev){A1<-A1 + scale_x_continuous(breaks=(min(exonRange)-6):max(exonRange),labels=c(rep("",6),paste(Index[exonRange])),limits=c(min(exonRange)-6.75,max(exonRange)))
            }else{A1<-A1 + scale_x_continuous(breaks=min(exonRange):(max(exonRange)+6),labels=c(paste(Index[exonRange]),rep("",6)),limits=c(min(exonRange),max(exonRange)+6.75))}
        }else{A1<-A1 + scale_x_continuous(breaks=exonRange,labels=paste(Index[exonRange]))}


        ##############SECOND PLOT###########

        genes_sel = unique(bed.file[exonRange,4])
        temp<-cbind(1:nrow(bed.file),bed.file)[exonRange,]
        len<-table(temp$name)
        mp<-tapply(exonRange,temp[,5],mean)
        mp<-mp[genes_sel]
        len<-len[genes_sel]
        Genes<-data.frame(c(genes_sel),c(mp),c(len-.5),1)
        names(Genes)=c("Gene","MP","Length","Ind")

        if(!singlechr){
            levels(Genes$Gene)=c(levels(Genes$Gene),newchr)
            if(prev){
                fakegene=c(newchr,min(exonRange)-5,3.5,1)
            }else{
                fakegene=c(newchr,max(exonRange)+5,3.5,1)
            }
            Genes<-rbind(Genes,fakegene)
            Genes$MP=as.numeric(Genes$MP)
            Genes$Length=as.numeric(Genes$Length)
        }

        GenesPlot<-ggplot(data=Genes, aes(x=MP,y=Ind,fill=Gene,width=Length,label=Gene)) +geom_tile() + geom_text() + theme_bw() + theme(legend.position="None",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),plot.margin=unit(c(.5,.5,.5,.55),"cm")) + ylab(" ") + xlab(" ")
        GenesPlot<-GenesPlot + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

        #######THIRD PLOT#############

        Totals<-rowSums(ExomeCount[exonRange,c(Sample,refs_sample)])
        ratio = (ExomeCount[exonRange,Sample]/Totals)/models[[Sample]][1]
        mins <- vector(length=length(exonRange))
        maxs <-vector(length=length(exonRange))
        for(i in 1:length(exonRange)){
            temp = qbetabinom(p=0.025,Totals[i],models[[Sample]][2],models[[Sample]][1])
            mins[i] = (temp/Totals[i])/models[[Sample]][1]
            temp = qbetabinom(p=0.975,Totals[i],models[[Sample]][2],models[[Sample]][1])
            maxs[i] = (temp/Totals[i])/models[[Sample]][1]
        }       

        CIData<-data.frame(exonRange,ratio,mins,maxs)
        names(CIData)<-c("Exon","Ratio","Min","Max")
        CIPlot<-ggplot(CIData,aes(x=Exon,y=Ratio))+geom_ribbon(aes(ymin=Min,ymax=Max),fill="grey")+geom_point(col="blue",cex=3.5) + theme_bw() +xlab("")+ylab("Observed/Expected")
temp = cnv.calls[cnv.calls$sample==Sample,]
        if(sum(temp$start.p%in%exonRange |temp$end.p%in%exonRange)>0){
            temp = temp[temp$start.p%in%exonRange|temp$end.p%in%exonRange,]
            for(i in 1:nrow(temp)){
                start.temp = temp[i,]$start.p
                end.temp = temp[i,]$end.p
                CIPlot<-CIPlot + geom_point(data=CIData[CIData$Exon%in%start.temp:end.temp,], aes(x=Exon,y=Ratio),color="red",cex=3.5)
            }
        }

        if(!singlechr){
            if(prev){CIPlot<- CIPlot + scale_x_continuous(breaks=(min(exonRange)-6):max(exonRange),labels=c(rep("",6),paste(Index[exonRange])),limits=c(min(exonRange)-6.75,max(exonRange)))
            }else{CIPlot<-CIPlot + scale_x_continuous(breaks=min(exonRange):(max(exonRange)+6),labels=c(paste(Index[exonRange]),rep("",6)),limits=c(min(exonRange),max(exonRange)+6.75))}
        }else{CIPlot<-CIPlot + scale_x_continuous(breaks=exonRange,labels=paste(Index[exonRange]))}

        #########WRTITE PLOTS TO FILE###########
    
        cnv_genes_sample=cnv.calls_plot[cnv.calls_plot$sample==Sample,]$Gene
        if(sum(cnv_genes_sample==Gene)==1){
            pdf(file=paste(plotFolder,"/",Sample,"_",Gene,".pdf",sep=""),useDingbats=FALSE)
        }else{
            cnv_genes_sample_index=which(cnv.calls_plot$sample==Sample & cnv.calls_plot$Gene==paste(Gene,collapse=", "))
            pdf(file=paste(plotFolder,"/",Sample,"_",paste(Gene,collapse="_"),"_",which(cnv_genes_sample_index==call_index),".pdf",sep=""),useDingbats=F)
        }
    
        if(sum(cnv_genes_sample==Gene)>2){print(paste("WARNING: more than 2 calls in ",Gene,", could affect plotting",sep=""))}

        grid.newpage()
        pushViewport(viewport(layout = grid.layout(6, 1)))

        print(A1,vp=viewport(layout.pos.row=1:3,layout.pos.col=1))
        print(GenesPlot,vp=viewport(layout.pos.row=4,layout.pos.col=1))
        print(CIPlot,vp=viewport(layout.pos.row=5:6,layout.pos.col=1))

        dev.off()

    }

}



##############################################################

if(!is.null(exon_numbers)){
save(ExomeCount,bed.file,counts,fasta,sample.names,bams,cnv.calls,cnv.calls_ids,refs,models,exon_numbers,exons,file=paste(out,".RData",sep=""))
}else{
save(ExomeCount,bed.file,counts,fasta,sample.names,bams,cnv.calls,cnv.calls_ids,refs,models,exon_numbers,file=paste(out,".RData",sep=""))
}
}
print("END makeCNVCalls.R")
