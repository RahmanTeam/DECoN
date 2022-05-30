library(shiny)
library(ggplot2)
library(reshape)
library(grid)
library(ExomeDepth)


load('Data.RData')          #wrapper script copies data to this hard path


# Define server logic required to draw a histogram
shinyServer(function(input, output) {



    bed.file.disp<-cbind(1:nrow(bed.file),bed.file)
    if(ncol(bed.file.disp==4)){colnames(bed.file.disp)=c("Exon Index","Chromosome","Start","Stop")}
    if(ncol(bed.file.disp==5)){colnames(bed.file.disp)=c("Exon Index","Chromosome","Start","Stop","Gene")}
   # bed.file.disp=cbind(bed.file.disp,ExomeCount$names)
   # colnames(bed.file.disp)[ncol(bed.file.disp)]="Gene.Exon"	
    names_files<-data.frame(sample.names,bams)
    colnames(names_files)<-c("Sample Name","File")


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

    bed.file.disp=cbind(bed.file.disp,Index)
    colnames(bed.file.disp)[ncol(bed.file.disp)]="Custom.exon"
    
#####Data tab
    output$bed.file<-renderDataTable({
	if(input$readBed%%2==0){bed.file.disp}
    })

    output$bamlist <- renderDataTable({
        if(input$readBams%%2==0){names_files}
    })

    output$fasta <- renderPrint({
	if(input$readfasta%%2==0){fasta}
    })

#####Coverage Evaluation tab
    output$downloadCov<-downloadHandler(
        filename=function(){
	    paste('Coverage-',Sys.Date(),'.txt',sep='')
	},
	content = function(file){
	    write.table(ExomeCount[,sample.names],file,sep="\t",quote=F,row.names=F)
	}
    )


    output$FailedSamples<-renderDataTable({         
            
        ReadDepths<-ExomeCount[,sample.names]                                           #extracts just the read depths
        Sample<-vector()
        Details<-vector()
        Correlation<-vector()
        MedCov<-vector()
        Corr<-cor(ReadDepths)                                                           #calculates correlation matrix
        MaxCorr<-apply(Corr,1,function(x)max(x[x!=1]))                              #finds the maximum correlation for each sample
        SampleMedian<-apply(ReadDepths,2,median)
        
        for(i in 1:length(MaxCorr)){                                                    #tests correlation and median coverage for each sample; if below threshold, adds to list of fails and adds reason to Details list
            if(MaxCorr[i]<input$SampleCorLimit[1]|SampleMedian[i]<input$SampleCovLimit[1]){
                Sample<-c(Sample,sample.names[i])
                Correlation<-c(Correlation,MaxCorr[i])
                MedCov<-c(MedCov,SampleMedian[i])
                if(MaxCorr[i]<input$SampleCorLimit[1] ){
                    if(SampleMedian[i]<input$SampleCovLimit[1]){
                        Details<-c(Details,"Low correlation and low median coverage ")
                    }else{Details<-c(Details,"Low correlation")}
                }else{
                    Details<-c(Details,"Low median coverage")
                }
            }
        }

        Failures<-data.frame(Sample,Details,Correlation,MedCov)
        names(Failures)[4]="Median Coverage"

        Failures

        })



        output$FailedExons<-renderDataTable({
            
            ReadDepths<-ExomeCount[,sample.names] 
            ExonMedian<-apply(ReadDepths,1,median)

            Exon<-which(ExonMedian<input$ExonCovLimit[1])
            Custom.exon<-Index[Exon]
            Gene<-bed.file[Exon,4]
            Coverage<-ExonMedian[Exon]

            Failures<-data.frame(Exon,Gene,Custom.exon,Coverage)
            names(Failures)=c("Exon number (bed file)","Gene","Exon number (custom)","Coverage")
            Failures
        })



        output$PlotSampleHighlight<-renderUI({
               selectInput("SampHigh",choices=as.list(c(sample.names)),label="Select sample to highlight",multiple=F,selected=sample.names[1])        
        })
 


        output$PlotSamplesInput<-renderUI({
            if(input$ChooseFrom==1){    
                selectInput("PlotSamp",choices=as.list(refs[[input$SampHigh]]),label="Select samples to display",multiple=T,selected=paste(refs[[input$SampHigh]]))        
            }else if(input$ChooseFrom==2){    
                selectInput("PlotSamp",choices=as.list(sample.names),label="Select samples to display",multiple=T,selected=paste(sample.names))        
            }
        })


        output$PlotGenesInput<-renderUI({
                    selectInput("PlotGenes",choices=as.list(paste(unique(bed.file[,4]))),label="Select genes to plot",multiple=T,selected=paste(unique(bed.file[,4])[1:2]))
        })



        output$CovPlot<-renderPlot({

            if(input$plotType==1){

                exons<-which(bed.file[,4]%in%input$PlotGenes)
                Data<-cbind(exons,ExomeCount[exons,sample.names])
                Data1<-melt(Data,id=c("exons"))
                Data1$exons<-as.factor(Data1$exons)
                if(input$plotScale==2){Data1$value=log(Data1$value)}
                p<-ggplot(data=Data1[Data1$variable%in%input$PlotSamp,],aes(x=exons,y=value))  + geom_boxplot(width=0.75) + theme_bw() + xlab(NULL)+ ggtitle("")
                if(input$plotScale==1){p <- p + ylab("Coverage")}
                if(input$plotScale==2){p <- p + ylab("Log (Coverage)")}
                p<-p + geom_point(data=Data1[Data1$variable==input$SampHigh,],aes(x=exons,y=value),colour="blue",cex=2.5)
                p
            
            }else if(input$plotType==2){

                exons<-which(bed.file[,4]%in%input$PlotGenes)
                Data<-cbind(exons,ExomeCount[exons,sample.names])
                Data1<-melt(Data,id=c("exons"))
                Data1$exons<-as.factor(Data1$exons)
                if(input$plotScale==2){Data1$value=log(Data1$value)}
                p<-ggplot(data=Data1[Data1$variable%in%input$PlotSamp,],aes(x=exons,y=value))  + geom_point(width=0.75) + theme_bw()  + xlab(NULL)+ ggtitle("")
                if(input$plotScale==1){p <- p + ylab("Coverage")}
                if(input$plotScale==2){p <- p + ylab("Log (Coverage)")}
                for(gene in input$PlotGenes){
                    p<-p + geom_line(data=Data1[Data1$variable%in%input$PlotSamp & Data1$exons%in%which(bed.file[,4]==gene),],aes(x=exons,y=value,group=variable),colour="grey")
                }
                p<-p + geom_point(data=Data1[Data1$variable==input$SampHigh,],aes(x=exons,y=value),colour="blue",cex=2.5)
                for(gene in input$PlotGenes){
                    p<-p + geom_line(data=Data1[Data1$variable==input$SampHigh & Data1$exons%in%which(bed.file[,4]==gene),],aes(x=exons,y=value,group=variable),colour="blue")
                }
                p
            }

        })

	



#####CNV calls tab

    trim <- function (x) gsub("^\\s+|\\s+$", "", x)

    if(!exists("cnv.calls_ids") | !exists("cnv.calls_ids$Confidence")){         #Formats CNV calls if generated by old DECoN version
        if(!exists("cnv.calls$Gene")){
            cnv.calls$Gene=apply(cnv.calls,1,function(x)paste(unique(bed.file[x[4]:x[5],4]),collapse=", "))         
        }

        if(!exists("cnv.calls$Confidence")){
            Flag<-rep("HIGH",nrow(cnv.calls))

            a<-which(cnv.calls$correlation<.985)
            if(length(a)>0){Flag[a]="LOW"}

            a<-which(cnv.calls$reads.ratio<1.25 & cnv.calls$reads.ratio>.75)
            if(length(a)>0){Flag[a]="LOW"}

            Ncomp<-lapply(refs,length)
            a<-which(Ncomp<=3)
            if(length(a)>0){
                samples_a<-sample.names[a]
                b<-which(cnv.calls$sample%in%samples_a)
                Flag[b]="LOW"
            }

            a<-which(cnv.calls$Gene=="PMS2")
            if(length(a)>0){Flag[a]="LOW"}

            cnv.calls$Confidence=Flag
        }


        cnv.calls_ids=cbind(1:nrow(cnv.calls),cnv.calls)
        names(cnv.calls_ids)[1]="ID"


        for(i in 1:nrow(cnv.calls_ids)){

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

    }    

    output$downloadCNVs<-downloadHandler(
	filename=function(){
	        paste('CNVs-',Sys.Date(),'.txt',sep='')
	},
	content = function(file){
	    write.table(cnv.calls_ids,file,sep="\t",quote=F,row.names=F)
	}
    )

	
    genes.exons.short<-apply(cnv.calls,1,function(x)paste(ExomeCount$names[x[4]:x[5]],collapse=", "))	
    genes.exons.short<-strsplit(genes.exons.short,",")
    genes.exons.short<-lapply(genes.exons.short,trim)
    genes.exons.split<-vector(length=nrow(cnv.calls_ids))
    
    for(i in 1:nrow(cnv.calls_ids)){
        j=cnv.calls_ids$ID[i]
        single.gene=cnv.calls_ids$Gene[i]
        all.genes=unlist(lapply(strsplit(genes.exons.short[[j]],".",fixed=T),function(x)x[1]))
        genes.exons.split[i]=paste(genes.exons.short[[j]][single.gene==all.genes],collapse=", ")
    }
        
   # CNVs<-cbind(cnv.calls_ids[,c("ID","sample","start.p","end.p","nexons","Gene")],genes.exons.split,cnv.calls_ids[,c("type","reads.ratio","correlation","Confidence")])
   # names(CNVs)<-c("CNV identifier","Sample","First exon","Last exon","Number of exons","Gene","Gene.Exon","Type","Read ratio","Correlation","Confidence")
    

   # CNVs<-cbind(cnv.calls_ids[,c("ID","sample","start.p","end.p","nexons","Gene")],genes.exons.split,cnv.calls_ids[,c("type","reads.ratio","correlation")])
   # names(CNVs)<-c("CNV ID","Sample","First exon","Last exon","Number of exons","Gene","Gene.Exon","Type","Read ratio","Correlation")
    
    CNVs<-cbind(cnv.calls_ids[,c("ID","sample","start.p","end.p","nexons","Gene")],Index[cnv.calls_ids$start.p],Index[cnv.calls_ids$end.p],cnv.calls_ids[,c("type","reads.ratio","correlation")])
    names(CNVs)<-c("CNV ID","Sample","First exon (BED file)","Last exon (BED file)","Number of exons","Gene","First exon (custom)","Last exon (custom)","Type","Read ratio","Correlation")
    



    output$CNVcalls<-renderDataTable(
        CNVs,options=list(iDisplayLength=10)
    )

    output$selVar<-renderUI({
        selectInput("selVar1","Use the CNV ID given in the table above",c('None',1:nrow(cnv.calls)))
    })

    output$minEx <- renderUI({
        if(input$selVar1!="None"){
            numericInput("minEx1",min=1,max=nrow(bed.file)-1,value=max(cnv.calls[input$selVar1,]$start.p-5,1),label="First exon")
        }
    })

    output$maxEx <- renderUI({
        if(input$selVar1!="None"){
            numericInput("maxEx1",min=2,max=nrow(bed.file),value=min(cnv.calls[input$selVar1,]$end.p+5,nrow(bed.file)),label="Last")
        }
    })

    output$plot<-renderPlot({
        if(input$selVar1=="None"){
            plot(NULL,xlim=c(1,10),ylim=c(0,1000))
        }else{
            Sample<-cnv.calls[input$selVar1,]$sample
            exonRange<-input$minEx1:input$maxEx1
	    VariantExon<- unlist(mapply(function(x,y)x:y,cnv.calls[cnv.calls$sample==Sample,]$start.p,cnv.calls[cnv.calls$sample==Sample,]$end.p))       
            if(input$chSamp==1){
                refs_sample<-refs[[Sample]]
                Data<-cbind(ExomeCount[exonRange,c(Sample,refs_sample)],exonRange) 
                if(input$chScale==2){Data[,-ncol(Data)]=log(Data[,-ncol(Data)])}
               # if(input$chScale==3){Data_temp<-data.frame(cbind(t(apply(Data[,1:(ncol(Data)-1)],1,function(x)x-median(x))),exonRange));names(Data_temp)<-names(Data);Data<-Data_temp}
               # if(input$chScale==4){Data_temp<-data.frame(cbind(apply(Data[,1:(ncol(Data)-1)],2,function(x)x-median(x)),exonRange));names(Data_temp)<-names(Data);Data<-Data_temp}
               # if(input$chScale==5){Data_temp<-t(apply(Data[,1:(ncol(Data)-1)],1,function(x)x-median(x)));Data_temp<-apply(Data_temp,2,function(x)x-median(x));Data_temp<-data.frame(Data_temp,exonRange);names(Data_temp)<-names(Data);Data<-Data_temp}
                Data1<-melt(Data,id=c("exonRange"))
		testref<-rep("gray",nrow(Data1))
		testref[Data1$variable==Sample]="blue"
		Data1<-data.frame(Data1,testref)
		levels(Data1$variable)=c(levels(Data1$variable),"VAR")
		Data1$testref<-as.factor(Data1$testref)
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
                A1<-A1 + geom_point(cex=2.5,lwd=1.5) 
                A1<-A1 + scale_colour_manual(values=new_cols)  
                A1<-A1 + geom_line(data=subset(Data1,testref=="Reference Sample"),lty="dashed",lwd=1.5,col="grey") 
                A1<-A1 + geom_point(data=subset(Data1,testref=="Reference Sample"),cex=2.5,col="grey")   
                A1<-A1 + geom_line(data=subset(Data1,testref=="Test Sample"),lty="dashed",lwd=1.5,col="blue")  
                A1<-A1 + geom_point(data=subset(Data1,testref=="Test Sample"),cex=2.5,col="blue") 
                A1<-A1 + geom_point(data=subset(Data1,testref=="Affected exon"),cex=3.5,col="red") 
                 if(input$chScale==1){
                    A1<-A1 + ylab("Coverage") + xlab("")
                }
                if(input$chScale==2){
                    A1<-A1 + ylab("Log (Coverage)") + xlab("")
                }
                A1<-A1 + theme_bw() + theme(legend.position="none",axis.text.x=element_blank())
                A1<-A1 + scale_x_continuous(breaks=exonRange)#,labels=paste(exonRange))


                Data2<-Data1[Data1$testref=="Affected exon",]
                if(nrow(Data2)>1){
                    for(i in 1:(nrow(Data2)-1)){
                        if((Data2$exonRange[i]+1)==Data2$exonRange[i+1]){ A1<-A1 + geom_line(data=Data2[i:(i+1),],aes(x=exonRange,y=value,group=1),lwd=1.5,col="red")}
                    } 
                }
                
                print(A1)
        
                }else if(input$chSamp==2){
                    refs_sample<-sample.names[sample.names!=Sample]
                    Data<-cbind(ExomeCount[exonRange,c(Sample,refs_sample)],exonRange)
                    if(input$chScale==2){Data[,-ncol(Data)]=log(Data[,-ncol(Data)])}
                    #if(input$chScale==3){Data_temp<-data.frame(cbind(t(apply(Data[,1:(ncol(Data)-1)],1,function(x)x-median(x))),exonRange));names(Data_temp)<-names(Data);Data<-Data_temp}
                    #if(input$chScale==4){Data_temp<-data.frame(cbind(apply(Data[,1:(ncol(Data)-1)],2,function(x)x-median(x)),exonRange));names(Data_temp)<-names(Data);Data<-Data_temp}
                    #if(input$chScale==5){Data_temp<-t(apply(Data[,1:(ncol(Data)-1)],1,function(x)x-median(x)));Data_temp<-apply(Data_temp,2,function(x)x-median(x));Data_temp<-data.frame(Data_temp,exonRange);names(Data_temp)<-names(Data);Data<-Data_temp}
		    Data1<-melt(Data,id=c("exonRange"))
                    testref<-rep("gray",nrow(Data1))
		    testref[Data1$variable==Sample]="blue"
		    Data1<-data.frame(Data1,testref)
		    levels(Data1$variable)=c(levels(Data1$variable),"VAR")
		    Data1$testref<-as.factor(Data1$testref)
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
                    A1<-A1 + geom_point(cex=2.5,lwd=1.5) 
                    A1<-A1 + scale_colour_manual(values=new_cols)  
                    A1<-A1 + geom_line(data=subset(Data1,testref=="Reference Sample"),lty="dashed",lwd=1.5,col="grey") 
                    A1<-A1 + geom_point(data=subset(Data1,testref=="Reference Sample"),cex=2.5,col="grey")   
                    A1<-A1 + geom_line(data=subset(Data1,testref=="Test Sample"),lty="dashed",lwd=1.5,col="blue")  
                    A1<-A1 + geom_point(data=subset(Data1,testref=="Test Sample"),cex=2.5,col="blue") 
                    A1<-A1 + geom_point(data=subset(Data1,testref=="Affected exon"),cex=3.5,col="red") 
                    if(input$chScale==1){
                        A1<-A1 + ylab("Coverage") + xlab("")
                    }
                    if(input$chScale==2){
                        A1<-A1 + ylab("Log (Coverage)") + xlab("")
                    }   
                    A1<-A1 + theme_bw() + theme(legend.position="none",axis.text.x=element_blank())
                    A1<-A1 + scale_x_continuous(breaks=exonRange)#,labels=paste(exonRange))


                Data2<-Data1[Data1$testref=="Affected exon",]
                if(nrow(Data2)>1){
                    for(i in 1:(nrow(Data2)-1)){
                        if((Data2$exonRange[i]+1)==Data2$exonRange[i+1]){ A1<-A1 + geom_line(data=Data2[i:(i+1),],aes(x=exonRange,y=value,group=1),lwd=1.5,col="red") }
                    } 
                }
                
                print(A1)
        
 

                }else if(input$chSamp==3){
                    Data1<-data.frame(rep(Sample,length(exonRange)),ExomeCount[exonRange,Sample],exonRange)
                    names(Data1)<-c("variable","value","exonRange")       
                    if(input$chScale==2){Data1$value=log(Data1$value)} 
                #    if(input$chScale==4){Data1$value = Data1$value - median(ExomeCount[,Sample])}
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
		    levels(Data1$testref)=c("Test Sample","Affected exon")
		    new_cols=c("blue","red")
                	
                    A1<-ggplot(data=Data1,aes(x=exonRange,y=value,group=variable,colour=testref))
                    A1<-A1 + geom_point(cex=2.5,lwd=1.5) 
                    A1<-A1 + scale_colour_manual(values=new_cols) 
                    A1<-A1 + geom_line(data=subset(Data1,testref=="Test Sample"),lty="dashed",lwd=1.5,col="blue") 
                    A1<-A1 + geom_point(data=subset(Data1,testref=="Affected exon"),cex=3.5,col="red") 
                    if(input$chScale==1){
                        A1<-A1 + ylab("Coverage") + xlab("")
                    }
                    if(input$chScale==2){
                        A1<-A1 + ylab("Log (Coverage)") + xlab("")
                    } 
                    A1<-A1 + theme_bw() + theme(legend.position="none",axis.text.x=element_blank())
                    A1<-A1 + scale_x_continuous(breaks=exonRange)#,labels=paste(exonRange))

                Data2<-Data1[Data1$testref=="Affected exon",]
                if(nrow(Data2)>1){
                for(i in 1:(nrow(Data2)-1)){
                        if((Data2$exonRange[i]+1)==Data2$exonRange[i+1]){ A1<-A1 + geom_line(data=Data2[i:(i+1),],aes(x=exonRange,y=value,group=1),lwd=1.5,col="red") }
                    } 
                }
 
                print(A1)

            }
        }
   
    })


    
    output$genes<-renderPlot({
        if(input$selVar1=="None"){
            par(mar=rep(0,4))
            plot.new()
        }else{
            exonRange<-input$minEx1:input$maxEx1
            genes_sel = unique(bed.file[exonRange,4])
            temp<-cbind(1:nrow(bed.file),bed.file)[exonRange,]
            len<-table(temp$name)
            mp<-tapply(exonRange,temp[,5],mean)
            mp<-mp[genes_sel]
            len<-len[genes_sel]
            Genes<-data.frame(genes_sel,as.vector(mp),as.vector(len-.5),1)
            names(Genes)=c("Gene","MP","Length","Ind")
   
            if(!is.null(exon_numbers)){

             qplot(data=Genes,MP,Ind,fill=Gene,geom="tile",width=Length,label=Gene) + geom_text() + theme_bw() + theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),plot.margin=unit(c(.5,.5,.5,.55),"cm")) + ylab(" ") + xlab("Custom Numbering") + scale_x_continuous(breaks=exonRange,labels=paste(Index[exonRange]))
       

            }else{
 
            qplot(data=Genes,MP,Ind,fill=Gene,geom="tile",width=Length,label=Gene) + geom_text() + theme_bw() + theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),plot.margin=unit(c(.5,.5,.5,.55),"cm")) + ylab(" ") + xlab(" ")
        }
    }
    })



    output$CIplot<-renderPlot({
        
        if(input$selVar1=="None"){
            plot(NULL,xlim=c(1,10),ylim=c(0,1000))
        }else{
            Sample<-cnv.calls[input$selVar1,]$sample
            exonRange<-input$minEx1:input$maxEx1
            refs_sample<-refs[[Sample]]
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
            CIPlot<-ggplot(CIData,aes(x=Exon,y=Ratio))+geom_ribbon(aes(ymin=Min,ymax=Max),fill="grey")+geom_point(col="blue",cex=3.5) + theme_bw()+xlab("Exon Index")+ylab("Observed/Expected") + scale_x_continuous(breaks=exonRange,labels=paste(exonRange))

            temp = cnv.calls[cnv.calls$sample==Sample,]
            if(sum(temp$start.p%in%exonRange |temp$end.p%in%exonRange)>0){
                temp = temp[temp$start.p%in%exonRange|temp$end.p%in%exonRange,]
                for(i in 1:nrow(temp)){
                    start.temp = temp[i,]$start.p
                    end.temp = temp[i,]$end.p
                    CIPlot<-CIPlot + geom_point(data=CIData[CIData$Exon%in%start.temp:end.temp,], aes(x=Exon,y=Ratio),color="red",cex=3.5) 
                }
            }
#	    my_grob = grobTree(textGrob("Confidence Region", x=0.8,  y=0.9, hjust=0,gp=gpar(col="black", fontsize=12, fontface="bold")))
#	    CIPlot<-CIPlot +  annotation_custom(my_grob)
            print(CIPlot)
        }
    })

})

