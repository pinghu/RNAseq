rm(list=ls())
#args <- commandArgs(trailingOnly = TRUE)
#filename <- args[1]
filename="GSS2910_11_deseq.sig_stat.txt"

test1=read.csv(filename,sep='\t')
NormalizedCnt=test1[, 46:89]
d=dim(NormalizedCnt)
geneID=test1[,1]
geneName=test1[,159]
geneSyn=test1[,160]
samplenames <- colnames(NormalizedCnt)
splitname<-strsplit(samplenames, "[.]")
Clen=length(samplenames)
trt=rep("NA", Clen)
study=rep("NA", Clen)
dup=rep("NA", Clen)
for(mm in  1:Clen ){
  trt[mm]=splitname[[mm]][2]
  study[mm]=splitname[[mm]][4]
  dup[mm]=splitname[[mm]][3]
}
group <- as.factor(paste0(study,".",trt))
samplenames=paste0(trt, ".", dup)
colnames(NormalizedCnt)=samplenames
rownames(NormalizedCnt)=geneID
meta<- data.frame(matrix("NA",nrow=length(samplenames), ncol=4))
colnames(meta)=c("ID","trt", "dup","study")
meta$ID=samplenames
meta$trt=group
meta$study=study
meta$dup=dup
###########################################
library("DESeq2")
y=round(NormalizedCnt)
#ddsMat <- DESeqDataSetFromMatrix(countData = y,colData = meta, design = ~ trt)
#keep <- rowSums(counts(ddsMat) >= 10) >= 3
#dds <- ddsMat[keep,]
#dds <- DESeq(dds, quiet=TRUE)
#glimmaMA(dds)
#library(Glimma)

#glMDSPlot(log2(NormalizedCnt+1), labels=group, 
#          groups=group, launch=TRUE)
#glimmaMDS(log2(NormalizedCnt+1), labels=group, 
 #         groups=group, launch=TRUE)

library(ggplot2)
library(tidyverse)

d<-dim(NormalizedCnt)
for (i in 1:d[1]){
  i=9
  ann=paste(geneID[i], geneName[i], geneSyn[i])
  gene0<-as.numeric(NormalizedCnt[i,])
  gene<-log(gene0)
  mydata0=data.frame(gene0,gene, group, trt,study, dup)
  #rownames(NormalizedCnt)[i]
  png(paste0(geneID,".png"))
  ggplot(mydata0, aes(x=group, y=gene, fill=group)) +
    geom_boxplot() +geom_point( aes(fill = group),shape = 21, size = 2) +
    coord_flip() +
    ggtitle(paste(ann, "log(count)")) +
    xlab("") + ylab("")+theme_bw()+theme(legend.position = "none")
  dev.off()
}

sessionInfo()