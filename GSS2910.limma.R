#############################################
#choose sample 2-73
#
#############################################
rm(list=ls())
library(limma)
library(Glimma)
library(edgeR)

args <- commandArgs(trailingOnly = TRUE)
#filename <- args[1]
filename="gene.count.GSS2910"
cnts=read.csv(filename,sep='\t')
d=dim(cnts)
geneID=cnts[,1]

samplenames <- colnames(cnts)[2:d[2]]
splitname<-strsplit(samplenames, "[._]")
Clen=length(samplenames)
trt=rep("NA", Clen)
id=rep("NA", Clen)
dup=rep("NA", Clen)
study=rep("NA", Clen)
for(mm in  1:Clen ){
  trt[mm]=paste0(splitname[[mm]][2], splitname[[mm]][3])
  study[mm]=splitname[[mm]][1]
  id[mm]=splitname[[mm]][5]
  dup[mm]=splitname[[mm]][4]
}
group <- as.factor(trt)
samplenames=paste0(trt, ".", dup)

cnts = cnts[,2:d[2]]
colnames(cnts)=samplenames
rownames(cnts)=geneID
sum(duplicated(geneID))####check if there is any geneid is duplicated

x <- DGEList(counts=cnts, group=group)


#####基于序列的CPM（counts per million）、log-CPM、FPKM（fragments per kilobase of transcript per million），和基于转录本数目的RPKM（reads per kilobase of transcript per million）。
##-----remove low expression data if cpm<=1 in more than 3 samples-----
cpm <- cpm(x)
###use edge R function, average count to be added to each observation to avoid taking log of zero. Used only if log=TRUE.
lcpm <- cpm(x, log=TRUE, prior.count=2)
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)

summary(lcpm)
#####删除低表达基因
table(rowSums(x$counts==0)==Clen)


keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)


##----generate before after figure for filtering-----------------
png("1filter.png", height=3600, width=3800, res=300)
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)

col1=brewer.pal(8, "Dark2")
col2=brewer.pal(9, "Set1")
col3=brewer.pal(9, "Pastel1")
col4=brewer.pal(8, "Accent")
col5=brewer.pal(8, "Set3")

nsamples <- ncol(x)

col <-  c(col1,col2,col3,col4,col5)[1:nsamples]

par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
dev.off()

## 归一化基因表达分布----normalize------------------------------------------
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

## ----normalizemodifieddata--
###I could not see big difference and do not understand why this is compare before and after?--------------------------
x2 <- x
x2$samples$norm.factors <- 1
png("2normalize.png", height=3000, width=8000, res=300)
library(limma)
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")
dev.off()

## 对样本的无监督聚类----MDS1
library(RColorBrewer)
col1=brewer.pal(8, "Dark2")
col2=brewer.pal(9, "Set1")
col3=brewer.pal(9, "Pastel1")
col4=brewer.pal(8, "Accent")
col5=brewer.pal(8, "Set3")
png("3mds.group.png")
lcpm <- cpm(x, log=TRUE)
col.group <- group

levels(col.group) <-  c(col1,col2,col3,col4,col5)[1:nlevels(col.group)]
col.group <- as.character(col.group)
plotMDS(lcpm, labels=samplenames, col=col.group)
#plotMDS(lcpm, labels=samplenames, groups=group,legend="all")
title(main="A. Sample groups")
dev.off()

## ----GlimmaMDSplot----if LAUCH= TRUE will launch into html-----------------------------

library(Glimma)
glMDSPlot(lcpm, labels=samplenames, 
          groups=group, launch=FALSE, legend="all")

## 差异表达分析----design---------------------------------------
###创建设计矩阵和对比

design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design

## ----contrasts----------------
contr.matrix <- makeContrasts(
   Thymol15vsCtl = Thymol15 - Control0, 
   Thymol15vsDMSO = Thymol15 - DMSO15,
   Thymol30vsCtl = Thymol30 - Control0, 
   Thymol30vsDMSO = Thymol30 - DMSO30, 
   Thymol60vsCtl = Thymol60 - Control0, 
   Thymol60vsDMSO = Thymol60 - DMSO60,
   Thymol60vs30 = Thymol60 - Thymol30,
   Thymol60vs15 = Thymol60 - Thymol15,
   Thymol30vs15 = Thymol30 - Thymol15,
   DMSO60vs30 = DMSO60 - DMSO30,
   DMSO60vs15 = DMSO60 - DMSO15,
   DMSO30vs15 = DMSO30 - DMSO15,
   DMSO15vsCtl = DMSO15 - Control0,
   DMSO30vsCtl = DMSO30 - Control0,
   DMSO60vsCtl = DMSO60 - Control0,
   levels = colnames(design))
contr.matrix

## 从表达计数数据中删除异方差----voom----当操作DGEList对象时，voom从x中自动提取文库大小和归一化因子，以此将原始计数转换为log-CPM值。在voom中，对于log-CPM值额外的归一化可以通过设定normalize.method参数来进行。


png("4voom.png", height=3600, width=3800, res=300)
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean−variance trend")
dev.off()
##拟合线性模型以进行比较 --检查DE基因数量--decidetests------
summary(decideTests(efit))

## ----treat- request fold change to be more than 2------------------------
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)


###for ArgB ###5 genes
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
####if input gene name in the $genes matrix, can be viewed here
#head(tfit$genes$SYMBOL[de.common], n=20)
## ----venn----can not do more than 5 sets============
png("5venn_1_5.png", height=1800, width=3600, res=300)
par(mfrow=c(1,2))
#####
vennDiagram(dt[,1:5], circle.col=c("turquoise", "salmon", "gold", "green", "red"), main="efit-lft=1")
dt_e <-decideTests(efit)
vennDiagram(dt_e[,1:5], circle.col=c("turquoise", "salmon", "gold", "green", "red"), main="efit")
dev.off()


write.fit(tfit, dt, file="6vfit-lft1.txt")
write.fit(efit, decideTests(efit), file="7efit.txt")
## 从上到下检查单个DE基因使用topTreat函数可以列举出使用treat得到的结果中靠前的DE基因（对于eBayes的结果可以使用topTable函数）。默认情况下，topTreat将基因按照校正p值从小到大排列，并为每个基因给出相关的基因信息、log-FC、平均log-CPM、校正t值、原始及经过多重假设检验校正的p值。列出前多少个基因的数量可由用户指定，如果设为n=Inf则会包括所有的基因。基因Cldn7和Rasef在basal与LP和basal于ML的比较中都位于DE基因的前几名。


#####----toptables---------printout log fold change data-----------this need to be changed

for (i in 1: dim(contr.matrix)[2]){
    result1 <- topTreat(tfit, coef=i, n=Inf)
    tname=colnames(contr.matrix)[i]
    write.csv(result1, file = paste0("tfit.",tname,".csv"))
    result2 <- topTreat(efit, coef=i, n=Inf)
    write.csv(result2, file = paste0("efit.",tname,".csv"))
    png(paste0("tfit",tname,".fc.png"), height=1800, width=1800, res=300)
    plotMD(tfit, column=i, status=dt[,i], main=colnames(tfit)[i], xlim=c(-8,13))
    dev.off()
    png(paste0("efit.",tname,".fc.png"), height=1800, width=1800, res=300)
    plotMD(efit, column=i, status=dt[,i], main=colnames(efit)[i], xlim=c(-8,13))
    dev.off()
    glMDPlot(tfit, coef=i, status=dt, main=colnames(tfit)[i],
         side.main="GENE", counts=x$counts, groups=group,
	 folder=paste0("tfit.", tname), launch=FALSE)
    glMDPlot(efit, coef=i, status=dt, main=colnames(efit)[i],
         side.main="GENE", counts=x$counts, groups=group,
	 folder=paste0("efit.",tname), launch=FALSE)

    #library(gplots)
    #mycol <- colorpanel(1000,"blue","white","red")
    #png(paste0("tfit.",tname,".heat.png"), height=3600, width=3600, res=300)
    #heatmap.2(v$E[result1$adj.P.Val<=0.05 ,], scale="row",labCol=group,
    #col=mycol, trace="none", density.info="none", margin=c(8,6),
    #lhei=c(2,10), dendrogram="column")
    #dev.off()

   # png("efit.KvsF.heat.png", height=3600, width=3600, res=300)
   # heatmap.2(v$E[result2$adj.P.Val<=0.05 ,], scale="row",labCol=group,
   # col=mycol, trace="none", density.info="none", margin=c(8,6), lhei=c(2,10),
   # dendrogram="column")
   # dev.off()

}

####if there is mapping files , can use map to run GSEA analysis
###check example
#load(system.file("extdata", "mouse_c2_v5p1.rda", package = "RNAseq123"))
#idx <- ids2indices(Mm.c2,id=rownames(v))
#cam.BasalvsLP <- camera(v,idx,design,contrast=contr.matrix[,1])
#head(cam.BasalvsLP,5)
#cam.BasalvsML <- camera(v,idx,design,contrast=contr.matrix[,2])
#head(cam.BasalvsML,5)

sessionInfo()
