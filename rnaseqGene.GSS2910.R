rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
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

meta<- data.frame(matrix("NA",nrow=length(samplenames), ncol=4))
colnames(meta)=c("ID","trt", "dup","id")
meta$ID=samplenames
meta$trt=trt
meta$id=id
meta$dup=dup
###########################################
library("DESeq2")
y=round(cnts)
ddsMat <- DESeqDataSetFromMatrix(countData = y,colData = meta, design = ~ trt)

#nrow(ddsMat)###14408
#keep <- rowSums(counts(dds)) > 1
# at least 3 samples with a count of 10 or higher

keep <- rowSums(counts(ddsMat) >= 10) >= 3
dds <- ddsMat[keep,]

###nrow(dds)###11287
#install.packages("hexbin")
###################################################################
#DESeq2 offers two transformations for count data that stabilize the variance across the mean: the variance stabilizing transformation (VST) for negative binomial data with a dispersion-mean trend (Anders and Huber 2010), implemented in the vst function, and the regularized-logarithm transformation or rlog (Love, Huber, and Anders 2014).For genes with high counts, both the VST and the rlog will give similar result to the ordinary log2 transformation of normalized counts. For genes with lower counts, however, the values are shrunken towards a middle value. The VST or rlog-transformed data then become approximately homoskedastic (more flat trend in the meanSdPlot), and can be used directly for computing distances between samples, making PCA plots, or as input to downstream methods which perform best with homoskedastic data.The VST is much faster to compute and is less sensitive to high count outliers than the rlog. The rlog tends to work well on small datasets (n < 30), potentially outperforming the VST when there is a wide range of sequencing depth across samples (an order of magnitude difference). We therefore recommend the VST for medium-to-large datasets (n > 30). You can perform both transformations and compare the meanSdPlot or PCA plots generated, as described below.
##########################################################
## vst: variance stablizing transformation----------------------------
###In the above function calls, we specified blind = FALSE, which means that differences between cell lines and treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment. The experimental design is not used directly in the transformation, only in estimating the global amount of variability in the counts. For a fully unsupervised transformation, one can set blind = TRUE (which is the default).
##############################################
library("vsn")
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
#head(assay(vsd), 3)
#colData(vsd)
############################################################################
## ----rlog is extremely time consuming, consider to get away from it. 
##########################################################################
#rld <- rlog(dds, blind = FALSE)
#head(assay(rld), 3)
######################################################################
## ----transformplot, fig.width = 6, fig.height = 2.5------------------------
############################################################
library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)
new_cnt<-counts(dds, normalized=TRUE)
###################################
#These plot is not as good as the limma data
#I do not like the new count glimma plot, did not seperate the sample well. 
#glMDSPlot(new_cnt, labels=group, 
#          groups=trt, launch=TRUE)

library(Glimma)
glMDSPlot(log2(cnts+1), labels=group, 
          groups=trt, launch=TRUE)
##############################################
write.csv(new_cnt, file = "Deseq.normalized.count.csv")
#write.csv(vsd, file = "Deseq.vsd.csv")
#write.csv(rld, file = "Deseq.rld.csv")
write.csv(cnts, file = "Desq.cleaned.count.csv")

####################################
#df <- bind_rows(
#  as_data_frame(log2(new_cnt[, 1:2]+1)) %>%
#         mutate(transformation = "log2(x + 1)"),
#  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
#  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

#colnames(df)[1:2] <- c("x", "y")  
#jpeg('vst-rlog-log2.jpg')
#ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
#  coord_fixed() + facet_grid( . ~ transformation)  
#dev.off()
###########################################################################
## https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
#install.packages("pheatmap")
########################################################
sampleDists <- dist(t(assay(vsd)))
sampleDists
library("pheatmap")
library("RColorBrewer")

## ----distheatmap, fig.width = 6.1, fig.height = 4.5------------------------
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
my_sample_col <- data.frame(sample = group)
row.names(my_sample_col) <- colnames(vsd)
jpeg('sampleDistance.jpg', width = 3000, height = 3000,  res = 300)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
	 #annotation_row = my_sample_col,
         col = colors)
dev.off()
#######################################################
## Another option for calculating sample distances is to use the Poisson Distance (Witten 2011), implemented in the PoiClaClu package. This measure of dissimilarity between counts also takes the inherent variance structure of counts into consideration when calculating the distances between samples. The PoissonDistance function takes the original count matrix (not normalized) with samples as rows instead of columns, so we need to transpose the counts in dds.

#install.packages("PoiClaClu")
###################################################
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
## ----poisdistheatmap, fig.width = 6.1, fig.height = 4.5--------------------
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- colnames(dds)
colnames(samplePoisDistMatrix) <- NULL

jpeg('samplePoisonDistance.jpg', width = 3000, height = 3000,  res = 300)
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,	 
         col = colors)
dev.off()
###################################################################
## ----plotpca, fig.width=6, fig.height=4.5----------------------------------
###############################################################
jpeg('PCA.jpg')
plotPCA(vsd, intgroup = "trt")
dev.off()
## --------------------------------------------------------------------------
pcaData <- plotPCA(vsd, intgroup = "trt", returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))

## ----ggplotpca, fig.width=6, fig.height=4.5--------------------------------
jpeg('PCA2.jpg')
ggplot(pcaData, aes(x = PC1, y = PC2, color = trt)) +
geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()+
  ggtitle("PCA with VST data")
dev.off()

################################################
#devtools::install_github("willtownes/glmpca")
####################################################
library("glmpca")
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
jpeg('PCA3.jpg')
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = group)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
dev.off()


################################################
##MDS plot
############################################
mds <- as.data.frame(colData(vsd))  %>%
         cbind(cmdscale(sampleDistMatrix))
jpeg('MDS.jpg')
ggplot(mds, aes(x = `1`, y = `2`, color = trt)) +
geom_point(size = 3) + coord_fixed()+ ggtitle("MDS with VST data")
dev.off()
####################################################
mdsPois <- as.data.frame(colData(dds)) %>%
   cbind(cmdscale(samplePoisDistMatrix))
jpeg('MDS2.jpg')
ggplot(mdsPois, aes(x = `1`, y = `2`, color = trt)) +
geom_point(size = 3) + coord_fixed()+ ggtitle("MDS with PoissonDistances")
dev.off()


########################Analysis##################################
results <- DESeq(dds)
for (i in unique(trt)){
  for (j in unique(trt)){
    if(i != j){
    	 comparison1 <- results(results, pAdjustMethod = "BH", contrast = c("trt", i,j))
	 comparison1a <- cbind(rownames(dds),comparison1)
	 colnames(comparison1a)[1] <- "Genes"
	 XX=2**comparison1a$log2FoldChange
	 comparison1a$TrueFC <- sapply(XX,function(x){ifelse(x>=1, x, -1/x)})
	 YY=cbind(rownames(dds), comparison1a$padj, comparison1a$TrueFC)
	 write.table(YY, file = paste0(i,"_vs_", j, ".stat"), sep="\t")
	 write.table(comparison1a, file = paste0(i, "_vs_", j, ".txt"), sep="\t")
	 print (paste0(i, "_vs_", j," p<=0.05 ",  sum(comparison1a$pvalue<=0.05, na.rm=TRUE), " padj<=0.05 ", sum(comparison1a$padj<=0.05, na.rm=TRUE)))
    }
  }
}


sessionInfo()