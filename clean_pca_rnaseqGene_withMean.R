rm(list=ls())
library("DESeq2")
library("dplyr")
library(geosphere)
library(ggplot2)
library(ggpubr)
library("RColorBrewer")
library(rstatix)
library("vsn")
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
outname <-args[2]
#filename="Porphyromonas_gingivalis_ATCC33277.ann.xls"
#outname="test"
filename="count.Pg.avg10"
outname="testpg"
drawSortedBarPlot <- function(testData, filename) {
  # Ensure testData is a data frame
  if (!is.data.frame(testData)) {
    testData <- as.data.frame(testData)
  }
  
  # Make sure testData has the expected columns
  if(!("Sample" %in% names(testData)) || !("similarityScore" %in% names(testData))) {
    stop("testData must have 'Sample' and 'similarityScore' columns")
  }
  
  # Convert similarityScore to numeric if it's not already
  testData$similarityScore <- as.numeric(as.character(testData$similarityScore))
  
  # In case of any NA introduced by conversion (e.g., non-numeric values present), handle or warn
  if(any(is.na(testData$similarityScore))) {
    warning("NAs introduced by coercion to numeric in similarityScore")
  }
  
  color_vector<- ifelse(grepl("Ctl",testData$Sample), "blue",
                        ifelse(grepl("SnF|SnCl",testData$Sample), "red", "steelblue"))
  testData$Col =color_vector
  
  sorted_df <- testData[order(testData$similarityScore), ]
  p <- ggplot(sorted_df, aes(y = reorder(Sample, similarityScore), x = similarityScore)) +
    geom_bar(stat = "identity", fill = sorted_df$Col) +
    ylab("Rank") +
    xlab("Distance Score") +
    ggtitle(paste(filename)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 12))
  
  # Save the plot as a JPEG image file
  filename_with_extension <- paste0(filename, ".rank.jpg")
  jpeg(filename_with_extension, width = 285, height = 380)
  print(p)
  dev.off()
  
  # Use rank with ties.method = "min" to handle scores that are the same
  sorted_df$order <- rank(sorted_df$similarityScore, ties.method = "min")
  
  # Specify your desired output file name
  output_tab_delimited_file <- paste0(filename, ".rank")
  
  # Write the specific columns to a tab-delimited text file
  write.table(sorted_df[, c("Sample", "order", "similarityScore")], 
              file = output_tab_delimited_file, 
              sep = "\t", 
              row.names = FALSE, 
              col.names = TRUE)
}

calculateSimilarityScores <- function(distance_dict) {
  # Flatten the distance dictionary to a numeric vector
  myDist <- unlist(distance_dict)
  
  # Find minimum and maximum distances
  min_dist <- min(myDist)
  max_dist <- max(myDist)
  
  # Calculate similarity scores based on distances from -1 to 1
  #similarityScore <- 2 * (myDist - min_dist) / (max_dist - min_dist) - 1
  # Calculate similarity scores based on distances from 0 to 1
  similarityScore <- (myDist - min_dist) / (max_dist - min_dist)
  return(similarityScore)
}

cnts=read.csv(filename,sep='\t')
d=dim(cnts)
geneID=cnts[,1]
samplenames <- colnames(cnts)[2:d[2]]
splitname<-strsplit(samplenames, "[.]")
Clen=length(samplenames)
trt=rep("NA", Clen)
id=rep("NA", Clen)
dup=rep("NA", Clen)
for(mm in  1:Clen ){
  trt[mm]=splitname[[mm]][1]
  id[mm]=splitname[[mm]][3]
  dup[mm]=splitname[[mm]][2]
}
group <- as.factor(trt)
samplenames=paste0(trt, ".", id)
cnts = cnts[,2:d[2]]
colnames(cnts)=samplenames
rownames(cnts)=geneID
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

keep <- rowSums(counts(ddsMat) >= 10) >= 3
dds <- ddsMat[keep,]

###################################################################
#DESeq2 offers two transformations for count data that stabilize the variance across the mean: the variance stabilizing transformation (VST) for negative binomial data with a dispersion-mean trend (Anders and Huber 2010), implemented in the vst function, and the regularized-logarithm transformation or rlog (Love, Huber, and Anders 2014).For genes with high counts, both the VST and the rlog will give similar result to the ordinary log2 transformation of normalized counts. For genes with lower counts, however, the values are shrunken towards a middle value. The VST or rlog-transformed data then become approximately homoskedastic (more flat trend in the meanSdPlot), and can be used directly for computing distances between samples, making PCA plots, or as input to downstream methods which perform best with homoskedastic data.The VST is much faster to compute and is less sensitive to high count outliers than the rlog. The rlog tends to work well on small datasets (n < 30), potentially outperforming the VST when there is a wide range of sequencing depth across samples (an order of magnitude difference). We therefore recommend the VST for medium-to-large datasets (n > 30). You can perform both transformations and compare the meanSdPlot or PCA plots generated, as described below.
##########################################################
## vst: variance stablizing transformation----------------------------
###In the above function calls, we specified blind = FALSE, which means that differences between cell lines and treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment. The experimental design is not used directly in the transformation, only in estimating the global amount of variability in the counts. For a fully unsupervised transformation, one can set blind = TRUE (which is the default).
##############################################

vsd <- vst(dds, blind = FALSE)
dds <- estimateSizeFactors(dds)
new_cnt<-counts(dds, normalized=TRUE)

###################################################################
## ----plotpca, fig.width=6, fig.height=4.5----------------------------------
###############################################################



pcaData <- plotPCA(vsd, intgroup = "trt", returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))



# Calculate the reference point
reference_point <- pcaData %>%
  filter(trt == "Ctl") %>%  # Filter for control group
  summarise(
    avg_PC1 = mean(PC1),  # Calculate average of PC1
    avg_PC2 = mean(PC2)   # Calculate average of PC2
  )

# Calculate the average PC1 and PC2 for each treatment group
avePCAData <- pcaData %>%
  group_by(trt) %>%
  summarise(
    avePC1 = mean(PC1, na.rm = TRUE), # Calculate the average of PC1, remove NA values
    avePC2 = mean(PC2, na.rm = TRUE)  # Calculate the average of PC2, remove NA values
  )

jpeg(paste0(outname, '.PCA2.jpg'), width=600, height=380)
ggplot(pcaData, aes(x = PC1, y = PC2, color = trt)) +
geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()+theme_bw()+
  stat_chull(aes(color=trt, fill=trt), alpha=0.1, geom="polygon") +
  geom_text(data = avePCAData, aes(label = trt, x = avePC1, y = avePC2), nudge_x =2, nudge_y = 2) + # Label first points
  ggtitle(paste0 (outname, " PCA with VST data"))
dev.off()

jpeg(paste0(outname, '.PCA0.jpg'), width=600, height=380)
ggplot(pcaData, aes(x = PC1, y = PC2, color = trt)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()+theme_bw()+
  stat_chull(aes(color=trt, fill=trt), alpha=0.1, geom="polygon") +
  ggtitle(paste(outname, "PCA with VST data"))
dev.off()

distance_dict <- list()
for (i in 1:dim(avePCAData)[1]) {
  myDist<-distm (avePCAData[i,2:3], reference_point, fun = distGeo)
  #myDist <- rdist(as.matrix(result$rotation[i,1:2, drop = FALSE]), as.matrix(reference_point))
  myName=as.character(avePCAData$trt[i])
  distance_dict[myName] <- myDist
}
pca_similarityScore=calculateSimilarityScores(distance_dict)
myResult=cbind(names(distance_dict),unlist(distance_dict), pca_similarityScore )
colnames(myResult)=c("Sample", "Distance", "similarityScore")
write.table(myResult, file = paste0(outname, ".diffRank.xls"), sep = "\t", row.names = FALSE)
drawSortedBarPlot(myResult, outname)


########################Analysis##################################
dds$trt <- relevel(dds$trt, ref = "Ctl")
dds <- DESeq(dds)
ref_level <- "Ctl"
trt_levels <- levels(dds$trt)
trt_levels <- trt_levels[trt_levels != ref_level]
length(trt_levels)
for (i in 1:length(trt_levels)) {
  group_name <- trt_levels[i]
  coef_name <- paste0("trt", i)
  if (coef_name %in% resultsNames(dds)) {
    comparison1 <- results(dds, name = coef_name, pAdjustMethod = "BH")
     #  comparison1 <- results(dds, pAdjustMethod = "BH", contrast = c("trt",i,"Ctl"))
       comparison1a <- comparison1
       colnames(comparison1a)[1] <- "Genes"
       XX=2**comparison1a$log2FoldChange
       comparison1a$TrueFC <- sapply(XX,function(x){ifelse(x>=1, x, -1/x)})
       YY=cbind(rownames(dds), comparison1a$padj, comparison1a$TrueFC)
       write.table(YY, file = paste0(outname, ".", group_name,"_vs_Ctl.padjstat"), sep="\t", row.names = FALSE)
       ZZ=cbind(rownames(dds), comparison1a$pvalue, comparison1a$TrueFC)
       write.table(ZZ, file = paste0(outname, ".", group_name,"_vs_Ctl.pstat"), sep="\t", row.names = FALSE)
       comparison1a <- cbind(rownames(dds),comparison1)
       write.table(comparison1a, file = paste0(outname, ".", group_name, "_vs_Ctl.txt"), sep="\t", row.names = FALSE)
       sigN=sum(comparison1a$pvalue<=0.05, na.rm=TRUE)
       adjSigN=sum(comparison1a$padj<=0.05, na.rm=TRUE)
       sigUpN=sum(comparison1a$pvalue<=0.05 & comparison1a$log2FoldChange>0, na.rm=TRUE)
       sigDownN=sum(comparison1a$pvalue<=0.05 & comparison1a$log2FoldChange <0, na.rm=TRUE)
       adjSigUpN = sum(comparison1a$padj<=0.05 & comparison1a$log2FoldChange >0, na.rm=TRUE)
       adjSigDownN=sum(comparison1a$padj<=0.05 & comparison1a$log2FoldChange <0, na.rm=TRUE)
       totalN=dim(cnts)[1]
       sigP=round(100 * sigN / totalN, 2)
       sigUpP=round(100 * sigUpN / totalN, 2)
       sigDownP=round(100 * sigDownN / totalN, 2)
       adjSigP=round(100 * adjSigN / totalN, 2)
       adjSigUpP=round(100 * adjSigUpN / totalN, 2)
       adjSigDownP=round(100 * adjSigDownN / totalN, 2)
       print (paste0(filename, ".",i,"|", outname, "|", group_name, "_vs_Ctl", "|", "p<=0.05:", sigN ,"|",  "padj<=0.05:", adjSigN,
                     "|", "SigUp p<=0.05 lfc>0:", sigUpN, "|", "SigDown p<=0.05 lfc<0:", sigDownN, "|",
                     "|", "adjSigUp adjp<=0.05 lfc>0:",  adjSigUpN, "|", "adjSigDown p<=0.05 lfc<0:", adjSigDownN,
                     "|Total Gene Number", totalN, "|Sig%:", sigP, "|SigUp%:", sigUpP, "|SigDown%:",sigDownP, 
                     "|adjSig%:", adjSigP, "|adjSigUp%:", adjSigUpP, "|adjSigDown%:", adjSigDownP ))
}
       
}
   
norm_counts <- counts(dds, normalized = TRUE)
write.table(norm_counts, file = paste0(filename, ".normalized_counts.tsv"), sep = "\t", quote = FALSE, col.names = NA)
# Extract sample groups (from colData)
sample_groups <- colData(dds)$trt  # Adjust "condition" to your actual column name
group_levels <- unique(sample_groups)

# Calculate group means and SEs
group_stats <- apply(norm_counts, 1, function(x) {
  tapply(x, sample_groups, function(values) {
    mean_val <- mean(values)
    se_val <- sd(values) / sqrt(length(values))
    c(mean = mean_val, se = se_val)
  })
})

# Convert group_stats (nested list) to tidy data frame
group_stats_df <- do.call(rbind, lapply(names(group_stats), function(gene) {
  group_entries <- group_stats[[gene]]
  # Extract group names and their stats
  do.call(rbind, lapply(names(group_entries), function(group) {
    vals <- group_entries[[group]]
    data.frame(
      gene = gene,
      group = group,
      mean = vals["mean"],
      se = vals["se"],
      stringsAsFactors = FALSE
    )
  }))
}))


# View a few lines
head(group_stats_df)
write.csv(group_stats_df, paste0(filename, ".group_means_and_se.csv"), row.names = FALSE)


subset_df <- group_stats_df %>%
  filter(gene == "pgi:PG_1542", group %in% c("Ctl", "SnF2_H"))


library(ggplot2)
library(dplyr)

# Rename groups and factor levels
subset_df <- subset_df %>%
  mutate(group = case_when(
    group == "SnF2_H" ~ "SnF2",
    group == "Ctl" ~ "Control",
    TRUE ~ group
  ))

# Set custom factor level order
subset_df$group <- factor(subset_df$group, levels = c("Control", "SnF2"))

# Create the bar plot with error bars
p1 <- ggplot(subset_df, aes(x = group, y = mean, fill = group)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  labs(
    title = "pgi:PG_1542 prtC Collagenase",
    x = "Group",
    y = "Mean Normalized Gene Expression"
  ) +
  scale_fill_manual(values = c("Control" = "deeppink4", "SnF2" = "blue3")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.background = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# Save the plot
png("PG1542_barplot.png", width = 800, height = 1200, res = 300)
print(p1)
dev.off()

write.csv(subset_df, "PG1542_barplot.csv", row.names = FALSE)

  




 
   
   
 

sessionInfo()

