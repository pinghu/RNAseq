remove(list=ls())
library(ComplexHeatmap)
filename="example_gene_data"
outname="example"
args <- commandArgs(trailingOnly = TRUE)
#filename <- args[1]
#outname <- args[2]
#pathID <- args[1]
#filename = paste0(pathID, ".data")
#outname = pathID
A<-read.table(filename,header=TRUE, sep="\t",stringsAsFactors=F)
d<-dim(A)
rownames(A) <- paste0( A[,107], ".",A[,1])
A <- A[order(rownames(A)), ]
lfc <- A[, 2:17]

#rownames(lfc) <- rownames(A) 
## Remove ".lfc" from column names
colnames(lfc) <- sub("\\.log2FoldChange$", "", colnames(lfc))
p <- A[,19:34]
adjp<-A[,36:51]
lfcse <-A[,53:68]
basemean <-A[,70:85]
#pathname <-A[2,95]
#rownames(p)<- rownames(lfc)
#rownames(adjp) <- rownames(lfc)
colnames(p) <-sub("\\.pvalue$", "", colnames(p))
colnames(adjp) <-sub("\\.padj$", "", colnames(adjp))
colnames(lfcse) <-sub("\\.lfcSE$", "", colnames(lfcse))
colnames(basemean) <-sub("\\.baseMean$", "", colnames(basemean))


# Create an empty matrix with the same dimensions as p
p_mat <- matrix("", nrow = nrow(p), ncol = ncol(p))

# # Set the row and column names
rownames(p_mat) <- rownames(p)
colnames(p_mat) <- colnames(p)
p_mat[!is.na(p) & p <= 0.05] <- "*"
p_mat[!is.na(adjp) & adjp <= 0.05] <- "**"
# 
# # Create row annotation with row names
# data <- as.matrix(lfc)
# d=dim(data)
# row_avg <- rowMeans(data)
# col_avg <- colMeans(data)
# 
# library(circlize)
# 
# color_scale = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
# 
# row_ha <- rowAnnotation(LFC = anno_barplot(row_avg))
# column_ha <- HeatmapAnnotation( LFC= anno_barplot(col_avg))
# bacteria <-as.factor(A[, 110])
# ht_list<-Heatmap(data, name = "Log2FoldChange", top_annotation = column_ha, 
#                  left_annotation = row_ha,  row_names_gp = gpar(fontsize = 7),  
#                  column_names_gp = gpar(fontsize = 9),
#                  cluster_rows = FALSE, cluster_columns  = FALSE,
#                  #rect_gp = gpar(col = "white", lwd = 2),
#                  row_split =bacteria,#column_title = paste(pathID,pathname),
#                  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
#                    grid.text(p_mat[i, j], x, y)
#                  },
#                  col = color_scale)
# 
# myheight=3
# if (d[1] >20){
#   myheight=0.16*d[1]
# }
# jpeg(paste0(outname,".jpg"), width = 6, height = myheight, units = "in", res = 300)
# draw(ht_list)#, heatmap_legend_side = "top")
# dev.off()




#http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
library(ggplot2)
# df <- ToothGrowth
# df$dose <- as.factor(df$dose)
# head(df)
# 
# #+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(ToothGrowth, varname="len", 
                    groupnames=c("supp", "dose"))
# Convert dose to a factor variable
df2$dose=as.factor(df2$dose)
head(df2)

#p_color = rep("blue", dim(lfc)[2])
p_color <- matrix("gray", nrow = nrow(p), ncol = ncol(p))
p_color[!is.na(p) & p <= 0.05] <- "orange"
p_color[!is.na(adjp) & adjp <= 0.05] <- "red"
i=1
df3 = cbind(colnames(lfc),as.numeric(lfc[i,]), as.numeric( lfcse[i,]),
          as.numeric( p[i,]),as.numeric( adjp[i,]), as.numeric(basemean[i,]), p_color[i,], p_mat[i,])
colnames(df3)=c("treatment", "log2FoldChange", "lfcSE", "pValue", "adjP", "basemean", "p_color", "p_sig")
df3 <- as.data.frame(df3)
df3$log2FoldChange <- as.numeric(df3$log2FoldChange)
df3$lfcSE <- as.numeric(df3$lfcSE)
library(ggplot2)
# Assuming df3 is already prepared and includes p_sig column
# Plotting code

p <-ggplot(df3, aes(x=treatment, y=log2FoldChange, fill=p_color)) +
  geom_bar(stat="identity", 
           #aes(fill=ifelse(log2FoldChange > 0, "red", "blue")), 
           color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=log2FoldChange - lfcSE, ymax=log2FoldChange + lfcSE), 
                width=.2, position=position_dodge(.9)) +
  geom_text(aes(label=p_sig, y=ifelse(log2FoldChange > 0, log2FoldChange + lfcSE+0.1, log2FoldChange - lfcSE-0.1)),
            position=position_dodge(.9))
# Print the plot
print(p)


# Default bar plot
p<- ggplot(df3, aes(x=treatment, y=log2FoldChange, fill=p_color)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width=.2,
                position=position_dodge(.9))
print(p)
# Finished bar plot
p+labs(title="Tooth length per dose", x="Dose (mg)", y = "Length")+
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00'))

# Keep only upper error bars
ggplot(df2, aes(x=dose, y=len, fill=supp)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=len, ymax=len+sd), width=.2,
                position=position_dodge(.9))

# Default line plot
p<- ggplot(df2, aes(x=dose, y=len, group=supp, color=supp)) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
                position=position_dodge(0.05))
# Note that, you can chose to keep only the upper error bars
print(p)
# Finished line plot
p+labs(title="Tooth length per dose", x="Dose (mg)", y = "Length")+
  theme_classic() +
  scale_color_manual(values=c('#999999','#E69F00'))

# Use geom_pointrange
ggplot(df2, aes(x=dose, y=len, group=supp, color=supp)) +
  geom_pointrange(aes(ymin=len-sd, ymax=len+sd))
# Use geom_line()+geom_pointrange()
ggplot(df2, aes(x=dose, y=len, group=supp, color=supp)) +
  geom_line()+
  geom_pointrange(aes(ymin=len-sd, ymax=len+sd))

p <- ggplot(df, aes(x=dose, y=len)) +
  geom_dotplot(binaxis='y', stackdir='center')
# use geom_crossbar()
p + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                 geom="crossbar", width=0.5)
# Use geom_errorbar()
p + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color="red", width=0.2) +
  stat_summary(fun.y=mean, geom="point", color="red")

# Use geom_pointrange()
p + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="pointrange", color="red")