
#-------------------------------------------------------------------#
# Author:  Dionne Swift
# Date:    September 22, 2022
# Purpose: Plots results 
#-------------------------------------------------------------------#
#
#------------------------------------------------------------------#

library(dplyr); library(ggplot2); library(tidyverse);library(ggbeeswarm); library(gridExtra)
rm(list=ls())

#---- Reading Data ----#
path <- "C:/Users/swift.dp/OneDrive - Procter and Gamble/Documents/DMS/Microbiome/GSS Projects/GSS2911_2910/"
filename <- "GSS2910_11_deseq.sig_stat.txt"
cnt_wide <- read.csv(paste0(path,filename),sep='\t')
geneID   <-cnt_wide[,1]
geneName <-cnt_wide[,159]
geneSyn  <-cnt_wide[,160]
ngene <- dim(cnt_wide)[1]

cnt_wide        <- cnt_wide[,c(1,46:89)]
names(cnt_wide) <-gsub("NormalizedCount.","",names(cnt_wide))

#---- Creating Long Data ----#
cnt_lng <- cnt_wide %>% 
  pivot_longer( cols      = !gene,
                names_to     = "trt_rep_study",
                values_to = "norm_cnt")  %>%
     separate(trt_rep_study, c("trt","rep","study"),"[.]")                              

cnt_lng$time <- as.numeric(gsub(".*?([0-9]+).*", "\\1", cnt_lng$trt))
cnt_lng$ntrt <- gsub("*?([0-9]+)", "", cnt_lng$trt)

#---- Check---#
# table(cnt_lng$ntrt)

#----- Datasets for plotting----#
cnt_lng_2910 <- subset(cnt_lng, study=="GSS2910")
cnt_lng_2911 <- subset(cnt_lng[,-6], study=="GSS2911")
cnt_lng_all  <- cnt_lng  %>%
                 mutate(ntime=ifelse(trt == "DMSOD", 0,
                                ifelse(trt ==  "ThymolI" ,15,
                                  ifelse(trt == "ThymolC", 30,
                                    ifelse(trt == "ThymolS", 60,time)))),
                        cond=ifelse(study=="GSS2910","Solution","Surface"))
                                 

#------------------------------#
#   All Data 
#------------------------------#

#------------------- lineplots--------#

path2 <- "Z:/NextSeq500/GSS2911_RS_E coli_Thymol_Surface/bioinformatics/Plots/"



for (i in 1:ngene){
 # i = 9
  ann            <- paste(geneID[i], geneName[i], geneSyn[i])
  mydata0        <- cnt_lng_all[cnt_lng_all$gene==geneID[i],]
  mydata0$logcnt <-log(mydata0$norm_cnt)
  
  df_mean <-  mydata0 %>% 
    dplyr::group_by(ntrt, ntime) %>% 
    dplyr::summarize(average = mean(logcnt)) %>%
    ungroup()
  cntl_mean = df_mean$average[df_mean$ntrt=="Control"]

  minval = min(mydata0$logcnt); maxval=max(mydata0$logcnt)
  
  solu <- mydata0[mydata0$cond=="Solution",]
  solu[nrow(solu) + 1,] <-data.frame(gene=geneID[i],trt="DMSO",  rep="1",study="GSS2910",norm_cnt=0,time=0,ntrt="DMSO",  ntime=0,cond="Solution",logcnt=cntl_mean)
  solu[nrow(solu) + 1,] <-data.frame(gene=geneID[i],trt="Thymol",rep="1",study="GSS2910",norm_cnt=0,time=0,ntrt="Thymol",ntime=0,cond="Solution",logcnt=cntl_mean)
  solu$ntrt   <- factor(solu$ntrt)
  solu$ntime  <- factor(solu$ntime)
  
  surf <- mydata0[mydata0$cond=="Surface",]
  surf$ntrt   <- factor(surf$ntrt)
  surf$ntime  <- factor(surf$ntime)
  
   
  p1= ggplot(solu, aes(x=ntime, y=logcnt)) +
        geom_point(aes(col=ntrt,shape=ntrt))  +   ylim(minval, maxval) +
          stat_summary(aes(x=ntime, y=logcnt,col=ntrt), geom = "point", fun = mean) + 
            stat_summary(aes(x=ntime, y=logcnt,col=ntrt,group=ntrt), geom = "line",  fun = mean) + 
              labs(x="Time (mins)", y="Log(Normalized Count)",col="Treatment",shape="Treatment")+theme_bw() +
                theme(legend.position="bottom")
  
  
  p2= ggplot(surf, aes(x=ntime, y=logcnt)) +
        geom_point(aes(shape=ntrt,col=ntrt)) +  ylim(minval, maxval) +
          stat_summary(aes(x=ntime, y=logcnt,group=1), geom = "line",  fun = mean) + 
            scale_x_discrete(labels= c("DMSOD","ThymolI","ThymolC","ThymolS"))  +
              labs(x="", y="Log(Normalized Cnt)",col="Treatment",shape="Treatment")+theme_bw() +
                theme(legend.position="bottom",axis.title.y=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y=element_blank())
  
  png(paste0(path2,geneID[i],".png"))
  grid.arrange(p1, p2,ncol=2,top=paste(ann, "- log(count)"))
  dev.off()
  
}
sessionInfo()

