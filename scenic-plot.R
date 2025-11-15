Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

library(Seurat)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(doParallel)
library(foreach)

load("sc.m.rdata")

mydbs <- c( "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather","hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
names(mydbs) <- c("10kb","500bp")



data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=1,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "epi")
saveRDS(scenicOptions, "int/scenicOptions.rds")


######################################################################


regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")


#####################################################################################
### H   M1   S1
colnames(regulonAUC)[1:10]


sc.hms= sc.epi 
table(sc.hms$celltype.3)


Idents(sc.hms)="celltype.3"
sc.CD8_Naive= sc.hms 

X=intersect(colnames(regulonAUC),rownames(sc.CD8_Naive@meta.data)) 
regulonAUC.hms=regulonAUC[,rownames(sc.CD8_Naive@meta.data)]

cellInfo.hms=sc.CD8_Naive@meta.data

rss <- calcRSS(AUC=getAUC(regulonAUC.hms), cellAnnotation=cellInfo.hms[colnames(regulonAUC.hms), "group1" ])

rss=na.omit(rss)
colnames(rss)
rss=rss[, c("NC","Tumor1","Tumor2")]
#可视化细胞特异性TF
rssPlot <- 
  plotRSS(
    rss ,
    zThreshold = 1.7,
    cluster_columns = FALSE,
    order_rows = TRUE,
    thr=0.1,
    varName = "group1",
    col.low = '#330066',
    col.mid = '#66CC66',
    col.high = '#FFCC33')


rssPlot

dev.off()
B_rss <- as.data.frame(rss)#rss特异性TF结果

colnames(B_rss)
celltype <- colnames(B_rss)
rssRanklist <- list()
library(ggrepel)
for(i in 1:length(celltype)) {
  
  data_rank_plot <- cbind(as.data.frame(rownames(B_rss)),
                          as.data.frame(B_rss[,celltype[i]]))#提取数据
  
  colnames(data_rank_plot) <- c("TF", "celltype")
  data_rank_plot=na.omit(data_rank_plot)#去除NA
  data_rank_plot <- data_rank_plot[order(data_rank_plot$celltype,decreasing=T),]#降序排列
  data_rank_plot$rank <- seq(1, nrow(data_rank_plot))#添加排序
  
  p <- ggplot(data_rank_plot, aes(x=rank, y=celltype)) + 
    geom_point(size=1, shape=16, color="#1F77B4",alpha =0.4)+
    geom_point(data = data_rank_plot[1:6,],
               size=1, color='#DC050C')+ #选择前6个标记，自行按照需求选择
    theme_bw()+
    theme(axis.title = element_text(colour = 'black', size = 12),
          axis.text = element_text(colour = 'black', size = 10),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(x='Regulons Rank', y='Specificity Score',title =celltype[i])+
    geom_text_repel(data= data_rank_plot[1:6,],
                    aes(label=TF), color="black", size=3, fontface="italic", 
                    arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.2,
                    point.padding = 0.3, segment.color = 'black', 
                    segment.size = 0.3, force = 1, max.iter = 3e3)
  rssRanklist[[i]] <- p
}


library(cowplot)

plot_grid(rssRanklist[[1]],rssRanklist[[2]],rssRanklist[[3]]  ,ncol=3)

tf.gene =  c(rssRanklist[[1]][["data"]]$TF[1:40],rssRanklist[[2]][["data"]]$TF[1:40] )

write.csv(c(rssRanklist[[1]][["data"]]$TF[1:40],rssRanklist[[2]][["data"]]$TF[1:40]  
            ,rssRanklist[[3]][["data"]]$TF[1:40] ),file = "NC.tumor.csv" )





####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################



#===============================================================================

cellInfo=sc.epi@meta.data
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$celltype.3),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))



ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[1:80,], name="Regulon activity")

################
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$group1),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[1:120,], name="Regulon activity")


################
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$orig.ident),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[1:80,], name="Regulon activity")


####################################################################







