#文章： https://www.cell.com/cancer-cell/fulltext/S1535-6108(22)00548-7#sectitle0030
#代码：https://github.com/ruoyan-li/RCC-spatial-mapping/blob/main/code/infercnv.R

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())


library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)

load("scRNA_harmony.rdata")


########################################################################################################################
#BiocManager::install("infercnv")
library("infercnv")

library("infercnv")
Sys.setenv(JAGS_HOME= "JAGS-4.3.1")

pos=read.table("human.gene.positions")
pos1=distinct(pos,V7,.keep_all = TRUE)
rownames(pos1)=pos1$V7
pos2=select(pos1,V7,V2,V3,V4)

write.table(pos2, 'geneLocate.txt', row.names=F, col.names=F, sep='\t')

scRNA1=scRNA_harmony[,sample(colnames(scRNA_harmony),500)]
exprMatrix <- as.matrix(GetAssayData(scRNA1, slot='counts'))
cellAnnota <- subset(scRNA1@meta.data, select='seurat_clusters')
groupFiles='groupFiles.txt'
dim(exprMatrix)
write.table(cellAnnota,file =" groupFiles.txt",sep = '\t',col.names = F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=exprMatrix,
                                    annotations_file=" groupFiles.txt",
                                    delim="\t",
                                    gene_order_file= "geneLocate.txt",
                                    ref_group_names=NULL)


infercnv_obj = infercnv::run(infercnv_obj,write_expr_matrix=T,
                             cutoff=0.1, # use 1 for smart-seq, 0.1 for 10x-genomics 
                             out_dir=  'cnv1/' ,  # dir is auto-created for storing outputs
                             cluster_by_groups=F,   #  cluster_by_groups：先区分细胞来源，再做层次聚类
                             hclust_method="ward.D2",denoise=T,HMM=T, plot_steps=F)


######################################################################################3
table(scRNA1$celltype)
cellAnnota <- subset(scRNA1@meta.data, select='celltype')
groupFiles='groupFiles1.txt'
dim(exprMatrix)
write.table(cellAnnota,file =" groupFiles1.txt",sep = '\t',col.names = F)

table(cellAnnota$celltype)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=exprMatrix,
                                    annotations_file=" groupFiles1.txt",
                                    delim="\t",
                                    gene_order_file= "geneLocate.txt",
                                    ref_group_names="Macrophage")
?run
infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1, # use 1 for smart-seq, 0.1 for 10x-genomics 
                              out_dir=  'cnv.REF/' ,  # dir is auto-created for storing outputs
                              cluster_by_groups=T,   #  cluster_by_groups：先区分细胞来源，再做层次聚类
                              hclust_method="ward.D2", plot_steps=F,denoise=TRUE,
                              HMM=F,  ##特别耗时间
                              num_threads=30,write_expr_matrix=T
)


#################################################################
cnv_score_table = data.table::fread("cnv.REF/infercnv.observations.txt", 
                                    data.table = F) %>% 
  column_to_rownames(var = 'V1')



library(scales)
cnvScore <- function(data){
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    rescale(to=c(-1, 1)) %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}

cnv_score <- cnvScore(cnv_score_table)
library(ggpubr)

cnv_score <- cbind(cnv_score, cellAnnota[row.names(cnv_score),])
names(cnv_score) <- c('cnv_score', 'celltype')
# 绘制箱线图
color <- ggsci::pal_aaas()(10)
ggplot(cnv_score, aes(reorder(celltype, cnv_score), cnv_score, color = celltype)) +
  geom_boxplot() +
  # scale_color_manual(values = color) +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "NA") +
  labs(x = '', y = 'CNV Scores', title = '') +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1)) +
  stat_compare_means()

library(ggplot2)
library(ggridges)

ggplot(cnv_score, aes(x = cnv_score, y = celltype)) + geom_density_ridges()

ggplot(cnv_score, aes(x = cnv_score, y = celltype, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)+  
  theme_bw(base_size=10)+
  theme(axis.text = element_text(colour = 'black'))

########################################################################################################

library(RColorBrewer)
infercnv::plot_cnv(infercnv_obj2, #上两步得到的infercnv对象
                   plot_chr_scale = T, #画染色体全长，默认只画出（分析用到的）基因
                   output_filename = "better_plot",output_format = "pdf", #保存为pdf文件
                   custom_color_pal =  color.palette(c("#8DD3C7","white","#BC80BD"), c(2, 2))) #改颜色




