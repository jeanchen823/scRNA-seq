Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(cowplot)
library(ggplot2)
load("scRNA_harmony.rdata")

##BiocManager::install("AUCell")
library(AUCell)
library(ggplot2)
library(Seurat)
##mn BiocManager::install("clusterProfiler")
library(clusterProfiler)

sc.id=sample(colnames(scRNA_harmony),1500)
sc2=scRNA_harmony[,sc.id]
##install.packages("doParallel")
##install.packages("doRNG")
exp=GetAssayData(sc2,slot = "data")
cells_rankings <- AUCell_buildRankings(exp,  nCores=6, plotStats=TRUE) 

cells_rankings

c2 <- read.gmt("c2.cp.kegg_medicus.v2023.2.Hs.symbols .gmt") 
geneSets <- lapply(unique(c2$term), function(x){print(x);c2$gene[c2$term == x]})
names(geneSets) <- unique(c2$term)
?AUCell_calcAUC
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,nCores =1,
                            aucMaxRank=nrow(cells_rankings)*0.1)



grep("OX",rownames(cells_AUC@assays@data$AUC),value = T)

geneSet <- "KEGG_MEDICUS_REFERENCE_BETA_OXIDATION" 
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
sc2$AUC <- aucs
df<- data.frame(sc2@meta.data, sc2@reductions$umap@cell.embeddings)
colnames(df)
class_avg <- df %>%
  group_by(celltype) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )

ggplot(df, aes(umap_1, umap_2))  +
  geom_point(aes(colour  = AUC)) + viridis::scale_color_viridis(option="D") +
  ggrepel::geom_label_repel(aes(label = celltype),
                            data = class_avg,
                            size = 3,
                            label.size = 1,
                            segment.color = NA
  )+   theme(legend.position = "none") + theme_bw() + facet_grid(.~orig.ident)

colnames(df)

library(viridis)
ggplot(data.frame(sc2@meta.data,sc2@reductions$umap@cell.embeddings), aes(umap_1, umap_2, color=AUC)
) + geom_point( size=1.5
) + scale_color_viridis(option="A")  + theme_light(base_size = 26) + facet_grid(.~orig.ident)

ggplot(data.frame(cd8.seurat@meta.data, cd8.seurat@dr$umap@cell.embeddings), aes(UMAP1, UMAP2, color=AUC)
) + geom_point( size=1.5
) + scale_color_viridis(option="A")  + theme_light(base_size = 26) + facet_grid(.~Type)


##install.packages("msigdbr")
library(msigdbr)
msigdbr_species()
msigdbr_collections()

m_df<- msigdbr(species = "Homo sapiens",  category = "C2", subcategory = "KEGG")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)



#############################################################################################
###  gsva
library('GSEABase')
library(GSVA)
library(msigdbr)
m_df<- msigdbr(species ="Homo sapiens",  category = "H" )
geneSets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
colnames(scRNA_harmony@meta.data)
table(scRNA_harmony$celltype)
Idents(scRNA_harmony)="celltype"
exp=AverageExpression(scRNA_harmony) 
exp=exp[["RNA"]]
exp=as.matrix(exp)

GSVA_hall <- gsva(expr=as.matrix(exp), 
                  gset.idx.list=geneSets, 
                  mx.diff=T, # 数据为正态分布则T，双峰则F
                  kcdf="Gaussian", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"，
                  parallel.sz=4) # 并行线程数目
head(GSVA_hall)
pheatmap::pheatmap(GSVA_hall[1:20,], #热图的数据
                   cluster_rows = T,#行聚类
                   cluster_cols =T,#列聚类，可以看出样本之间的区分度
                   
                   show_colnames=T,
                   scale = "row", #以行来标准化，这个功能很不错
                   color =colorRampPalette(c("blue", "white", "red"))(100))




exp=AverageExpression(scRNA_harmony,add.ident = "orig.ident") 
exp=exp[["RNA"]]

counts2=exp[,c(1:6)]

GSVA_hall <- gsva(expr=as.matrix(counts2), 
                  gset.idx.list=geneSets, 
                  mx.diff=T, # 数据为正态分布则T，双峰则F
                  kcdf="Gaussian", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"，
                  parallel.sz=4) # 并行线程数目
head(GSVA_hall)

pheatmap::pheatmap(GSVA_hall, #热图的数据
                   cluster_rows = T,#行聚类
                   cluster_cols =F,#列聚类，可以看出样本之间的区分度
                   
                   show_colnames=T,
                   scale = "row", #以行来标准化，这个功能很不错
                   color =colorRampPalette(c("#FF7744", "white","#AAAAAA","#0044BB"))(100))

#################################################################################################

##install.packages("msigdbr")

##BiocManager::install("GSVA")
library('GSEABase')
library(GSVA)
library(msigdbr)
m_df<- msigdbr(species ="Homo sapiens",  category = "H" )
geneSets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

table(scRNA_harmony$celltype,scRNA_harmony$orig.ident)

Idents(scRNA_harmony)="celltype"
sc.b=subset(scRNA_harmony,ident="Endothelial_cells")
table(sc.b$orig.ident)
Idents(sc.b)="orig.ident"
sc.b=subset(sc.b,ident=c("sample2","sample21"))
table(sc.b$orig.ident)
exp=GetAssayData(sc.b,slot = "data")
exp=as.matrix(exp)


GSVA_hall <- gsva(expr=exp, 
                  gset.idx.list=geneSets, 
                  mx.diff=T, # 数据为正态分布则T，双峰则F
                  kcdf="Gaussian", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"，
                  parallel.sz=4) # 并行线程数目
head(GSVA_hall)

## limma差异通路分析

#BiocManager::install('limma')
library(limma)

group <- factor(sc.b@meta.data$orig.ident, levels = c( 'sample2','sample21'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(GSVA_hall)
design
# Tunor VS Normal
compare <- makeContrasts(sample2 - sample21, levels=design)
fit <- lmFit(GSVA_hall, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
head(Diff)


dat_plot <- data.frame(id = row.names(Diff),
                       t = Diff$t)

library(stringr)
dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","")

dat_plot$threshold = factor(ifelse(dat_plot$t  >-1, ifelse(dat_plot$t >= 1 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))

dat_plot <- dat_plot %>% arrange(t)

dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)

library(ggplot2)

##install.packages("ggthemes")
library(ggthemes)
# install.packages("ggprism")
library(ggprism)
p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#36638a','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
  geom_hline(yintercept = c(-1,1),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('t value of GSVA score, S2 VS S1') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p

low1 <- dat_plot %>% filter(t < -21) %>% nrow()
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
high0 <- dat_plot %>% filter(t < 1) %>% nrow()
high1 <- nrow(dat_plot)

p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                   hjust = 0,color = 'black') + # 小于-1的为黑色标签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 大于1的为黑色标签

p

ggsave("gsva_bar.pdf",p,width = 8,height  = 8)

#######################################################################################
table(scRNA_harmony$orig.ident)
table(scRNA_harmony$celltype,scRNA_harmony$orig.ident)
degdf <- FindMarkers(scRNA_harmony,ident.1 = "sample2",ident.2 = "sample21", 
                     logfc.threshold = 0.5,group.by = "orig.ident",
                     ident=1,subset.ident = "Chondrocytes",min.pct = 0.4)

degdf1=filter(degdf,p_val_adj<0.01)
##BiocManager::install("org.Mm.eg.db")

library(org.Hs.eg.db)
library(clusterProfiler)
degs.list=rownames(degdf)
erich.go.BP = enrichGO(gene =degs.list,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
dotplot(erich.go.BP,showCategory = 10 )
barplot(erich.go.BP,showCategory = 8)
erich.go.BP=erich.go.BP@result
write.table(erich.go.BP,"7.29.erich.go.BP.deg.con.hf.txt",sep = "\t",col.names = NA)


######################################################################
library(msigdbr)
msigdbr_species()
msigdbr_collections()
##pathway=read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt")
library(msigdbr)
m_df.h <- msigdbr(species ="Homo sapiens",  category = "H" )
m_df.h = m_df.h[,c(3,4)] 

?enricher
res=enricher( degs.list,TERM2GENE =m_df.h )

dotplot(res)


### kegg
m_df.kegg<- msigdbr(species ="Homo sapiens",  category = "C2", subcategory = "KEGG" )
m_df.kegg = m_df.kegg[,c(3,4)] 

res.kegg=enricher( degs.list,TERM2GENE =m_df.kegg )

dotplot( res.kegg)


############################################################################################################################################

colnames(erich.go.BP)
###install.packages("GOplot")
library(GOplot)


go1=erich.go.BP[1:3,c(1,2,8,6)]
rownames(go1)=go1$ID
go=go1
###install.packages("stringr")
library(stringr)
go$geneID=str_replace_all(go$geneID,"/",",")
names(go)=c('ID','Term','Genes','adj_pval')
go$Category="BP"



x1=strsplit(go$Genes[1],split=",",fixed=T)
x2=strsplit(go$Genes[2],split=",",fixed=T)
x3=strsplit(go$Genes[3],split=",",fixed=T)
g1=c(x1[[1]],x2[[1]],x3[[1]])


genedata1=degdf[g1,]   
genedata1$ID=rownames(genedata1)
genedata2=data.frame(ID=genedata1$ID,logFC=genedata1$avg_log2FC)
#genedata=data.frame(ID=degs.list,logFC=degdf[degs.list,]$avg_log2FC)


circ <- circle_dat(go,genedata2)
##条形图
GOBar(subset(circ, category == 'BP'))
#气泡图
GOBubble(circ, labels = 3)

GOCircle(circ, nsub = 3)


chord<-chord_dat(circ, genedata2)

GOChord(chord, gene.order = 'logFC')

##基因与GO Term的热图(GOHeat)
GOHeat(chord, nlfc = 1, fill.col = c('red', 'yellow', 'green'))




#################################################################3
############################################################################################################################################
library(org.Mm.eg.db)
library(clusterProfiler)
degs.list=rownames(degdf)
erich.go.BP = enrichGO(gene =degs.list,
                       OrgDb = org.Mm.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
dotplot(erich.go.BP,showCategory = 10 )
barplot(erich.go.BP,showCategory = 8)
erich.go.BP=erich.go.BP@result
write.table(erich.go.BP,"erich.go.BP.deg.con.hf.txt",sep = "\t",col.names = NA)

kegg=read.table("go.txt",header = T,sep = "\t")

k = data.frame(kegg)
library(ggplot2)
library(dplyr)
before <- as.numeric(sub("/\\d+$", "", k$GeneRatio))
after <- as.numeric(sub("^\\d+/", "", k$GeneRatio))
k$GeneRatio = before /after
font.size =10

k %>% 
  ## 对进行p值排序
  arrange(p.adjust) %>% 
  ##指定富集的通路数目
  dplyr::slice(1:14) %>% 
  ## 开始ggplot2 作图，其中fct_reorder调整因子level的顺序
  ggplot(aes(GeneRatio,forcats::fct_reorder(Description,Count)))+ 
  ## 画出点图
  geom_point(aes(color=p.adjust, size = Count)) +
  ## 调整颜色，guide_colorbar调整色图的方向
  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  ## 调整泡泡的大小
  scale_size_continuous(range=c(3, 8))+
  ## 如果用ylab("")或出现左侧空白
  labs(y=NULL) +
  ## 如果没有这一句，上方会到顶
  ggtitle("")+
  ## 设定主题
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))




erich.go.CC = enrichGO(gene =degs.list,
                       OrgDb = org.Mm.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
dotplot(erich.go.CC,showCategory = 8)

barplot(erich.go.CC,showCategory = 8)

erich.go.MF = enrichGO(gene =degs.list,
                       OrgDb = org.Mm.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
dotplot(erich.go.MF,showCategory = 8)

barplot(erich.go.MF,showCategory = 8)

keytypes(org.Hs.eg.db)


DEG.entrez_id = mapIds(x = org.Mm.eg.db,
                       keys =  degs.list,
                       keytype = "SYMBOL",
                       column = "ENTREZID")

## install.packages("R.utils")

library(R.utils)

R.utils::setOption("clusterProfiler.download.method","auto")

erich.kegg.res <- enrichKEGG(gene = DEG.entrez_id,
                             organism = "mmu",
                             keyType = "kegg")

dotplot(erich.kegg.res)
kegg.res.df=erich.kegg.res@result
write.table(kegg.res.df,"erich.kegg.res.txt",sep = "\t",col.names = NA)




