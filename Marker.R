Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(cowplot)




load("scRNA_harmony.rdata")
###########################################################################################################
pbmc3k.final=scRNA_harmony

library(ggplot2)
library(patchwork)

features <- c("LYZ", "CCL5", "IL32", "GAPDH")

RidgePlot(pbmc3k.final, features = features, ncol = 2)

VlnPlot(pbmc3k.final, features = features, ncol = 2,pt.size=0)
VlnPlot(subset(pbmc3k.final,LYZ > 0), features ="LYZ",pt.size=0) 

FeaturePlot(pbmc3k.final, features = features,reduction = "tsne")
pbmc3k.final=RunTSNE(pbmc3k.final,reduction = "harmony",dims = 1:20)
FeaturePlot(pbmc3k.final, features = features,reduction = "tsne")

DotPlot(pbmc3k.final, features = features,scale=F)
DotPlot(pbmc3k.final, features = features,scale=F) + RotatedAxis()
DoHeatmap( pbmc3k.final, features = features, size = 3)
DoHeatmap( pbmc3k.final, features = features, size = 3,slot="data")

DoHeatmap(subset(pbmc3k.final, downsample = 100), features = features, size = 3)

FeaturePlot(pbmc3k.final, features = "MS4A1")

FeaturePlot(pbmc3k.final, features = "MS4A1", min.cutoff = 1, max.cutoff = 2)


FeaturePlot(pbmc3k.final, features = c("MS4A1", "GAPDH"), min.cutoff = "q10", max.cutoff = "q90")

FeaturePlot(pbmc3k.final, features = c("MS4A1", "CD79A"), blend = TRUE)

FeaturePlot(pbmc3k.final, features = c("MS4A1", "CD79A"), split.by = "orig.ident")


plot1 <- DimPlot(pbmc3k.final)
# Create scatter plot with the Pearson correlation value as the title
plot2 <- FeatureScatter(pbmc3k.final, feature1 = "LYZ", feature2 = "CCL5")
# Combine two plots
plot1 + plot2

baseplot <- DimPlot(pbmc3k.final, reduction = "umap")
# Add custom labels and titles
baseplot + labs(title = "huage PBMCs")+xlab("x1")
BiocManager::install("ggmin")
baseplot + ggmin::theme_powerpoint()
#########################################################################################################

dir=c("CAFs.1/","CAFs.2/","dapi1/","dapi2/")
dir=c("CAFs.1/","CAFs.2/","dapi1/","dapi2/")
names(dir) = c('CAF1',  'CAF2', 'DapiNeg1',"DapiNeg2")      



counts <- Read10X(data.dir =dir)
scRNA1 = CreateSeuratObject(counts,min.cells = 3, min.features = 200)

scRNA1[["percent.mt"]] <- PercentageFeatureSet(scRNA1, pattern = "^mt-")



VlnPlot(scRNA1,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt" ), 
        
        pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
        ncol = 3) 

scRNA2 <- subset(scRNA1, subset = nFeature_RNA > 500 & 
                   percent.mt < 20 &   nCount_RNA > 1000)



scRNA3=SplitObject(scRNA2,split.by = "orig.ident")

ifnb=scRNA2
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$orig.ident)
ifnb

##Perform analysis without integration
# run standard anlaysis workflow
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb)
ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
ifnb <- FindClusters(ifnb, resolution = 1, cluster.name = "unintegrated_clusters")
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))

##############################################################################
##############################################################################
##############################################################################

## Perform integration
ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose =T)

# re-join layers after integration
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)

ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca")
# Visualization
DimPlot(ifnb, reduction = "umap", group.by = c("orig.ident","seurat_clusters"))

###################################################################
scRNA1=ifnb

DefaultAssay(scRNA1) <- "RNA"
##BiocManager::install("SingleR")
library(SingleR)
refdata <- ref_Mouse
testdata=LayerData(scRNA1, assay="RNA", layer='data')

clusters <- scRNA1@meta.data$seurat_clusters

cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = FALSE)
scRNA1@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA1@meta.data[which(scRNA1@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}



DimPlot(scRNA1, reduction = "umap", group.by = "celltype",label = T)
DimPlot(scRNA1, reduction = "umap",label = T)

save(scRNA1,file = "scRNA1.rdata")

load("scRNA1.rdata")
#################################################################################################
##scRNA1@active.ident=scRNA1@meta.data$celltype

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)
DimPlot(scRNA1,group.by = "seurat_clusters",label=T,cols = my36colors,reduction = "umap")
DimPlot(scRNA1,  label = T,reduction = "umap")


DimPlot(scRNA1,group.by = "celltype",reduction = "tsne")
DimPlot(scRNA1,group.by = "celltype",reduction = "umap")
DimPlot(scRNA1 )
DimPlot(scRNA1 ,group.by = "celltype")


colnames(scRNA1@meta.data)
Idents(scRNA1)= "celltype"
DimPlot(scRNA1)
Idents(scRNA1)="orig.ident"
DimPlot(scRNA1)

Idents(scRNA1)= "celltype"
degs <- FindAllMarkers(scRNA1, logfc.threshold = 0.5,
                       test.use = "roc", 
                       return.thresh = 0.25, 
                       min.pct = 0.3, only.pos = T) 

degs_sig <- degs %>% 
  filter(pct.1 > 0.3 &
           power > 0.25) %>%
  filter(cluster != "other") %>%
  arrange(cluster, -power)  

# select degs for heatmap
degs_top50 <- degs_sig %>% 
  #  filter(cluster!="other") %>%
  group_by(cluster) %>% 
  top_n(50, power) %>%
  top_n(50, avg_diff) %>%
  arrange(cluster, -power)


avgData <- scRNA1@assays$RNA@data[degs_top50$gene,] %>% 
  apply(1, function(x){
    tapply(x, scRNA1$celltype, mean) # ExpMean
  }) %>% t


phData <- MinMax(scale(avgData), -2, 2) # z-score
rownames(phData) <- 1:nrow(phData)

##install.packages("pheatmap")
library(pheatmap)

phres <- pheatmap(
  phData, 
  color = colorRampPalette(c("darkblue", "white", "red3"))(99), #配色
  scale = "row",
  cluster_rows = F, #不按行聚类
  cluster_cols = T, #按列聚类
  clustering_method = "complete",
  show_rownames = F, #显示cluster名
  annotation_row = data.frame(cluster = degs_top50$cluster), 
)  

phData1= phData[, c("Fibroblasts","Adipocytes","Macrophages","Granulocytes","Monocytes","T cells")] 


phres <- pheatmap(
  phData1, 
  color = colorRampPalette(c("darkblue", "white", "red3"))(99), #配色
  scale = "row",
  cluster_rows = F, #不按行聚类
  cluster_cols = F, #按列聚类
  clustering_method = "complete",
  show_rownames = F, #显示cluster名
  annotation_row = data.frame(cluster = degs_top50$cluster), 
) 

phres
########################################################################################################
###################################################################################################


load("scRNA1.rdata")

Idents(scRNA1)="celltype"

scRNA1=scRNA2

VlnPlot(scRNA1,features = "Gapdh")
scRNA1@active.ident=factor(scRNA1@active.ident,levels = c( "T cells","Fibroblasts" ,
                                                           "Granulocytes", "Macrophages",   "Monocytes","Adipocytes"))
scRNA1$celltype1=scRNA1@active.ident

DimPlot(scRNA1, reduction = "umap", group.by = "celltype1",label = T)
DimPlot(scRNA1  )
VlnPlot(scRNA1,features = "Gapdh")


DotPlot(scRNA1,features = c("Cfd","S100a8","Gapdh"))



scRNA1@active.ident=factor(scRNA1@active.ident,levels = c("Adipocytes","Fibroblasts" , "Granulocytes", "Macrophages",   "Monocytes","T cells"))
table(scRNA1@active.ident)


marker=FindAllMarkers(scRNA1,logfc.threshold = 0.5)

write.csv(marker,file = "marker.csv")

top10 <-  marker  %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

scRNA1 <-  ScaleData(scRNA1,features = rownames(scRNA1))
DoHeatmap(scRNA1,features = top10$gene,label = F,slot = "scale.data")
table(scRNA1@active.ident)
DoHeatmap( subset(scRNA1,downsample=300) ,features = top10$gene,label = F,slot = "scale.data")


table(scRNA1$celltype)
scRNA2=scRNA1

table(scRNA1@active.ident)
scRNA1=subset(scRNA1, downsample = 150)
library(pheatmap)

colanno=scRNA1@meta.data 
colanno$barcode=rownames(colanno)
colanno=colanno%>%arrange(celltype)
rownames(colanno)=colanno$barcode
colanno$barcode=NULL



colanno1=colanno[,6]
colanno1=as.data.frame(colanno1)
rownames(colanno1)=rownames(colanno)
colnames(colanno1)="celltype"

rowanno=top10
colnames(top10)
rowanno=rowanno%>%arrange(cluster)

mat4=scRNA1[["RNA"]]@scale.data[rowanno$gene,rownames(colanno1)]
mat4[mat4>=2.5]=2.5
mat4[mat4 < (-1.5)]= -1  

pheatmap(mat4,cluster_rows = F,cluster_cols = F,
         show_colnames = F,
         annotation_col = colanno1,
         gaps_row=as.numeric(cumsum(table(rowanno$cluster))[-6]),
         gaps_col=as.numeric(cumsum(table(colanno$celltype))[-6]) 
         
)



pheatmap(mat4,cluster_rows = F,cluster_cols = F,
         show_colnames = F,
         annotation_col = colanno1,
         gaps_row=as.numeric(cumsum(table(rowanno$cluster))[-6]),
         gaps_col=as.numeric(cumsum(table(colanno$celltype))[-6]),
         color = colorRampPalette(c("darkblue", "white", "red3"))(99)
         
)



pheatmap(mat4,cluster_rows = F,cluster_cols = F,
         show_colnames = F,
         annotation_col = colanno1,
         gaps_row=as.numeric(cumsum(table(rowanno$cluster))[-6]),
         gaps_col=as.numeric(cumsum(table(colanno$celltype))[-6]) ,
         filename="marker.heatmap.pdf",width=11,height = 7
)


################################################################################################


######################################################################################################
VlnPlot(scRNA1, features = top10$gene[1:3], pt.size = 0, ncol = 1)+
  scale_x_discrete("")+
  theme(
    axis.text.x.bottom = element_blank()
  )


top10 <-  marker  %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


library(reshape2)
vln.df=as.data.frame(scRNA1[["RNA"]]@data[top10$gene,])
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("barcode","exp")

anno=scRNA1@meta.data
anno$barcode=rownames(anno)


vln.df=inner_join(vln.df,anno,by="barcode")
##vln.df$gene=factor(vln.df$gene,levels = top10$gene[1:7]) #为了控制画图的基因顺序

vln.df%>%ggplot(aes(celltype,exp))+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(vln.df$gene~.,scales = "free_y")+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_classic()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )



###########################################################################################################
## 
top10 <-  marker  %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
DotPlot(scRNA1, features = top10$gene,cols = c("lightgrey", "#FF00FF"))+
  scale_x_discrete("huage")+scale_y_discrete("")+theme_classic()+RotatedAxis()
# 
?DotPlot
bubble.df=as.matrix(scRNA1[["RNA"]]@data[ top10$gene,])
bubble.df=t(bubble.df)
bubble.df=as.data.frame(scale(bubble.df))
bubble.df$CB=rownames(bubble.df)
meta=scRNA1@meta.data[,c(1,6)]
bubble.df=merge(bubble.df,meta,by.x = "CB",by.y = 0)
bubble.df$CB=NULL

celltype_v=c()
gene_v=c()
mean_v=c()
ratio_v=c()
for (i in unique(bubble.df$celltype)) {
  bubble.df_small=bubble.df%>%filter(celltype==i)
  for (j in top10$gene) {
    exp_mean=mean(bubble.df_small[,j])
    exp_ratio=sum(bubble.df_small[,j] > min(bubble.df_small[,j])) / length(bubble.df_small[,j])
    celltype_v=append(celltype_v,i)  ##Add elements to a vector.
    gene_v=append(gene_v,j)
    mean_v=append(mean_v,exp_mean)
    ratio_v=append(ratio_v,exp_ratio)
  }
}

plotdf=data.frame(
  celltype=celltype_v,
  gene=gene_v,
  exp=mean_v,
  ratio=ratio_v
)
library(RColorBrewer)
plotdf$celltype=factor(plotdf$celltype,levels = sort(unique(plotdf$celltype)))
plotdf$gene=factor(plotdf$gene,levels = rev(as.character(top10$gene[1:10])))
plotdf$exp=ifelse(plotdf$exp>3,3,plotdf$exp)
plotdf%>%ggplot(aes(x=celltype,y=gene,size=ratio,color=exp))+geom_point()+
  scale_x_discrete(" ")+scale_y_discrete(" ")+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  scale_size_continuous(limits = c(0,1))+theme_classic()+
  theme(
    axis.text.x.bottom = element_text(hjust =1, vjust = 1, angle = 90)
  )






################################################################################
scRNA1=scRNA2

FeaturePlot(scRNA1,features = "Mzb1",reduction = "tsne",pt.size = 1,split.by = "orig.ident")
FeaturePlot(scRNA1,features = "S100a9",reduction = "tsne",pt.size = 1 )

FeaturePlot(scRNA1,features = "Rbp1",reduction = "tsne",pt.size = 1)+
  scale_x_continuous("")+scale_y_continuous("")+
  theme_bw()+ 
  theme(  
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),  
    axis.ticks = element_blank(),axis.text = element_blank(),  
    legend.position = "none",  
    plot.title = element_text(hjust = 0.5,size=14)  
  )

#### 

mat1=as.data.frame(scRNA1[["RNA"]]@data["Rbp1",])
colnames(mat1)="exp"
mat2=Embeddings(scRNA1,"tsne")
mat3=merge(mat2,mat1,by="row.names")


mat3%>%ggplot(aes(tSNE_1,tSNE_2))+geom_point(aes(color=exp))+
  scale_color_gradient(low = "grey",high ="#1B9E77")+theme_bw()+ 
  theme(  
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),  
    axis.ticks = element_blank(),axis.text = element_blank(),  
    legend.position = "none",  
    plot.title = element_text(hjust = 0.5,size=14)  
  )+
  scale_x_continuous("")+scale_y_continuous("")+labs(title= "Rbp1")




library(ggpubr)

VlnPlot( scRNA1 , "Rbp1",group.by = "orig.ident")+ 
  theme_bw()+stat_compare_means(method = 't.test')+geom_boxplot()

VlnPlot(subset(scRNA1,Rbp1 > 0 ), "Rbp1",group.by = "orig.ident")+ 
  theme_bw()+stat_compare_means(method = 't.test')+geom_boxplot()

sc.FURIN=subset(scRNA1,Rbp1 > 0 )
FURIN.exp=sc.FURIN@assays$RNA@data["Rbp1",]
FURIN.exp.m= sc.FURIN@meta.data
table(sc.FURIN$orig.ident)
FURIN.exp.m$exp=FURIN.exp
FURIN.exp.m$Cohort=factor(FURIN.exp.m$orig.ident,
                          levels = c("DapiNeg1","DapiNeg2","CAF1",'CAF2'))
write.csv(FURIN.exp.m,file = "FURIN.exp.m.csv")
colnames(FURIN.exp.m)
library(ggplot2)
library("ggsignif")
my_comparisons=list(c("DapiNeg1","CAF1"),c("DapiNeg2",'CAF2'),
                    c("DapiNeg2" ,"CAF1"))

ggplot(FURIN.exp.m,aes(Cohort,exp,fill=Cohort)) +
  geom_violin()+theme_classic()+
  theme(plot.title=element_text(size = 15,face="bold"),
        axis.text.x=element_text(size=15,face="bold",angle=25,hjust = 1,vjust = 1),
        axis.text.y=element_text(size=18,face="bold"),
        axis.title.x=element_text(size = 0,vjust = 1,face="bold"),
        axis.title.y=element_text(size = 15,face="bold"))+
  labs(y= 'Expression',x="Group")+
  geom_signif(comparisons = my_comparisons,step_increase = 0.1,
              map_signif_level = F,test = t.test,size=1,textsize = 6)+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  ggtitle("Rbp1")+
  theme(plot.title = element_text(size =18,hjust = 0.5))+
  guides(fill=guide_legend(title = " "))+ ##更改legend名字
  ##如果是color分组 需要 guides(color=guide_legend(title = "Group")) 
  theme(axis.line.x=element_line(linetype=1,color="black",size=1))+
  theme(axis.line.y=element_line(linetype=1,color="black",size=1))+
  theme(legend.text = element_text(size = 15) ) ##legend字体变大


