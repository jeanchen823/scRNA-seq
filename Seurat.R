Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())


## install.packages('remotes')
## remotes::install_version("Seurat", version = "3.2.0")
install.packages('Seurat')
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
rm(list=ls())

scRNA.counts=Read10X("GSE152048_BC21.matrix")

scRNA = CreateSeuratObject(scRNA.counts ,
                           min.cells = 3,project="sampl.1", 
                           min.features = 40)


gene=scRNA@assays[["RNA"]]@data@Dimnames[[1]]
bc=scRNA@assays[["RNA"]]@data@Dimnames[[2]]
bc[1:10]
gene[1:100]

phe=scRNA@meta.data
count=scRNA@assays[["RNA"]]@counts
z=scRNA@assays[["RNA"]]@counts@Dimnames[[2]]

x1=scRNA@assays$RNA@data

count=scRNA@assays$RNA@counts
scRNA@meta.data

table(scRNA@meta.data$orig.ident)        

x=scRNA[[]]
x=PercentageFeatureSet(scRNA, pattern = "^MT-")
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
scRNA@meta.data$per.mt=PercentageFeatureSet(scRNA, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
#head(scRNA@meta.data)


a=c("a1","a2","a3")
b=c("b1","a2","m1","m2","a1")
match(a,b)


####Feature、count、线粒体基因、红细胞基因占比可视化。
violin <- VlnPlot(scRNA,
                  features = c("nFeature_RNA", "nCount_RNA",
                               "percent.mt","percent.HB"), 
                  
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) 
violin 

VlnPlot(scRNA,
        features = c("nFeature_RNA", "nCount_RNA",
                     "percent.mt","percent.HB"), 
        
        pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
        ncol = 4)

violin

ggsave("vlnplot_before_qc.pdf",dpi = 300, plot = violin, width = 12, height = 6) 
ggsave("vlnplot_before_qc.png", plot = violin, width = 12, height = 6)  
plot1=FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2=FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3=FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
plot1
plot2
plot3
pearplot

scRNA1 <- subset(scRNA, subset = nFeature_RNA > 500 & 
                   percent.mt < 20 & percent.HB < 1 & nCount_RNA > 1000)
scRNA
scRNA1
scRNA1 <- NormalizeData(scRNA1, normalization.method = "LogNormalize", 
                        scale.factor = 10000)

save(scRNA1,file='scRNA1.Rdata')


scRNA1 <- FindVariableFeatures(scRNA1, selection.method = "vst", 
                               nfeatures = 2000) 

top10 <- head(VariableFeatures(scRNA1), 10) 

plot1 <- VariableFeaturePlot(scRNA1) 

plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5) 
plot <- CombinePlots(plots = list(plot1, plot2),legend="bottom") 

plot 


scale.genes <-  rownames(scRNA1)
scRNA1 <- ScaleData(scRNA1, features = scale.genes)
scale.genes <-  VariableFeatures(scRNA1)
scRNA1 <- ScaleData(scRNA1, features = scale.genes)


scRNA1@assays[["RNA"]]@layers[["counts"]][1:10,1:10]

scRNA1@assays[["RNA"]]@layers[["data"]][1:10,1:10]

scRNA1@assays[["RNA"]]@layers[["scale.data"]][1:10,1:10]


cc.genes
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(scRNA1))
#细胞周期评分
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = VariableFeatures(scRNA1))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = VariableFeatures(scRNA1))
scRNA1 <- CellCycleScoring(object=scRNA1,  g2m.features=g2m_genes,  s.features=s_genes)
scRNAa <- RunPCA(scRNA1, features = c(s_genes, g2m_genes))
p <- DimPlot(scRNAa, reduction = "pca", group.by = "Phase")
p

VlnPlot(scRNAa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB","G2M.Score","S.Score"), ncol = 6)
ggsave("cellcycle_pca.png", p, width = 8, height = 6)

#scRNAb <- ScaleData(scRNA1, vars.to.regress = c("S.Score", "G2M.Score"))



scRNA1 <- RunPCA(scRNA1 , features = VariableFeatures(object = scRNA1) )
plot1 <- DimPlot(scRNA1, reduction = "pca", group.by="orig.ident") 
y1=scRNA1@reductions[["pca"]]@cell.embeddings

ElbowPlot(scRNA1, ndims=50, reduction="pca") 

pc.num=1:20
scRNA1 <- FindNeighbors(scRNA1, dims = 1:20) 
scRNA1 <- FindClusters(scRNA1, resolution = 1)
table(scRNA1$seurat_clusters)

scRNA1<-BuildClusterTree(scRNA1)
PlotClusterTree(scRNA1)


scRNA1 = RunTSNE(scRNA1, dims = 1:20)
embed_tsne <- Embeddings(scRNA1, 'tsne')
write.csv(embed_tsne,'embed_tsne.csv')
plot1 = DimPlot(scRNA1, reduction = "tsne") 
plot1
DimPlot(scRNA1, reduction = "tsne" ) 
ggsave("tSNE.pdf", plot = plot1, width = 8, height = 7)

#UMAP
scRNA1 <- RunUMAP(scRNA1, dims = 1:20)
embed_umap <- Embeddings(scRNA1, 'umap')
write.csv(embed_umap,'embed_umap.csv') 
plot2 = DimPlot(scRNA1, reduction = "umap") 
DimPlot(scRNA1, reduction = "umap",label = T) 

