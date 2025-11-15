Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(DoubletFinder)
set.seed(1)


data_directory=c("BC2/" ,"BC21/")
project_name=c("sample2", "sample21")


samples <- project_name

sample1 <- make_seurat_object_and_doublet_removal(data_directory[1], samples[1])



###  多个样本合并 
seu_list <- sample1
for (i in 2:length(samples)){
  
  
  
  sc.i = make_seurat_object_and_doublet_removal(data_directory[i], samples[i])
  seu_list=merge(seu_list,sc.i)
  
}

table(seu_list$orig.ident)


scRNA_harmony=seu_list
scRNA_harmony  <- NormalizeData(scRNA_harmony ) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
library(harmony)
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")

scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution =0.5)

scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:20)

DimPlot(scRNA_harmony , reduction = "umap",label = T) 
DimPlot(scRNA_harmony, reduction = "umap", split.by ='orig.ident')

DimPlot(scRNA_harmony, reduction = "umap", group.by='orig.ident')
table(scRNA_harmony$orig.ident)  

save(scRNA_harmony,file = "scRNA_harmony.rdata")


##############################################################################
###  锚定整合

ifnb=seu_list
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



############################################################
############################################################
scRNA_harmony=JoinLayers(scRNA_harmony)

options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
BiocManager::install("MAST")
library(MAST)

table(scRNA_harmony@meta.data$seurat_clusters)

sc.s=sample(colnames(scRNA_harmony),)

sc.s=subset(scRNA_harmony,downsample=100)
table(sc.s@meta.data$seurat_clusters)
markers <- FindAllMarkers(sc.s, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#Macrophage
#### 2  C1QA  C1QB  
scRNA_harmony=RenameIdents(scRNA_harmony,"2"="Macrophage")

scRNA_harmony=RenameIdents(scRNA_harmony,"7"="Endothelial cell","4"="Fibroblast",
                           "1"="B cell","5"="Myeloid","3"="NK")

table(scRNA_harmony@active.ident)
DimPlot(scRNA_harmony,label = T)


########################################################
##  方法二

table(scRNA_harmony$seurat_clusters)

Idents(scRNA_harmony)="seurat_clusters"

DotPlot(scRNA_harmony,features =c("LYZ","NKG7","CD3E",'CD3D','VCAN',
                                  'CD19', 'CD79A',"PECAM1" ))


scRNA_harmony=RenameIdents(scRNA_harmony,"7"="Endothelial cell","4"="Fibroblast",
                           "1"="B cell" ,"0"="Myeloid","8"="Myeloid","6"="Myeloid",
                           "2"="Myeloid", "5"="Myeloid","3"="NK_T")


DimPlot(scRNA_harmony,label = T)

genes_to_check <- list(Tcells = c("CD2","CD3D","CD3E","NKG7"),
                       Bcells = c("CD79A","MS4A1"),
                       SMC = c("ACTA2","CNN1","PLN"),
                       MonoMacro = c("CD14","CD68","CD163"),
                       Dendriticcells = c("CD1C","CLEC9A","CLEC10A","LILR4"),
                       Mastcells = c("TPSAB1","TPSB2","CPA3"),
                       Fibroblasts = c("PDGFRA","FN1","DCN","LUM"),
                       Pericyte = c("MYH11","RGS5","MCAM","NOTCH3"),
                       Endothelialcells = c("PECAM1","VWF","PLVAP","CD34"),
                       Neutrophils = c('CSF3R','FCGR3B','CXCL8')) 


p_all_markers <- DotPlot(scRNA_harmony, features =genes_to_check, assay = "RNA") + 
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1),axis.title = element_blank()) 
p_all_markers

##############################################################################

sc1=subset(scRNA_harmony,ident=c(1:10))
DimPlot(sc1)
table(scRNA_harmony$seurat_clusters)
#第三种方法用SingleR鉴定细胞类型
BiocManager::install("SingleR")

BiocManager::install("celldex")

library(SingleR)

HumanPrimaryCellAtlasData=celldex::HumanPrimaryCellAtlasData()
save(HumanPrimaryCellAtlasData,file = "HumanPrimaryCellAtlasData.rdata")
BlueprintEncodeData=celldex::BlueprintEncodeData()
save(BlueprintEncodeData,file = "BlueprintEncodeData.rdata")
DatabaseImmuneCellExpressionData=celldex::DatabaseImmuneCellExpressionData()
save(DatabaseImmuneCellExpressionData,file = "DatabaseImmuneCellExpressionData.rdata")
ImmGenData=celldex::ImmGenData()
save(ImmGenData,file = "ImmGenData.rdata")

MonacoImmuneData=celldex::MonacoImmuneData()
save(MonacoImmuneData,file = "MonacoImmuneData.rdata")


MouseRNAseqData=celldex::MouseRNAseqData()
save(MouseRNAseqData,file = "MouseRNAseqData.rdata")

NovershternHematopoieticData=celldex::NovershternHematopoieticData()
save(NovershternHematopoieticData,file = "NovershternHematopoieticData.rdata")


load("ref_Human_all.RData")

refdata <- HumanPrimaryCellAtlasData

testdata <- GetAssayData(scRNA_harmony, slot="data")

clusters <- scRNA_harmony@meta.data$seurat_clusters

table(refdata$label.main)

cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels)

write.csv(celltype,"celltype_singleR.csv",row.names = FALSE)
##把singler的注释写到metadata中 有两种方法
###方法一
scRNA_harmony@meta.data$celltype ="NA"
for(i in 1:nrow(celltype)){
  scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}


###  第一次循环  
i = 1
celltype$celltype[1]


celltype$ClusterID[1]
which(scRNA_harmony@meta.data$seurat_clusters == celltype$ClusterID[i])

which(scRNA_harmony@meta.data$seurat_clusters ==  "0")
a=c("a1","a2","a3")
b=c("b1","a2","m1","m2","a1")
which(b=="a1")


DimPlot(scRNA_harmony, group.by="celltype", label=T, label.size=5)
###方法二：
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels ) 
scRNA_harmony@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
DimPlot(scRNA_harmony, group.by="singleR", label=T, label.size=5)

##鉴定结果展示

DimPlot(scRNA_harmony, group.by="celltype", label=T, label.size=5, reduction='umap')


##############  
library(SingleR)
library(celldex)
library(Seurat)

# read in seurat object

colon <-  scRNA_harmony

# set up reference, and define cell types to use
ref <- celldex::HumanPrimaryCellAtlasData()
types_to_use <- c("DC","Epithelial_cells","B_cell","Neutrophils","T_cells","Monocyte",
                  "Endothelial_cells","Neurons","Macrophage","NK_cell",
                  "BM","Platelets","Fibroblasts","Astrocyte","Myelocyte","Pre-B_cell_CD34-","Pro-B_cell_CD34+","Pro-Myelocyte")
ref <- ref[,(colData(ref)$label.main %in% types_to_use)]

# Run singleR
singler.pred <- SingleR(test = as.SingleCellExperiment(colon), ref = ref, labels = ref$label.fine)

# Add to seurat object, and then you cna plot the results if desired
colon <- AddMetaData(colon, metadata = singler.pred$labels, col.name = "SingleR.labels")














