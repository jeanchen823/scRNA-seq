Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

library(data.table)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
#install.packages("arrow")
#BiocManager::install(c("AUCell", "RcisTarget"))
#BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost
#BiocManager::install(c("zoo", "mixtools", "rbokeh"))
#BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
#BiocManager::install(c("doMC", "doRNG"))
#devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
#devtools::install_github("aertslab/SCENIC") 

library(patchwork)
library(SCENIC)
library(harmony)

load("sc.m.rdata")
Idents(sc.epi)="celltype.3"
table(sc.epi$celltype.3)

sc.s=subset(sc.epi,downsample=150)
save(sc.s,file = "sc.s.rdata") 
####################################################################################################################################################


scRNAsub= sc.s

exprMat <- as.matrix(scRNAsub@assays$RNA@data)

##https://resources.aertslab.org/cistarget/

data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

mydbDIR <- "F:/ref.data/"
mydbs <- c( "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather","hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
names(mydbs) <- c("10kb","500bp")
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=30,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "m")
saveRDS(scenicOptions, "int/scenicOptions.rds")

##基因过滤

genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
table(genesKept)
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)

runCorrelation(exprMat_filtered, scenicOptions)

exprMat_filtered_log <- log2(exprMat_filtered+1)

runGenie3(exprMat_filtered, scenicOptions, nParts = 20)
head("int/1.4_GENIE3_linkList.Rds")
head(readRDS("int/1.4_GENIE3_linkList.Rds") )

scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 123

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions) #1. 获取共表达模块
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)  #2. 获取regulons
?runSCENIC_2_createRegulons

##==regulon活性评分与可视化==##

load("scRNA_harmony.Rdata")
mydbDIR <- "D:/ref.data/"
mydbs <- c( "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather","hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
names(mydbs) <- c("10kb","500bp")

scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=1,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "m")
library(foreach)
exprMat_all <- as.matrix(scRNA_harmony@assays$RNA@counts)
exprMat_all <- log2(exprMat_all+1)
rm(scRNA_harmony)
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)


runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all) 

cellInfo <- data.frame(scRNA_harmony@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
colnames(cellInfo)[which(colnames(cellInfo)=="celltype_Monaco")] <- "celltype"
cellInfo <- cellInfo[,c("sample","cluster","celltype")]
saveRDS(cellInfo, file="int/cellInfo.Rds")

mydbDIR <- "D:/ref.data/"
mydbs <- c( "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather","hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
names(mydbs) <- c("10kb","500bp")

scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=1,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "os")


library(foreach)
nPcs <- c(5)

fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))

# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))

par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="celltype",   cex=.5)




tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")


library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)

AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat_all, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("JUN","MYC")],], plots="Expression")

#par(bg = "black")
par(mfrow=c(1,2))
regulonNames <- c( "JUN","MYC")
cellCol <- plotEmb_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)


regulonNames <- list( green=c("JUN"),
                      blue=c( "MYC"))
cellCol <- plotEmb_rgb(scenicOptions, regulonNames, aucType="Binary")

regulons <- loadInt(scenicOptions, "regulons")
regulons[c( "JUN","MYC")]

regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="YY1" & highConfAnnot==TRUE]
viewMotifs(tableSubset, options=list(pageLength=5)) 

motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="YY1"]
viewMotifs(tableSubset)



regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$celltype),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[1:20,], name="Regulon activity")


topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)

minPerc <- .7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$celltype), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
ComplexHeatmap::Heatmap(binaryActPerc_subset[1:20,], name="Regulon activity (%)", col = c("white","pink","red"))


topRegulators <- reshape2::melt(regulonActivity_byCellType_Binarized)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>minPerc),]
viewTable(topRegulators)



rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "celltype"])
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)


plotRSS_oneSet(rss, setName = "MSC")


library(Seurat)
#scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:16)
dr_coords <- Embeddings(scRNA_harmony, reduction="tsne")
par(mfrow=c(4,4))
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, "MYC"), plots = "AUC")


tfs <- c("HDAC2","RAD21","YY1", "SMARCA4")
par(mfrow=c(2,2))
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, tfs), plots = "AUC")
########################################################################


AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(scRNA_harmony, AUCmatrix)
scRNAauc@assays$integrated <- NULL
saveRDS(scRNAauc,'scRNAauc.rds')

BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(scRNA_harmony, BINmatrix)
scRNAbin@assays$integrated <- NULL
saveRDS(scRNAbin, 'scRNAbin.rds')

dir.create('scenic_seurat')
#FeaturePlot
colnames(scRNAauc@meta.data)[20:30]
p1 = FeaturePlot(scRNAauc, features="HCFC1_24g", label=T, reduction = 'tsne')
p2 = FeaturePlot(scRNAbin, features="HCFC1_24g", label=T, reduction = 'tsne')
p3 = DimPlot(scRNA_harmony, reduction = 'tsne', group.by = "celltype", label=T)
plotc = p1|p2|p3
plotc



#RidgePlot&VlnPlot
p1 = RidgePlot(scRNAauc, features = "NR3C1_1339g", group.by="celltype") + 
  theme(legend.position='none')
p2 = VlnPlot(scRNAauc, features ="NR3C1_1339g", pt.size = 0, group.by="celltype") + 
  theme(legend.position='none')
plotc = p1 + p2
plotc


library(pheatmap)
cellInfo <- readRDS("int/cellInfo.Rds")
celltype = subset(cellInfo,select = 'celltype')
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)
my.regulons <- rownames(AUCmatrix)[80:100]
myAUCmatrix <- AUCmatrix[rownames(AUCmatrix)%in%my.regulons,]
myBINmatrix <- BINmatrix[rownames(BINmatrix)%in%my.regulons,]
pheatmap(myAUCmatrix, show_colnames=F, annotation_col=celltype )
pheatmap(myBINmatrix, show_colnames=F, annotation_col=celltype,
         color = colorRampPalette(colors = c("white","black"))(100))


