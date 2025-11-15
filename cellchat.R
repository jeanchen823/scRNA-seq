Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

library(data.table)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)


load("scRNA_harmony.rdata")

####################################################################################################################################################

install.packages("ggalluvial")
install.packages('NMF')
devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")

devtools::install_github("sqjin/CellChat")

devtools::install_local("CellChat-master.zip") 


library(CellChat)
library(tidyverse)
library(ggalluvial)
library(Seurat)
library(data.table)
library(ggsci)


data.input <- GetAssayData(scRNA_harmony,   slot = "data")
load("data.input.rdata")
identity <- subset(scRNA_harmony@meta.data, select = "celltype")
cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "celltype")


CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

##
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)


unique(CellChatDB$interaction$annotation)
# use Secreted Signaling for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use # set the used database in the object


cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 10)

cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)

cellchat <- filterCommunication(cellchat, min.cells = 3)


df.net.1 <- subsetCommunication(cellchat,slot.name = "netP")
df.net.2 <- subsetCommunication(cellchat )

df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))


cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}



levels(cellchat@idents)            #查看细胞顺序
vertex.receiver = c(1, 2)          #指定靶细胞的索引
cellchat@netP$pathways             #查看富集到的信号通路
pathways.show <- "CCL"            #指定需要展示的通路


##层次图
vertex.receiver = seq(1,2) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = "MK",  
                    vertex.receiver = c(1,2,3),layout="hierarchy")
在层次图中，实体圆和空心圆分别表示源和目标。圆的大小与每个细胞组的细胞数成比例。线越粗，互作信号越强。
左图中间的target是我们选定的靶细胞。右图是选中的靶细胞之外的另外一组放在中间看互作。

##圈图
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling ="MK", layout = "circle")

##和弦图
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling ="MK", vertex.receiver = c(1 ), 
                    layout = "chord", vertex.size = groupSize)

##热图
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = "MK", color.heatmap = "Reds")

pathways.show <- "CCL"  
netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.MK <- extractEnrichedLR(cellchat, signaling = "CCL", geneLR.return = FALSE)

LR.show <- pairLR.MK[2,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,2) # a numeric vector
##层次图
netVisual_individual(cellchat, signaling ="CCL"  ,  pairLR.use = LR.show, vertex.receiver = vertex.receiver,layout="hierarchy")
##圈图
netVisual_individual(cellchat, signaling ="CCL"  , pairLR.use = LR.show, layout = "circle")
##和弦图
netVisual_individual(cellchat, signaling ="CCL"  , pairLR.use = LR.show, layout = "chord")


##############批量保存
pathway.show.all=cellchat@netP$pathways
levels(cellchat@idents)
vertex.receiver=c(1,2,3,4)
for (i in 1:length(pathway.show.all)) {
  
  netVisual(cellchat,signaling = pathway.show.all[i],out.format = c("pdf"),
            vertex.receiver=vertex.receiver,layout="circle")
  plot=netAnalysis_contribution(cellchat,signaling = pathway.show.all[i])
  ggsave(filename = paste0(pathway.show.all[i],".contribution.pdf"),
         plot=plot,width=6,height=4,dpi=300,units="in")
  
}

#####################################################################################
##气泡图
levels(cellchat@idents)
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:5), remove.isolate = FALSE)
##sources.use = 2 是值第二个细胞亚群
netVisual_bubble(cellchat, sources.use =c(1,3), targets.use = c(1:5), remove.isolate = FALSE)
##指定信号通路
cellchat@netP$pathways 
netVisual_bubble(cellchat, sources.use =c(1,3), targets.use =c(1:5),
                 signaling =  c("MK","CCL"), remove.isolate = FALSE)


pairLR  <- extractEnrichedLR(cellchat, signaling =c("MK","CCL"), geneLR.return = FALSE)
netVisual_bubble(cellchat, sources.use =c(1,3), targets.use =c(1:5),pairLR.use =pairLR , remove.isolate = FALSE)

#和弦图               
netVisual_chord_gene(cellchat, sources.use = 2, targets.use = c(1:5), lab.cex = 0.5,legend.pos.y = 30)


plotGeneExpression(cellchat, signaling = "MK")

plotGeneExpression(cellchat, signaling = "MK", enriched.only = FALSE)

plotGeneExpression(cellchat, signaling = "MK",type = "dot")

###################################################################################################
##可视化配体和受体
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = "MK", width = 8, height = 2.5, font.size = 10)
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, 
signaling = c("MK", "PARs"))
gg1 + gg2

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

########################################################################################################
#细胞通讯模式和信号网络

library(NMF)
library(ggalluvial)

selectK(cellchat, pattern = "outgoing")

nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
##river plot
netAnalysis_river(cellchat, pattern = "outgoing")
#气泡图
netAnalysis_dot(cellchat, pattern = "outgoing")
#信号输入细胞的模式识别
selectK(cellchat, pattern = "incoming")

#################################################################################################
#  信号网络聚类
##reticulate::py_install(packages = 'umap-learn')

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)

#############################################################################################
#不同分组之间的配对分析
table(scRNA_harmony@meta.data$orig.ident )
sc.sp=SplitObject(scRNA_harmony,split.by = "orig.ident")
sc.11=scRNA_harmony[,sample(colnames(sc.sp[["sample2"]]),1000)]
sc.3=scRNA_harmony[,sample(colnames(sc.sp[["sample21"]]),1000)]



cellchat.sc11 <- createCellChat(object =sc.11@assays$RNA@data, meta =sc.11@meta.data,  group.by ="celltype")
cellchat.sc3 <- createCellChat(object =sc.3@assays$RNA@data, meta =sc.3@meta.data,  group.by ="celltype")

dir.create("compare")
setwd("compare/")

cellchat=cellchat.sc11 
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.sc11 = cellchat
#################################
cellchat=cellchat.sc3
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.sc3 = cellchat
##############################################
cc.list=list(SC11=cc.sc11,SC3=cc.sc3)
cellchat=mergeCellChat(cc.list,cell.prefix = T,add.names = names(cc.list))
##可视化

compareInteractions(cellchat,show.legend = F,group = c(1,3),measure = "count")
compareInteractions(cellchat,show.legend = F,group = c(1,3),measure = "weight")

netVisual_diffInteraction(cellchat,weight.scale = T)
netVisual_diffInteraction(cellchat,weight.scale = T,measure = "weight")

netVisual_heatmap(cellchat)
netVisual_heatmap(cellchat,measure = "weight")

rankNet(cellchat,mode = "comparison",stacked = T,do.stat = T)
rankNet(cellchat,mode = "comparison",stacked =F,do.stat = T)


weight.max=getMaxWeight(cc.list,attribute = c("idents","count"))
netVisual_circle(cc.list[[1]]@net$count,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc11" )

netVisual_circle(cc.list[[2]]@net$count,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc3" )


table(scRNA_harmony@active.ident)
s.cell=c( "Macrophage", "Tissue_stem_cells","Monocyte")
count1=cc.list[[1]]@net$count[s.cell,s.cell]
count2=cc.list[[2]]@net$count[s.cell,s.cell]

netVisual_circle(count1,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc11" )

netVisual_circle(count2,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc3" )



