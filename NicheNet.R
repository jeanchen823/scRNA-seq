#代码 ：  https://github.com/ruoyan-li/RCC-spatial-mapping/blob/main/code/NichNet_Mac_RCC.R
#文章：  https://www.cell.com/cancer-cell/fulltext/S1535-6108(22)00548-7#sectitle0030 
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)

load("scRNA_harmony.rdata")

####################################################################################################################################################
##library(devtools)
##install_github("saeyslab/nichenetr")
library(nichenetr)
ligand_target_matrix = readRDS("ligand_target_matrix.rds")
ligand_target_matrix[1:5,1:5]

lr_network = readRDS("lr_network.rds")
head(lr_network)

weighted_networks = readRDS("weighted_networks.rds")
head(weighted_networks$lr_sig) 
scRNA_harmony@meta.data$celltype %>% table()

scRNA_harmony@meta.data$orig.ident %>% table()
DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident")


Idents(scRNA_harmony) <- "celltype"


sender_celltypes = c(  "Chondrocytes"  )

nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = scRNA_harmony, 
  receiver = "Macrophage", 
  condition_colname = "orig.ident", condition_oi = "sample2", condition_reference = "sample21", 
  sender = sender_celltypes, 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks )

nichenet_output %>% names()

nichenet_output$ligand_activities

nichenet_output$top_ligands


DotPlot(scRNA_harmony, features = nichenet_output$top_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

DotPlot(scRNA_harmony, features = nichenet_output$top_ligands %>% rev(), split.by = "orig.ident") + RotatedAxis()

VlnPlot(scRNA_harmony, features = c( "CTGF" ,"INHBA","LAMB2","COL4A1" ), split.by = "orig.ident", pt.size = 0, combine = T)

nichenet_output$ligand_target_heatmap

nichenet_output$ligand_target_heatmap + 
  scale_fill_gradient2(low = "whitesmoke",  high = "royalblue", breaks = c(0,0.0045,0.009)) + 
  xlab("Macrophage") + 
  ylab("Prioritized  cell ligands")

ligand.target.matrix <- nichenet_output$ligand_target_matrix


DotPlot(scRNA_harmony %>% subset(idents = "Macrophage"), 
        features = nichenet_output$top_targets, 
        split.by = "orig.ident") + RotatedAxis()

VlnPlot(scRNA_harmony %>% subset(idents = "Macrophage"), features = nichenet_output$top_targets[1:5], 
        split.by = "orig.ident", pt.size = 0, combine = T, ncol = 8)

nichenet_output$ligand_receptor_heatmap

ligand.receptor.matrix=nichenet_output$ligand_receptor_matrix


DotPlot(scRNA_harmony %>% subset(idents = "Macrophage"), 
        features = nichenet_output$top_receptors, 
        split.by = "orig.ident") + RotatedAxis()

VlnPlot(scRNA_harmony %>% subset(idents = "Macrophage"),  features = nichenet_output$top_receptors[1:5], 
        split.by = "orig.ident",pt.size = 0, combine = T, ncol = 8)


################################################################
