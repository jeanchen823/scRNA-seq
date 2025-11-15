Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

library(Seurat)
library(scatterplot3d)

load("sc.combined.rdata")

table(sc.combined@meta.data$new.celltype.2)
table(sc.combined@meta.data$type)

table(sc.combined@meta.data$orig.ident)



##  Fibroblasts

sc.imm=subset(sc.combined,idents = c("Epithelial","Endothelial","Stem_cell"))


Idents(sc.imm)="orig.ident"
table( sc.imm@active.ident) 
sc.imm@meta.data$orig.ident=factor(sc.imm@meta.data$orig.ident,levels = c(
  "GF_1","GF_2","GF_3","SPF_1","SPF_2","SPF_3","FMT_1","FMT_2","FMT_3"
))
Idents(sc.imm)="orig.ident"
table(sc.imm@active.ident)
?AverageExpression
x=AverageExpression(sc.imm ,add.ident="new.celltype.2")
x1=AverageExpression(sc.imm  )
x=x$RNA
x=as.matrix(x)
x1=x1$RNA
x1=as.matrix(x1)
pca <- prcomp(t(x))

pca.data=pca$x

rownames(pca.data)

pca.data=as.data.frame(pca.data)
pca.data=pca.data[,1:3]

pca.data$Type = c(rep("GF",9),rep("SPF",9),rep("FMT",9))

s1=strsplit(rownames(pca.data),split = "_",fixed = T)
type=sapply(s1, function(x){x[1]}   )


pca.data$cell.type = c(rep(c("Epithelial","Endothelial","Stem_cell"),9 ))



Type.p=c(rep(11,9),rep(16,9),rep(17,9))
cell.type.p = c(rep(c("blue","red","orange"),9 ))


colors.lib <- c("blue","red","orange")
shapes.lib = c(11,16,17)

getwd()

# 1. Source the function
source( "E:/super.lesson/lesson9/huage.3D.plot.R" )
# 2. 3D scatter plot


s3d <- scatterplot3d(pca.data[,c("PC1","PC2","PC3")],
                     
                     pch = Type.p,color = cell.type.p,
                     
                     cex.symbols = 1,grid=FALSE, box=FALSE, 
                     
                     main = "3D PCA plot")


legend("topright",
       c("Epithelial","Endothelial","Stem_cell"),
       fill=c('blue',"red","orange"),
       box.col=NA)

legend("topleft",
       c('GF','SPF','FMT'),
       pch = shapes.lib,
       box.col=NA)

# 3. Add grids
addgrids3d(pca.data[,c("PC1","PC2","PC3")], grid = c("xy", "xz", "yz"))

table(sc.combined@meta.data$new.celltype.2)

##  Goblet Paneth cells Enteroendocrine
Idents(sc.combined)="new.celltype.3"
sc.imm=subset(sc.combined,idents = c("Goblet","MSC","T_NK"))


Idents(sc.imm)="orig.ident"
x=AverageExpression(sc.imm ,add.ident="new.celltype.3")
x=x$RNA
pca <- prcomp(t(x))

pca.data=pca$x

rownames(pca.data)

pca.data=as.data.frame(pca.data)
pca.data=pca.data[,1:3]

pca.data$Type = c(rep("GF",9),rep("SPF",9),rep("FMT",9))
pca.data$cell.type = c(rep(c("Goblet","Paneth cells","Enteroendocrine"),9 ))



Type.p=c(rep(11,9),rep(16,9),rep(17,9))
cell.type.p = c(rep(c("blue","red","orange"),9 ))


colors.lib <- c("blue","red","orange")
shapes.lib = c(11,16,17)


getwd()


# 1. Source the function
source('plot.3d.R')
# 2. 3D scatter plot


s3d <- scatterplot3d(pca.data[,c("PC1","PC2","PC3")],
                     
                     pch = Type.p,color = cell.type.p,
                     
                     cex.symbols = 1,grid=FALSE, box=FALSE, 
                     
                     main = "3D PCA plot")




legend("topright",
       c("Goblet","Paneth cells","Enteroendocrine"),
       fill=c('blue',"red","orange"),
       box.col=NA)

legend("topleft",
       c('GF','SPF','FMT'),
       pch = shapes.lib,
       box.col=NA)

# 3. Add grids
addgrids3d(pca.data[,c("PC1","PC2","PC3")], grid = c("xy", "xz", "yz"))
