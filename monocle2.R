Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
load("scRNA_harmony.rdata")
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)

table(scRNA_harmony@meta.data$celltype)
table(scRNA_harmony@meta.data$seurat_clusters)
sc.t=scRNA_harmony[,rownames(subset(scRNA_harmony@meta.data,celltype=="T_cells"))]  

options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")
#BiocManager::install("monocle")
library(monocle)
Idents(scRNA_harmony)="celltype"
scRNA.Osteoclastic=subset(scRNA_harmony,ident=c("T_cells"))

##
Idents(scRNA.Osteoclastic)="orig.ident"
sc.sample21=subset(scRNA.Osteoclastic,ident="sample21")
table(scRNA.Osteoclastic$orig.ident)

######
data <- GetAssayData(scRNA.Osteoclastic,slot = "count")
data[1:20,1:20]

pd <- new('AnnotatedDataFrame', data = scRNA.Osteoclastic@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())



monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
print(head(fData(monocle_cds)))

HSMM=monocle_cds
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, disp.genes)
plot_ordering_genes(HSMM)


HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')


HSMM <- orderCells(HSMM)
plot_cell_trajectory(HSMM, color_by = "seurat_clusters")

plot_cell_trajectory(HSMM, color_by = "celltype")

plot_cell_trajectory(HSMM, color_by = "orig.ident")
plot_cell_trajectory(HSMM, color_by = "State")




plot_cell_trajectory(HSMM, color_by = "Pseudotime")

HSMM <- orderCells(HSMM,root_state = 5)
plot_cell_trajectory(HSMM, color_by = "Pseudotime")



plot_cell_trajectory(HSMM, color_by = "celltype")
plot_cell_trajectory(HSMM, color_by = "State") +
  facet_wrap(~State, nrow = 3)

plot_cell_trajectory(HSMM, color_by = "State") +
  facet_wrap(~orig.ident )

blast_genes <- row.names(subset(fData(HSMM),
                                gene_short_name %in% c("GAPDH", "RORA")))
plot_genes_jitter(HSMM[blast_genes,],
                  grouping = "State",
                  min_expr = 0.1)



HSMM_expressed_genes <-  row.names(subset(fData(HSMM),
                                          num_cells_expressed >= 10))
HSMM_filtered <- HSMM[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("YWHAB", "GAPDH", "TNNC1")))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "State")

plot_genes_in_pseudotime(cds_subset, color_by =  "celltype")



genes <- c("TNNT2", "TNNC1", "CDK1")
p1 <- plot_genes_jitter(HSMM[genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(HSMM[genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(HSMM[genes,], color_by = "State")
plotc <- p1|p2|p3
plotc 

to_be_tested <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% c("MYH3", "MEF2C", "CCNB2", "TNNT1")))
cds_subset <- HSMM[to_be_tested,]

diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res[,c("gene_short_name", "pval", "qval")]
plot_genes_in_pseudotime(cds_subset, color_by ="State")




marker_genes <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% c("MEF2C", "MEF2D", "MYF5",
                                                        "ANPEP", "PDGFRA","MYOG",
                                                        "TPM1",  "TPM2",  "MYH2",
                                                        "MYH3",  "NCAM1", "TNNT1",
                                                        "TNNT2", "TNNC1", "CDK1",
                                                        "CDK2",  "CCNB1", "CCNB2",
                                                        "CCND1", "CCNA1", "ID1")))

diff_test_res <- differentialGeneTest(HSMM[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 1))
plot_pseudotime_heatmap(HSMM[sig_gene_names,],
                        num_clusters = 4,
                        cores = 1,
                        show_rownames = T)






plot_cell_trajectory(HSMM, color_by = "State")


BEAM_res <- BEAM(HSMM, branch_point = 1, cores = 7)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res,
                                                  qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)


genes <- row.names(subset(fData(HSMM),
                          gene_short_name %in% c( "MEF2C", "CCNB2", "TNNT1")))

plot_genes_branched_pseudotime(HSMM[genes,],
                               branch_point = 1,
                               color_by = "State",
                               ncol = 1)




