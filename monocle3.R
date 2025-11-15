Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(monocle3)

load("scRNA_harmony.rdata")


data <- GetAssayData(scRNA_harmony,  slot = 'counts')
cell_metadata <- scRNA_harmony@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)

plot_cells(cds)
DimPlot(scRNA_harmony)
plot_cells(cds, color_cells_by="seurat_clusters")


plot_cells(cds, genes=c("GAPDH"))

?reduce_dimension
cds <- reduce_dimension(cds,umap.fast_sgd=T,cores = 3,reduction_method="tSNE")

plot_cells(cds, reduction_method="tSNE")

plot_cells(cds, reduction_method="tSNE", color_cells_by="celltype")


plot_cells(cds, color_cells_by="orig.ident", label_cell_groups=FALSE,reduction_method="tSNE")

cds = align_cds(cds, num_dim = 100, alignment_group ="orig.ident")
cds = reduce_dimension(cds,reduction_method="tSNE")
plot_cells(cds, color_cells_by="orig.ident", label_cell_groups=FALSE,,reduction_method="tSNE")

cds = reduce_dimension(cds )


cds = cluster_cells(cds, resolution=1e-5)
plot_cells(cds)

table(cds@clusters$UMAP$clusters)
cds@clusters$UMAP$partitions
plot_cells(cds, color_cells_by="partition", group_cells_by="partition")

plot_cells(cds, color_cells_by="celltype")

marker_test_res <- top_markers(cds, group_cells_by="partition", 
                               reference_cells=1000, cores=8)


top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="celltype",
                    ordering_type="maximal_on_diag",
                    max.size=3)





top_specific_markers = marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="celltype",
                    ordering_type="cluster_row_col",
                    max.size=3)

colData(cds)$assigned_cell_type <- as.character(partitions(cds))

cds_subset <- choose_cells(cds)


cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "celltype",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           group_label_size=4,cell_size=0.5)


cds = order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=F,
           graph_label_size=1.5,
           group_label_size=4,cell_size=0.5)

table(scRNA_harmony@meta.data$celltype)
# a helper function to identify the root principal points:
get_earliest_principal_node  <- function(cds, time_bin="T_cells"){
  cell_ids <- which(colData(cds)[, "celltype"] == time_bin)
  
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1.5)


cds_sub <- choose_graph_segments(cds)




#Working with 3D trajectories

cds_3d = reduce_dimension(cds, max_components = 3)
cds_3d = cluster_cells(cds_3d)
cds_3d = learn_graph(cds_3d)
cds_3d = order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))

cds_3d = order_cells(cds_3d)

cds_3d_plot_obj = plot_cells_3d(cds_3d, color_cells_by="celltype")
cds_3d_plot_obj


ciliated_genes = top_specific_markers$gene_id[5:10]
cds_subset = cds[rowData(cds)$gene_short_name %in% ciliated_genes,]

gene_fits = fit_models(cds_subset, model_formula_str = "~seurat_clusters")
fit_coefs
fit_coefs = coefficient_table(gene_fits)

emb_time_terms = fit_coefs %>% filter(term == "seurat_clusters1")
emb_time_terms

emb_time_terms = fit_coefs %>% filter(term == "seurat_clusters1")

emb_time_terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

plot_genes_violin(cds_subset[,], group_cells_by="celltype", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


neurons_cds <- cds[,grepl("Macrophage", colData(cds)$celltype, ignore.case=TRUE)]
plot_cells(neurons_cds, color_cells_by="partition")

pr_graph_test_res <- graph_test(neurons_cds, neighbor_graph="knn", cores=3)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

gene_module_df = find_gene_modules(neurons_cds[pr_deg_ids,], resolution=1e-2)

cell_group_df = tibble::tibble(cell=row.names(colData(neurons_cds)), cell_group=neurons_cds@colData@listData[["seurat_clusters"]])
agg_mat = aggregate_gene_expression(neurons_cds, gene_module_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) = stringr::str_c("Partition ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)


plot_cells(neurons_cds,
           genes=gene_module_df %>% filter(module %in% c(16,38,33,42)),
           group_cells_by="cluster",
           color_cells_by="cluster",
           show_trajectory_graph=FALSE)


plot_cells(cds,
           color_cells_by = "celltype",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)



trace('calculateLW', edit = T, where = asNamespace("monocle3"))

ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)



pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))

AFD_genes <-pr_deg_ids[20:22]
AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,
]
plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="celltype",
                         min_expr=0.5)


cds_subset <- choose_cells(cds)

subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))

gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=0.001)

agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, 
                                levels = row.names(agg_mat)[module_dendro$order])

plot_cells(cds_subset,
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)




