setwd("C:/study/eye/2022.7/data & code")
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)

load("C:/study/eye/2022.7/data & code/Tcell0423.rdata")
load("C:/study/eye/2022.7/data & code/CD4&CD8.rdata")
# Load the data

CD4<-subset(CD4,idents=c(1,2,4,5,8,9,10))###remore C0 and C7
CD4.cds <- as.cell_data_set(CD4)
CD4.cds <- cluster_cells(cds = CD4.cds, reduction_method = "UMAP",resolution=0.0005)
CD4.cds <- learn_graph(CD4.cds, use_partition = TRUE,close_loop=F)
CD4.cds <- order_cells(CD4.cds, reduction_method = "UMAP")

CD8.cds <- as.cell_data_set(CD8)
CD8.cds <- cluster_cells(cds = CD8.cds, reduction_method = "UMAP",resolution=0.0005)
CD8.cds <- learn_graph(CD8.cds, use_partition = TRUE,close_loop=F)
CD8.cds <- order_cells(CD8.cds, reduction_method = "UMAP")

# plot trajectories colored by pseudotime
plot_cells(
  cds = CD4.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

plot_cells(
  cds = CD8.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

###all Tcells
Tcell.cds <- as.cell_data_set(Tcell)
Tcell.cds <- cluster_cells(cds = Tcell.cds, reduction_method = "UMAP",resolution=0.0001)
Tcell.cds <- learn_graph(Tcell.cds, use_partition = TRUE,close_loop=F)
Tcell.cds <- order_cells(Tcell.cds, reduction_method = "UMAP")


plot_cells(
  cds = Tcell.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_groups_by_cluster=F,
  label_leaves=T,
  label_branch_points=T,
  group_label_size = 5,
  graph_label_size = 5
)


plot_cells(
  cds = Tcell.cds,
  color_cells_by = "ident",
  show_trajectory_graph = TRUE,
  label_groups_by_cluster=F,
  label_roots = T,
  label_leaves=T,
  label_branch_points=T,
  group_label_size = 0,
  graph_label_size = 5
)

#####CD4
pdf("C:/study/eye/2022.7/figure/Fig2/Trajectory-CD4-psudotime.pdf",6,4)
plot_cells(
  cds = CD4.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_groups_by_cluster=F,
  label_leaves=T,
  label_branch_points=F,
  group_label_size = 8,
  graph_label_size = 6
)
dev.off()

pdf("C:/study/eye/2022.7/figure/Fig2/Trajectory-CD4-cluster.pdf",5,4)
plot_cells(
  cds = CD4.cds,
  color_cells_by = "ident",
  show_trajectory_graph = TRUE,
  label_groups_by_cluster=F,
  label_leaves=T,
  label_branch_points=F,
  group_label_size = 8,
  graph_label_size = 6

)
dev.off()
#####CD8

pdf("C:/study/eye/2022.7/figure/Fig2/Trajectory-CD8-psudotime.pdf",6,4)
plot_cells(
  cds = CD8.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_groups_by_cluster=F,
  label_leaves=T,
  label_branch_points=F,
  group_label_size = 8,
  graph_label_size = 6
)
dev.off()
pdf("C:/study/eye/2022.7/figure/Fig2/Trajectory-CD8-cluster.pdf",5,4)
plot_cells(
  cds = CD8.cds,
  color_cells_by = "ident",
  show_trajectory_graph = TRUE,
  label_groups_by_cluster=F,
  label_leaves=T,
  label_branch_points=F,
  group_label_size = 8,
  graph_label_size = 6
)
dev.off()



##########################################seek DGE on trajectory(CD4)
library(dplyr)
CD4.cds <- estimate_size_factors(CD4.cds)
Track_genes <- graph_test(CD4.cds,neighbor_graph="principal_graph")
Track_genes <-Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
##挑选top10画图展示
Track_genes_sig <- Track_genes %>%top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
##基因表达趋势图

plot_genes_in_pseudotime(CD4.cds[Track_genes_sig,],color_cells_by="seurat_clusters",min_expr=0.5, ncol= 2)
###FeaturePlot图

p <- plot_cells(CD4.cds,genes=Track_genes_sig, show_trajectory_graph=FALSE,
                
                label_cell_groups=FALSE, label_leaves=FALSE)

p$facet$params$ncol <- 5
p

##寻找共表达基因模块
genelist <- pull(Track_genes,gene_short_name) %>% as.character()

gene_module <-find_gene_modules(CD4.cds[genelist,], resolution=1e-1, cores = 6)

write.csv(gene_module,"Genes_Module.csv", row.names = F)

cell_group <-tibble::tibble(cell=row.names(colData(CD4.cds)),
                            
                            cell_group=colData(CD4.cds)$seurat_clusters)

agg_mat <-aggregate_gene_expression(CD4.cds, gene_module, cell_group)

row.names(agg_mat) <-stringr::str_c("Module ", row.names(agg_mat))

p <- pheatmap::pheatmap(agg_mat,scale="column", clustering_method="ward.D2")

##########################################seek DGE on trajectory(CD8)
library(dplyr)
CD8.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(CD8)
CD8.cds <- estimate_size_factors(CD8.cds)
Track_genes <- graph_test(CD8.cds,neighbor_graph="principal_graph")
Track_genes <-Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
##挑选top10画图展示
Track_genes_sig <- Track_genes %>%top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
##基因表达趋势图

plot_genes_in_pseudotime(CD8.cds[Track_genes_sig,],color_cells_by="seurat_clusters",min_expr=0.5, ncol= 2)
###FeaturePlot图

p <- plot_cells(CD8.cds,genes=Track_genes_sig, show_trajectory_graph=FALSE,
                
                label_cell_groups=FALSE,?? label_leaves=FALSE)

p$facet$params$ncol <- 5
p

##寻找共表达基因模块
genelist <- pull(Track_genes,gene_short_name) %>% as.character()

gene_module <-find_gene_modules(CD8.cds[genelist,], resolution=1e-1, cores = 6)

write.csv(gene_module,"Genes_Module.csv", row.names = F)

cell_group <-tibble::tibble(cell=row.names(colData(CD8.cds)),
                            
                            cell_group=colData(CD8.cds)$seurat_clusters)

agg_mat <-aggregate_gene_expression(CD8.cds, gene_module, cell_group)

row.names(agg_mat) <-stringr::str_c("Module ", row.names(agg_mat))

p <- pheatmap::pheatmap(agg_mat,scale="column", clustering_method="ward.D2")


##########################################seek DGE on trajectory(Tcell)
library(dplyr)
Tcell.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(Tcell)
Tcell.cds <- estimate_size_factors(Tcell.cds)
Track_genes <- graph_test(Tcell.cds,neighbor_graph="principal_graph")
Track_genes <-Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
##挑选top10画图展示
Track_genes_sig <- Track_genes %>%top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
##基因表达趋势图

plot_genes_in_pseudotime(Tcell.cds[Track_genes_sig,],color_cells_by="seurat_clusters",min_expr=0.5, ncol= 2)
###FeaturePlot图

p <- plot_cells(Tcell.cds,genes=Track_genes_sig, show_trajectory_graph=FALSE,
                
                label_cell_groups=FALSE,?? label_leaves=FALSE)

p$facet$params$ncol <- 5
p

##寻找共表达基因模块
genelist <- pull(Track_genes,gene_short_name) %>% as.character()

gene_module <-find_gene_modules(Tcell.cds[genelist,], resolution=1e-1)

write.csv(gene_module,"Genes_Module.csv", row.names = F)

cell_group <-tibble::tibble(cell=row.names(colData(Tcell.cds)),
                            
                            cell_group=colData(Tcell.cds)$seurat_clusters)

agg_mat <-aggregate_gene_expression(Tcell.cds, gene_module, cell_group)

row.names(agg_mat) <-stringr::str_c("Module ", row.names(agg_mat))

p <- pheatmap::pheatmap(agg_mat,scale="column", clustering_method="ward.D2")
