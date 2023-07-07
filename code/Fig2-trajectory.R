#setwd("C:/study/eye/2022.7/data & code")
setwd("C:/study/eye/2023.3")
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
#create CD4 AND CD8

UMAPPlot(CD4,label=T)
UMAPPlot(CD8,label=T)


Tcell<-FindSubCluster(Tcell,graph.name="integrated_snn",cluster=1,subcluster.name = "subclusterCD8",resolution = 0.5,algorithm = 1)
Idents(Tcell)<-Tcell$subclusterCD8
#Tcellnew<-FindSubCluster(Tcell,graph.name="integrated_snn",cluster=8,subcluster.name = "subcluster",resolution = 0.5,algorithm = 1)
UMAPPlot(Tcell,label=T,group.by="subclusterCD8")

CD4<-subset(Tcell,idents=c(0,2,3,4,5,6,10,11))
CD4 <- FindNeighbors(CD4, dims = 1:30)
CD4 <- RunUMAP(CD4, dims = 1:30)

##re-sort the CD4 cluster order
my_levels <- c(0,2,3,4,5,6,10,11)
Idents(CD4) <- factor(x = Idents(CD4), levels = my_levels)
UMAPPlot(CD4,label=T)+scale_color_d3("category10")


UMAPPlot(CD4,label=T)+ scale_colour_manual(values = cols)
CD4.cds <- as.cell_data_set(CD4)
CD4.cds <- cluster_cells(cds = CD4.cds, reduction_method = "UMAP",resolution=0.00001)#0.001


#CD4.cds <- partitionCells(CD4.cds)

CD4.cds <- learn_graph(CD4.cds, use_partition = TRUE,close_loop=F)
CD4.cds <- order_cells(CD4.cds, reduction_method = "UMAP")

#CD4.cds <- reduce_dimension(CD4.cds)
#plot_cells(CD4.cds, label_groups_by_cluster=FALSE,  color_cells_by = "seurat_clusters")|
plot_cells(CD4.cds,color_cells_by = "ident",show_trajectory_graph = TRUE,label_cell_groups =F,label_groups_by_cluster=F,label_roots = F,label_leaves=F,label_branch_points=F,group_label_size = 5,graph_label_size = 5,trajectory_graph_segment_size=1.5)+ scale_colour_manual(values = cols)|  ##+scale_color_gradient(low = "navyblue",high = "gold")
plot_cells(CD4.cds,color_cells_by = "pseudotime",show_trajectory_graph = TRUE,label_groups_by_cluster=F,label_roots = F,label_leaves=F,label_branch_points=F,group_label_size = 5,graph_label_size = 5,trajectory_graph_segment_size=1.5)


Tcell<-FindSubCluster(Tcell,graph.name="integrated_snn",cluster=1,subcluster.name = "subclusterCD8",resolution = 0.5,algorithm = 1)
Idents(Tcell)<-Tcell$subclusterCD8
CD8<-subset(Tcell,idents=c("1_0","1_1","1_2","1_3","1_4","1_5","1_6","12"))
CD8 <- FindNeighbors(CD8, dims = 1:30)
CD8 <- RunUMAP(CD8, dims = 1:30)

##re-sort the CD8 cluster order
my_levels <- c("1_0","1_1","1_2","1_3","1_4","1_5","1_6","12")
Idents(CD8) <- factor(x = Idents(CD8), levels = my_levels)
UMAPPlot(CD8,label=T)+scale_color_d3("category10")
UMAPPlot(Tcell,label=T,group.by="subclusterCD8")
UMAPPlot(CD4,label=T,label.size=5)+scale_color_d3("category10")|UMAPPlot(CD8,label=T,label.size=5)+scale_color_d3("category10")

#get DGE for CD8
Idents(CD8)<-CD8$subclusterCD8
markersCD8<-FindAllMarkers(CD8,logfc.threshold = 1)
library(dplyr)
markersCD8 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(CD8, features = top10$gene) + NoLegend()
DotPlot(CD8, features = top10$gene[!duplicated(top10$gene)])+coord_flip()
DotPlot(CD8, features = c("CCR7", "SELL", "KLF2", "TCF7", "CCL5", "GZMA", "GZMB", "PRF1", "MKI67", "CD38", "IFNG", "TNF", "PDCD1","CXCR6", "ITGB1","KLRB1","C1QA","C1QB", "CTLA4","FOXP3","IL2RA"),assay = "RNA")+coord_flip()



CD8.cds <- as.cell_data_set(CD8)
CD8.cds <- cluster_cells(cds = CD8.cds, reduction_method = "UMAP",resolution=0.0005)#0.001
CD8.cds <- learn_graph(CD8.cds, use_partition = TRUE,close_loop=F)
CD8.cds <- order_cells(CD8.cds, reduction_method = "UMAP")

plot_cells(CD8.cds,color_cells_by = "ident",show_trajectory_graph = TRUE,label_cell_groups =F,label_groups_by_cluster=F,label_roots = F,label_leaves=F,label_branch_points=F,group_label_size = 5,graph_label_size = 5,trajectory_graph_segment_size=1.5)+ scale_colour_manual(values = cols)|  ##+scale_color_gradient(low = "navyblue",high = "gold")
plot_cells(CD8.cds,color_cells_by = "pseudotime",show_trajectory_graph = TRUE,label_groups_by_cluster=F,label_roots = F,label_leaves=F,label_branch_points=F,group_label_size = 5,graph_label_size = 5,trajectory_graph_segment_size=1.5)


##########################################seek DGE on trajectory(CD4)
library(dplyr)
CD4.cds <- estimate_size_factors(CD4.cds)
Track_genes <- graph_test(CD4.cds,neighbor_graph="principal_graph", cores=6)
Track_genes=cbind(Track_genes,rownames(Track_genes))
colnames(Track_genes)<-c(colnames(Track_genes)[1:5],"gene_short_name")
Track_genes <-Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
##挑选top10画图展示
Track_genes_sig <- Track_genes %>%top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
Track_genes_sig <- Track_genes %>%top_n(n=10, morans_I)  %>% as.character()
##基因表达趋势图

plot_genes_in_pseudotime(CD4.cds[Track_genes_sig,],color_cells_by="seurat_clusters",min_expr=0.5, ncol= 2,label_by_short_name = F)
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

###################################scVelo(unfinished)
library(scuttle)
sce <- as.SingleCellExperiment(Tcellnew)
sce <- logNormCounts(sce, assay.type=1)

library(scran)
dec <- modelGeneVar(sce)
top.hvgs <- getTopHVGs(dec, n=2000)
#We can plug these choices into the scvelo() function with our SingleCellExperiment object. By default, scvelo() uses the steady-state approach to estimate velocities, though the stochastic and dynamical models implemented in scvelo can also be used by modifying the mode argument.

library(velociraptor)
velo.out <- scvelo(sce, subset.row=top.hvgs, assay.X="data")
velo.out

library(scater)

set.seed(100)
sce <- runPCA(sce, subset_row=top.hvgs)
sce <- runTSNE(sce, dimred="PCA")

sce$velocity_pseudotime <- velo.out$velocity_pseudotime
plotTSNE(sce, colour_by="velocity_pseudotime")


plot_cells(CD4.cds,color_cells_by = "ident",show_trajectory_graph = TRUE,label_cell_groups =F,label_groups_by_cluster=F,label_roots = F,label_leaves=F,label_branch_points=F,group_label_size = 5,graph_label_size = 5,trajectory_graph_segment_size=1.5)+ scale_colour_manual(values = cols)|  ##+scale_color_gradient(low = "navyblue",high = "gold")
  plot_cells(CD4.cds,color_cells_by = "pseudotime",show_trajectory_graph = TRUE,label_groups_by_cluster=F,label_roots = F,label_leaves=F,label_branch_points=F,group_label_size = 5,graph_label_size = 5,trajectory_graph_segment_size=1.5)|
  plot_cells(CD8.cds,color_cells_by = "ident",show_trajectory_graph = TRUE,label_cell_groups =F,label_groups_by_cluster=F,label_roots = F,label_leaves=F,label_branch_points=F,group_label_size = 5,graph_label_size = 5,trajectory_graph_segment_size=1.5)+ scale_colour_manual(values = cols)|  ##+scale_color_gradient(low = "navyblue",high = "gold")
  plot_cells(CD8.cds,color_cells_by = "pseudotime",show_trajectory_graph = TRUE,label_groups_by_cluster=F,label_roots = F,label_leaves=F,label_branch_points=F,group_label_size = 5,graph_label_size = 5,trajectory_graph_segment_size=1.5)
