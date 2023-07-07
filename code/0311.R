setwd("C:/study/eye/source")
library(Seurat)
library(scRepertoire)
library(ggsci)
count<-Read10X('BD-2/count/filtered_feature_bc_matrix')
BD2<-CreateSeuratObject(count)

count<-Read10X('BD-3/count/filtered_feature_bc_matrix')
BD3<-CreateSeuratObject(count)
count<-Read10X('VKH-3/count/filtered_feature_bc_matrix')
VKH3<-CreateSeuratObject(count)

BD2 <- NormalizeData(BD2)
BD2 <- FindVariableFeatures(BD2)
BD2 <- ScaleData(BD2)
BD2<- RunPCA(BD2,npcs=50)
ElbowPlot(BD2)
BD2 <- FindNeighbors(BD2, dims = 1:20)
BD2 <- RunUMAP(BD2, dims = 1:20)
BD2 <- FindClusters(BD2)
UMAPPlot(BD2)
FeaturePlot(BD2,features = c("IGHM","CD3E","CD4",'CD8A'))

########################TCR integration
setwd("C:/study/eye/source/tcr/new")

files<-list.files(path = "C:/study/eye/source/tcr/new", pattern = '*.csv')

desired_length <- length(files) 
contig_list <- vector(mode = "list", length = desired_length)
#contig_list[[1]]<-read.csv("BD2.csv")


####strip barcode(use for scRNA-seq data)
combined <- combineTCR(contig_list,samples="BD2",cells ="T-AB")
BD2new <- RenameCells(BD2,add.cell.id = "BD2")
BD2new <- combineExpression(combined, BD2new, cloneCall="aa", proportion = T, cloneTypes=c(Small = 0.002, Medium = 0.01, Large = 0.05,SuperLarge=1))
clonalOverlay(BD2new,reduction = "umap",freq.cutpoint = 0.01, bins = 10,)|FeaturePlot(BD2new,features = c("IGHM","CD3E","CD4",'CD8A'))


#######################################################################BD3
BD3 <- NormalizeData(BD3)
BD3 <- FindVariableFeatures(BD3)
BD3 <- ScaleData(BD3)
BD3<- RunPCA(BD3,npcs=50)
ElbowPlot(BD3)
BD3 <- FindNeighbors(BD3, dims = 1:20)
BD3 <- RunUMAP(BD3, dims = 1:20)
BD3 <- FindClusters(BD3)
UMAPPlot(BD3)
FeaturePlot(BD3,features = c("IGHM","CD3E","CD4",'CD8A'))
########################TCR integration

#desired_length <- length(files) 
contig_list <- vector(mode = "list", length = 1)
contig_list[[1]]<-read.csv("tcr/BD-3.csv")


####strip barcode(use for scRNA-seq data)
combined1 <- combineTCR(contig_list,samples="BD-3",cells ="T-AB")
BD3new <- RenameCells(BD3,add.cell.id = "BD3")
BD3new <- combineExpression(combined1, Tcell, cloneCall="aa", proportion = T, cloneTypes=c(Small = 0.002, Medium = 0.01, Large = 0.05,SuperLarge=1))
clonalOverlay(BD3new,reduction = "umap",freq.cutpoint = 0.01, bins = 10,)|FeaturePlot(BD3new,features = c("IGHM","CD3E","CD4",'CD8A'))
BD3new <- highlightClonotypes(BD3new, 
                               cloneCall= "aa", 
                               sequence = rownames(head(sort(table(BD3new@meta.data[["CTaa"]]),T),n=15)))

DimPlot(BD3new, split.by = "highlight",ncol = 4,label = T) + 
  theme(plot.title = element_blank())

#######################################################################VKH3
VKH3 <- NormalizeData(VKH3)
VKH3 <- FindVariableFeatures(VKH3)
VKH3 <- ScaleData(VKH3)
VKH3<- RunPCA(VKH3,npcs=50)
ElbowPlot(VKH3)
VKH3 <- FindNeighbors(VKH3, dims = 1:20)
VKH3 <- RunUMAP(VKH3, dims = 1:20)
VKH3 <- FindClusters(VKH3)
UMAPPlot(VKH3)
FeaturePlot(VKH3,features = c("IGHM","CD3E","CD4",'CD8A'))
########################TCR integration
setwd("C:/study/eye/source/tcr/new")

files<-list.files(pattern='*.csv')

#desired_length <- length(files) 
contig_list <- vector(mode = "list", length = 1)
contig_list[[1]]<-read.csv("VKH3.csv")


####strip barcode(use for scRNA-seq data)
combined <- combineTCR(contig_list,samples="VKH3",cells ="T-AB")
VKH3new <- RenameCells(VKH3,add.cell.id = "VKH3")
VKH3new <- combineExpression(combined, VKH3new, cloneCall="aa", proportion = T, cloneTypes=c(Small = 0.002, Medium = 0.01, Large = 0.05,SuperLarge=1))
clonalOverlay(VKH3new,reduction = "umap",freq.cutpoint = 0.01, bins = 10,)|FeaturePlot(VKH3new,features = c("IGHM","CD3E","CD4",'CD8A'))
DimPlot(VKH3new, group.by = "cloneType")
VKH3new <- highlightClonotypes(VKH3new, 
                              cloneCall= "aa", 
                              sequence = rownames(head(sort(table(VKH3new@meta.data[["CTaa"]]),T),n=15)))

DimPlot(VKH3new, split.by = "highlight",ncol = 4) + 
  theme(plot.title = element_blank())

#####################investigate 6 samples together
integrated <- readRDS("C:/study/eye/2023.3/integrated.rds")
integrated <- NormalizeData(integrated)
integrated <- FindVariableFeatures(integrated)
integrated <- ScaleData(integrated)
integrated<- RunPCA(integrated,npcs=50)
ElbowPlot(integrated,50)
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- RunUMAP(integrated, dims = 1:30)
integrated <- FindClusters(integrated)

UMAPPlot(integrated,label=T)
UMAPPlot(integrated,split.by="Patient",ncol=3)
table(integrated$seurat_clusters,integrated$Patient)

FeaturePlot(integrated,features = c("CD19",'CD14',"CD3E","CD4",'CD8A','TRDC'),max.cutoff = 'q90',min.cutoff = 'q10')
FeaturePlot(integrated,features = c('CD79A','CD19','IGHM','CD14','FCGR3A','CD3E','TRAC','CD4','CD8A','TRDC','HBA1','CLEC9A', 'CLEC10A', 'FCGR3B','KLRB1','KLRD1'),max.cutoff = 'q90',min.cutoff = 'q10')

VlnPlot(integrated,features = c('CD79A','CD19','IGHM','CD14','FCGR3A','CD3E','TRAC','CD4','CD8A','TRDC','HBA1','CLEC9A', 'CLEC10A', 'FCGR3B','KLRB1','KLRD1'),stack=T)
DotPlot(integrated,features = c('CD79A','CD19','IGHM','CD14','FCGR3A','CD3E','TRAC','CD4','CD8A','TRDC','HBA1','CLEC9A', 'CLEC10A', 'FCGR3B','KLRB1','KLRD1'),assay = "RNA")

saveRDS(integrated,file = "integrated-new.rds")


######################################integrated- subset 14
integrated<-FindSubCluster(integrated,graph.name="integrated_snn",cluster=14,subcluster.name = "subcluster14",resolution = 0.5,algorithm = 1)
UMAPPlot(integrated,label=T,group.by="subcluster14")
UMAPPlot(integrated,label=T,split.by="subcluster14",ncol=5)
DotPlot(integrated,features = c('CD79A','CD19','IGHM','CD14','FCGR3A','CD3E','TRAC','CD4','CD8A','TRDC','TRGC1','HBA1','CLEC9A', 'CLEC10A', 'FCGR3B','KLRB1','KLRD1','MKI67'),group.by = "subcluster14",assay = "RNA")
##################################################################annotation
Tcell<-subset(integrated,idents = c(0,1,2,3,4,8,12,13,19))
unknown<-subset(integrated,idents = c(7,10,11,14,16))
Bcell<-subset(integrated,idents = c(9,17,18,20))
monocyte<-subset(integrated,idents = c(5,6))
DC<-subset(integrated,idents = c(15))
marker14<-FindMarkers(integrated,ident.1 = 14)

dim(Tcell)[2]+dim(unknown)[2]+dim(Bcell)[2]+dim(monocyte)[2]+dim(DC)[2]==dim(integrated)[2]#TRUE



#####cell cycle and regression
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

Tcell <- CellCycleScoring(Tcell, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Tcell <- ScaleData(Tcell, vars.to.regress = c("S.Score", "G2M.Score"))
UMAPPlot(Tcell)
Tcell$oldcluster<-Tcell$seurat_clusters
#Tcell <- NormalizeData(Tcell)
Tcell <- FindVariableFeatures(Tcell)
Tcell <- ScaleData(Tcell)
Tcell<- RunPCA(Tcell,npcs=50)
ElbowPlot(Tcell,50)
Tcell <- FindNeighbors(Tcell, dims = 1:30)
Tcell <- RunUMAP(Tcell, dims = 1:30)
Tcell <- FindClusters(Tcell,resolution = 0.5)

UMAPPlot(Tcell,label=T)|UMAPPlot(Tcell,group.by="Phase",label=T)
  saveRDS(Tcell,"Tcell0322.rds")
#########################remove propiferating Tcells for downstream
Tcell<-subset(Tcell,idents = c(0,1,2,3,4,5,6,8,9,10,11,12))
Tcell <- FindNeighbors(Tcell, dims = 1:30)
Tcell <- RunUMAP(Tcell, dims = 1:30)
##re-sort the Tcell cluster order
my_levels <- c(0,'1_0','1_1','1_2','1_3','1_4','1_5','1_6',2,3,4,5,6,8,9,10,11,12)
Idents(Tcell) <- factor(x = Idents(Tcell), levels = my_levels)
UMAPPlot(Tcell,label=T)+scale_color_d3("category20")

#########################subset CD4 CD8 and re-UMAP for trajectory and TCR
markersT<-FindAllMarkers(Tcell,group.by="subclusterCD8",logfc.threshold = 1)
library(dplyr)
markersT %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10

DoHeatmap(Tcell, features = top10$gene,group.colors=pal_d3('category20')(19),angle = 90,hjust=-0.5) + NoLegend()
DotPlot(Tcell, features = top10$gene[!duplicated(top10$gene)])+coord_flip()

#######integrate TCR
files<-list.files(path = "C:/study/eye/source/tcr/new", pattern = '*.csv')

desired_length <- length(files) 
contig_list <- vector(mode = "list", length = desired_length)
sample<-strsplit(files,split='.csv')
for (i in 1:length(contig_list)) {
  countname <- paste0("C:/study/eye/source/tcr/new/", files[i])
  contig_list[[i]]<-read.csv(countname)
}

####strip barcode(use for scRNA-seq data)
combined <- combineTCR(contig_list,samples=sample,cells ="T-AB")

Tcell <- RenameCells(Tcell,add.cell.id = Tcell$Patient)
#Tcell <- RenameCells(integrated,add.cell.id = integrated$Patient)
Tcell <- RenameCells(Tcell,new.names =substr(colnames(Tcell),start=1,stop = nchar(colnames(Tcell))-2))
Tcell <- combineExpression(combined, Tcell, cloneCall="aa", proportion = T, cloneTypes=c(Rarely = 0.005, Morderately = 0.02, Largely = 1))
Tcell <- combineExpression(combined, Tcell, cloneCall="aa", proportion = T, cloneTypes=c(Rarely-expanded=0.005, Morderately-expanded=0.02, Largely-expanded=1))

#now Tcell is all cells with TCR label
UMAPPlot(Tcell,split.by="cloneType")

###Cluster 9 is ambiguous, remove
Tcell<-subset(Tcell,idents = c(0,1,2,3,4,5,6,7,8,10,11,12))
Tcell <- FindVariableFeatures(Tcell)
Tcell <- ScaleData(Tcell)
Tcell<- RunPCA(Tcell,npcs=50)
ElbowPlot(Tcell,50)
Tcell <- FindNeighbors(Tcell, dims = 1:30)
Tcell <- RunUMAP(Tcell, dims = 1:30)
Tcell <- FindClusters(Tcell,resolution = 0.5)
################################################################################################
################################################################################################DRAW FIGURES
UMAPPlot(Tcell,label=T)|FeaturePlot(Tcell,features = c("CD4"),max.cutoff = 'q90',min.cutoff = 'q10')
FeaturePlot(Tcell,features = c("CD8A","CD8B"),max.cutoff = 'q90',min.cutoff = 'q10',blend=T)
UMAPPlot(Tcell,split.by="Patient",ncol=3)
Idents(Tcell)<-Tcell$seurat_clusters
UMAPPlot(Tcell,split.by="Patient",ncol=3)

table(Tcell$seurat_clusters,Tcell$Patient)
table(Tcell$seurat_clusters,Tcell$Disease)
###draw violin plot(SHJ)
source("C:/study/figure/ShjVlnplot.R")
source("C:/study/figure/ShjVlnplot-nolabel.R")
ShjVlnPlot(Tcell,features = c("TRAC","TRBC1", "TRBC2","TRDC","TRGC1","TRGC2", "CD3D", "CD4",  "CD8A","CD8B","CD19"))
ShjVlnPlot(Tcell,features = c("CCR7", "SELL", "KLF2", "TCF7", "GZMK", "CCL5", "GZMA", "GZMB", "PRF1", "MKI67", "CD38", "IFNG", "TNF", "PDCD1", "CTLA4","FOXP3","IL2RA","CXCR6", "ITGB1","KLRB1",'CD19'))
DotPlot(Tcell,assay = "RNA",group.by="subclustering",features = c('CD3E', 'CD4','CD8A','CD8B','SELL', 'CCR7', 'CD44', 'CD27','CD28','ITGA1', 'CD69', 'ITGAE', 'GZMK','TRAC', 'TRBC1', 'IFNG','FOXP3','BCL6','CXCR5','PDCD1','RORA','IL17RA','IL17F','RORC','KLRB1', 'TRDC', 'TRGC1','MKI67'))
Tcell$Disease<-substr(Tcell$Patient,start=1,stop = nchar(Tcell$Patient)-2)
saveRDS(Tcell,file = "Tcell.rds")

DotPlot(Tcell,assay = "RNA",features = c('CD3E', 'CD4','CD8A','CD8B','SELL', 'CCR7', 'CD44', 'CD27','CD28','ITGA1', 'CD69', 'ITGAE', 'GZMK','TRAC', 'TRBC1', 'IFNG','FOXP3','BCL6','CXCR5','PDCD1','RORA','IL17RA','IL17F','RORC','KLRB1', 'TRDC', 'TRGC1','MKI67'))

markers=c('CD3E', 'CD4','CD8A','CD8B', 'CD69', 'ITGAE','TRAC', 'TRBC1', 'IFNG','FOXP3',
         'SELL', 'CCR7', 'CD44', 'CD27' ,'CD28',
         'BCL6','CXCR5','PDCD1',
         'IL13','IL5','RORA',
        'IL17RA','IL17F','RORC',
        'KLRB1', 'TRDC', 'TRGC1', 'TRDV2', 'TRDV1', 'TRDV3' ,
        'ITGA1','NCAM1','KLRD1',"NCR1","GNLY","KLRC2","KLRC1",
        'ITGAX', 'ITGAM', 'FCN1', 'CD14', 'FCGR3A', 'CD68', 'FCGR1A', "FCGR2B", "C1QC", "C1QA",
        'CD274', 'XCR1', 'CLEC9A', 'CLEC10A', 'HLA-DRA', 'HLA-DRB1',
        "S100A8", "S100A9", "FCGR3B",
        'MKI67')
DotPlot(Tcell,assay = "RNA",features = markers,cluster.idents =F)+coord_flip()

UMAPPlot(Tcell,group.by="Patient")
