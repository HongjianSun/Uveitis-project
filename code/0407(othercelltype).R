setwd("C:/study/eye/source")
library(Seurat)
library(scRepertoire)

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

source("C:/study/figure/ShjVlnplot-nolabel.R")
ShjVlnPlot(integrated,features = c('CD79A','CD19','IGHM','IGHG1','IGHD','CD14','FCGR3A','CD3E','TRAC','CD4','CD8A','TRDC','HBA1','CLEC9A', 'CLEC10A', 'FCGR3B','KLRB1','KLRD1','NCAM1'))
DotPlot(integrated,assay="RNA",features = c('CD79A','CD19','IGHM','IGHG1','IGHD','CD14','FCGR3A','CD3E','TRAC','CD4','CD8A','TRDC','HBA1','CLEC9A', 'CLEC10A', 'FCGR3B','KLRB1','KLRD1','NCAM1'))

FeaturePlot(integrated,features = c("CD19",'CD14',"CD3E","CD4",'CD8A','TRDC'),max.cutoff = 'q90',min.cutoff = 'q10')
FeaturePlot(integrated,features = c('CD79A','CD19','IGHM','CD14','FCGR3A','CD3E','TRAC','CD4','CD8A','TRDC','HBA1','CLEC9A', 'CLEC10A', 'FCGR3B','KLRB1','KLRD1'),max.cutoff = 'q90',min.cutoff = 'q10')

VlnPlot(integrated,features = c('CD79A','CD19','IGHM','IGHG','IGHD','CD14','FCGR3A','CD3E','TRAC','CD4','CD8A','TRDC','HBA1','CLEC9A', 'CLEC10A', 'FCGR3B','KLRB1','KLRD1'),stack=T)
DotPlot(integrated,features = c('CD79A','CD19','IGHM','CD14','FCGR3A','CD3E','TRAC','CD4','CD8A','TRDC','HBA1','CLEC9A', 'CLEC10A', 'FCGR3B','KLRB1','KLRD1'),assay = "RNA")

saveRDS(integrated,file = "integrated-new.rds")

####azimuth annotation

integrated <- Azimuth::RunAzimuth(integrated, reference = "pbmcref")
p1 <- DimPlot(integrated, group.by = "predicted.celltype.l1", label = TRUE, label.size = 3) + NoLegend()
p2 <- DimPlot(integrated, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
p3 <- DimPlot(integrated, group.by = "seurat_clusters", label = TRUE, label.size = 3) + NoLegend()
p1 + p2+ p3
##################################################################annotation
Tcell<-subset(integrated,idents = c(0,1,2,3,4,8,12,13,19))
Bcell<-subset(integrated,idents = c(9,17,18,20))
Myeloid<-subset(integrated,idents = c(5,6,7,10,15))
NK<-subset(integrated,idents = c(11))
unknown<-subset(integrated,idents = c(14,16))
##################################################################look at Myeloid
Myeloid$Disease<-substr(Myeloid$Patient,start=1,stop = nchar(Myeloid$Patient)-2)

DefaultAssay(Myeloid)<-"integrated"
Myeloid <- FindVariableFeatures(Myeloid)
Myeloid <- ScaleData(Myeloid)
Myeloid<- RunPCA(Myeloid,npcs=50)
ElbowPlot(Myeloid,50)
Myeloid <- FindNeighbors(Myeloid, dims = 1:30)
Myeloid <- RunUMAP(Myeloid, dims = 1:30)
Myeloid <- FindClusters(Myeloid,resolution = 0.5)

Myeloid <- Azimuth::RunAzimuth(Myeloid, reference = "pbmcref")
UMAPPlot(Myeloid,label=T,cols=brewer.pal(11,"Paired"))
UMAPPlot(Myeloid,split.by="Patient",ncol=3,cols=brewer.pal(11,"Paired"))+
guides(color = guide_legend(override.aes = list(size=4), ncol=11) )
table(Myeloid$Patient,Myeloid$seurat_clusters)
Idents(Myeloid)<-Myeloid$Disease
Myeloidmarker_BD<-FindMarkers(Myeloid,ident.1 = "BD")

#MyeloidMarkers<-FindAllMarkers(Myeloid,assay = "RNA")
MyeloidMarkers<-read.csv("Myeloidmarker.csv")
library(dplyr)
MyeloidMarkers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(Myeloid, features = top10$gene) + NoLegend()
DotPlot(Myeloid, features = top10$gene[!duplicated(top10$gene)])+coord_flip()
#common markers
myeloidcommon<-c('ITGAX', 'ITGAM', 'FCN1', 'CD14', 'FCGR3A', 'CD68', 'FCGR1A', "FCGR2B", "C1QC", "C1QA",'CCR7', 'CD274', 'XCR1', 'CLEC9A', 'CLEC10A', 'HLA-DRA', 'HLA-DRB1',"S100A8", "S100A9", "FCGR3B",'IRF8', "SIGLEC1", "KIT", "JCHAIN", "IL3RA",'CD3E',"KLRD1",'IGHM')
DotPlot(Myeloid, features = myeloidcommon,assay = "RNA")+coord_flip()



###########0508 myeloid t-test
m=read.csv("figure/myeloid.csv")
rownames(m)<-m$X
m=m[,-1]
aa=matrix(1,11)
for (i in c(1:11)){
  x<-m[1:3,i]#BD
  y<-m[4:6,i]#VKHD
  a=t.test(x, y, alternative = "two.sided")
  aa[i]=a$p.value
}


##################################################################look at Tcell
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
ShjVlnPlot(Tcell,features = c("TRAC","TRBC1", "TRBC2", "CD3D", "CD4",  "CD8A","CD8B"))
ShjVlnPlot(Tcell,features = c("CCR7", "SELL", "KLF2", "TCF7", "GZMK", "CCL5", "GZMA", "GZMB", "PRF1", "MKI67", "CD38", "IFNG", "TNF", "PDCD1", "CTLA4","FOXP3","IL2RA","CXCR6", "ITGB1","KLRB1"))
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
