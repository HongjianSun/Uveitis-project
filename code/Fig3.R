library(Seurat)
library(scRepertoire)
library(readxl)

Tcell <- readRDS("C:/study/eye/2022.4/Tcell0603.rds")
##highlight clonotype
Tcell$highlight<-as.numeric(substr(Tcell$highlight,10,11)) #let the number be sorted in the right order.
table(Tcell$highlight,Tcell$file)

Seurat::DimPlot(Tcell, split.by = "highlight",ncol = 4,label = T,label.size = 2.5)


"""
BD1 VKH1 VKH2
1    0  243    0
2    0    0  230
3    0   96    0
4    0    0   82
5    0    0   80
6    0    0   75
7    0    0   73
8    0    0   68
9    0    0   62
10  61    0    0
11   0   58    0
12   0    0   57
13  57    0    0
14   0   54    0
15   0    0   52
16   0   49    0
"""

UMAPPlot(Tcell,split.by="CloneType")

#########
pdf("C:/study/eye/2022.7/figure/Fig3/3a-UMAP clonotype level.pdf",15,3)
UMAPPlot(`Tcell-4groupclonetypes`,split.by="cloneType",label=T)
dev.off()


clonalOverlay(Tcell0616withcongaclusters, reduction = "umap", 
              freq.cutpoint = 30, bins = 5, facet = "file",contour_var = "ndensity") + 
  guides(color = FALSE)
###using proportion Sun
source("C:/study/eye/2022.4/Sun-clonalOverlay.R")
clonalOverlay(Tcell0616withcongaclusters, reduction = "umap", 
              prop.cutpoint = 0.01, bins = 5, facet = "file",contour_var = "ndensity") + 
  guides(color = FALSE)


