setwd("C:/study/eye/2022.7/data & code")
library(Seurat)
load("C:/study/eye/2022.7/data & code/Tcell0423.rdata")
####S1 all cells heatmap


metadata<-cbind(Tcell0616withcongaclusters@meta.data[1:25],Tcell0616withcongaclusters@meta.data[47:70])

Cluster25<-subset(Tcell0616withcongaclusters,idents = c(2,5))
metadata<-cbind(Cluster25@meta.data[1:25],Cluster25@meta.data[47:70])
Tcell202207<-Tcell
Tcell202207@meta.data<-metadata


Tcell0616withcongaclusters$istrav1_2=grepl("TRAV1-2",Tcell0616withcongaclusters@meta.data[["CTgene"]])
table(Tcell0616withcongaclusters$istrav1_2,Tcell0616withcongaclusters$seurat_clusters)

istrav1_2<-list()
for (i in c(1:11983)){
  if (Tcell0616withcongaclusters$CTgene[i] %in% "TRAV1-2"){print("bingo")}
    
  }
vizGenes(combined, gene = "V", 
         chain = "TRB", 
         plot = "bar", 
         order = "variance", 
         scale = TRUE)

Tcell0616withcongaclusters$istrav1_2 <-
table<-table(Tcell0616withcongaclusters$CTgene %in% "TRAV1-2")
istrav1_2[i]=TRUE