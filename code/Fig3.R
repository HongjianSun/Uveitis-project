library(Seurat)
library(scRepertoire)
library(readxl)
library(ggplot2)
Tcell <- readRDS("C:/study/eye/2023.3/newest/Tcell.rds")
##highlight clonotype
Tcell$highlight<-as.numeric(substr(Tcell$highlight,10,11)) #let the number be sorted in the right order.
table(Tcell$highlight,Tcell$file)

Seurat::DimPlot(Tcell, split.by = "highlight",ncol = 4,label = T,label.size = 2.5)


table(Tcell$cloneType,Idents(Tcell),Tcell$Patient)
table(Tcell$cloneType,Idents(Tcell),Tcell$Disease)

#########
#pdf("C:/study/eye/2022.7/figure/Fig3/3a-UMAP clonotype level.pdf",15,3)
#UMAPPlot(Tcell,group.by="cloneType",label=T)
Idents(Tcell)<-Tcell$subclusterCD8
UMAPPlot(Tcell,split.by="cloneType",label=T)
#dev.off()


clonalOverlay(Tcell0616withcongaclusters, reduction = "umap", 
              freq.cutpoint = 30, bins = 5, facet = "file",contour_var = "ndensity") + 
  guides(color = FALSE)
###using proportion Sun
source("C:/study/eye/2022.4/Sun-clonalOverlay.R")
clonalOverlay(Tcell0616withcongaclusters, reduction = "umap", 
              prop.cutpoint = 0.01, bins = 5, facet = "file",contour_var = "ndensity") + 
  guides(color = FALSE)


####draw clonoindex plot
clonoindex<-read.csv("clonoindex.csv")
rownames(clonoindex) <- clonoindex$X
clonoindex$X <- NULL
cols2=as.matrix(cols)
rownames(cols2)=paste0("C",rownames(as.matrix(cols2)))

clonoindex2 <- data.frame(clonotype=factor(rep(c("rare","moderate","large"),ncol(clonoindex)*2),levels=c("rare","moderate","large")),disease=rep(c("VKHD","BD"),each =ncol(clonoindex)*3),
                          cluster=factor(rep(c(colnames(clonoindex),colnames(clonoindex)),each=nrow(clonoindex)/2),levels = rownames(as.matrix(cols2))),
                          index=unlist(clonoindex))

clonoindex2 <- data.frame(clonotype=factor(rep(c("rare","moderate","large"),ncol(clonoindex)*2),levels=c("rare","moderate","large")),disease=rep(c("VKHD","VKHD","VKHD","BD","BD","BD"),ncol(clonoindex)),
                          cluster=factor(rep(c(colnames(clonoindex)),each=nrow(clonoindex)),levels = rownames(as.matrix(cols2))),
                          index=unlist(clonoindex))
ggplot(clonoindex2, aes(x =clonotype , y = index, color = cluster, group = cluster)) + 
  geom_point() + geom_line() + 
  facet_grid(~ disease)+ scale_colour_manual(values = cols2)


clonalOverlay(CD4,reduction = "umap",freq.cutpoint = 0.005, bins = 10)+ scale_colour_manual(values = cols)|
  clonalOverlay(CD8,reduction = "umap",freq.cutpoint = 0.005, bins =10)+ scale_colour_manual(values = cols)