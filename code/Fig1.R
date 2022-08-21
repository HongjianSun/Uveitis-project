setwd("C:/study/eye/2022.7/data & code")
library(Seurat)
load("C:/study/eye/2022.7/data & code/All cells.rdata")
VlnPlot(integrated,features = c("CD3E"),pt.size = 0)
Tcells<-subset(integrated,idents = c(0,1,2,3,4,6,7,8,9,10,11))
####S1 all cells heatmap
integrated<-ScaleData(integrated,assay = "RNA")
markers=cbind(markers0[1:5],markers1[1:5],markers2[1:5],markers3[1:5],markers4[1:5],markers5[1:5],markers6[1:5],markers7[1:5],markers8[1:5],markers9[1:5],markers10[1:5],markers11[1:5])
markers <- as.vector(markers)
pdf("C:/study/eye/2022.7/figure/Fig Supplement/FigS1-heatmap.pdf",10,8)
DoHeatmap(integrated,features =markers,assay="RNA", size = 5)
dev.off()

pdf("C:/study/eye/2022.7/figure/Fig Supplement/FigS1-dotplot.pdf",8,6)
DotPlot(integrated,features =c("CD3D","CD3E","CD3G","CD19","CD14"),assay = "RNA")
dev.off()

