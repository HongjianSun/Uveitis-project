setwd("C:/study/eye/2023.3")
library(Seurat)
integrated <- readRDS("C:/study/eye/2023.3/newest/integrated.rds")
VlnPlot(integrated,features = c("CD3E"),pt.size = 0)
Tcells<-subset(integrated,idents = c(0,1,2,3,4,6,7,8,9,10,11))

integrated <- readRDS("integrated-new.rds")

#FeaturePlot(integrated,features = c("CD3E","CD19","CD14"),max.cutoff = 'q90',min.cutoff = 'q10',ncol=3)
FeaturePlot(integrated,features = c("CD3E","CD19","CD14"),max.cutoff = 'q90',min.cutoff = 'q10',ncol=3,order=T)
####S1 all cells heatmap
integrated<-ScaleData(integrated,assay = "RNA")

###top5 heatmap/dotplot

library(dplyr)
markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) -> top10
#DoHeatmap(integrated, features = top10$gene) + NoLegend()
DotPlot(integrated, features = top10$gene[!duplicated(top10$gene)],assay = "RNA")+coord_flip()
#RNA(all genes)
markersRNA <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1,assay = "RNA")
markersRNA %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) -> top10RNA
#DoHeatmap(integrated, features = top10RNA$gene) + NoLegend()
DotPlot(integrated, features = top10RNA$gene[!duplicated(top10RNA$gene)],assay = "RNA")+coord_flip()
My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  axis.text.y = element_text(size = 10))
DotPlot(integrated, features = top10$gene[!duplicated(top10$gene)])+coord_flip()+My_Theme
DotPlot(integrated, features = top10RNA$gene[!duplicated(top10RNA$gene)])+coord_flip()+My_Theme


#integrated=ScaleData(integrated,assay="RNA")
pdf("figure/AllcellsHeatmap5.pdf",10,12)
DoHeatmap(integrated,features =markers,assay="RNA", size = 5)
dev.off()


saveRDS(integrated,file="integrated-scale.rds")
markers=cbind(markers0[1:5],markers1[1:5],markers2[1:5],markers3[1:5],markers4[1:5],markers5[1:5],markers6[1:5],markers7[1:5],markers8[1:5],markers9[1:5],markers10[1:5],markers11[1:5])
markers <- as.vector(markers)
pdf("C:/study/eye/2022.7/figure/Fig Supplement/FigS1-heatmap.pdf",10,8)
DoHeatmap(integrated,features =markers,assay="RNA", size = 5)
dev.off()
###Dotplot
pdf("C:/study/eye/2022.7/figure/Fig Supplement/FigS1-heatmap.pdf",10,8)
#DoHeatmap(integrated,features =markers,assay="RNA", size = 5)

pdf("figure/AllcellsHeatmap5.pdf",10,12)
DoHeatmap(subset(integrated, downsample = 1000),assay="RNA", features = markers)
dev.off()



pdf("C:/study/eye/2023.3/figure/Fig Supplement/FigS1-dotplot.pdf",8,4)
DotPlot(integrated,features =c('CD3D','CD3E','CD3G','CD19','CD79A','IGHM','CD14','FGCR3A','CLEC9A', 'CLEC10A', 'FCGR3B','KLRB1','KLRD1'),assay = "RNA")
dev.off()


###new color mode
 
cols <-c('9'='#FFEC8B',
'17'='#FA8072',#
'18'='#fe9929',
'20'='#d95f0e',

'11'='#e7d4e8',#NK
'16'='#756bb1',


'0'='#9ecae1',
'1'='#6baed6',
'2'='#3182bd',
'3'='#08519c',
'12'='#c6dbef',
'13'='#eff3ff',#CD4

'4'='#fcae91',#CD8
'8'='#fb6a4a',
'19'='#a50f15',

'5'='#7fcdbb',
'6'='#005a32',
'7'='#698B22',
'10'='#d9f0a3',#Myeloid
'15'='#78c679',
'14'='#d9d9d9')
Idents(integrated)<-factor(Idents(integrated),levels = rownames(as.matrix(cols)))
UMAPPlot(integrated)+ scale_colour_manual(values = cols)
UMAPPlot(integrated,split.by="Patient",ncol=2)+ scale_colour_manual(values = cols)
UMAPPlot(integrated,group.by="Patient")
UMAPPlot(integrated,group.by="Disease")
#1e dotplot
DotPlot(integrated,assay="RNA",features = c('CD79A','CD19','IGHM','IGHG1','IGHD','CD14','FCGR3A','CD3E','TRAC','CD4','CD8A','TRDC','HBA1','CLEC9A', 'CLEC10A', 'FCGR3B','KLRB1','KLRD1','NCAM1'))

DotPlot(integrated,assay="RNA",cols = c("#FFFFFF","blue"),col.min = -0.5,
        features = c('CD79A','CD19','IGHM','IGHG1','IGHD','CD14','FCGR3A','CD3E','TRAC','CD4','CD8A','TRDC','HBA1','CLEC9A', 'CLEC10A', 'FCGR3B','KLRB1','KLRD1','NCAM1'),
        )
DotPlot(integrated,assay="RNA",cols = c("#FFFFFF","blue"),
        features = c('CD14','FCGR3A','CLEC9A', 'CLEC10A', 'FCGR3B','CD3E','TRAC','CD4','CD8A','TRDC','KLRB1','KLRD1','NCAM1','CD19','IGHM','IGHG1','IGHD','CD79A'),
)+theme(axis.text.x = element_text(angle = 20))


DotPlot(integrated,assay="RNA",
        features = c('CD14','FCGR3A','CLEC9A', 'CLEC10A', 'FCGR3B','CD3E','TRAC','CD4','CD8A','RORC','RORA','GNLY','NCAM1','NCR1','CD19','IGHM','IGHG1','IGHD','CD79A'),
)+theme(axis.text.x = element_text(angle = 20))

DotPlot(integrated,assay="RNA",
        features = c('CD14','FCGR3A','CLEC9A', 'CLEC10A', 'FCGR3B','CD3E','TRAC','TRBC1','TRBC2','CD4','CD8A','TRDC','TRGC1','TRGC2','RORC','KLRB1','KLRD1','CD19','IGHM','IGHG1','IGHD','CD79A','KIT', 'IL3RA','ITGA2','NCAM1','ZBTB16','CD44','IL2RB','CD69'),
)+theme(axis.text.x = element_text(angle = 20))
