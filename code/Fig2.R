setwd("C:/study/eye/2023.3")
library(Seurat)
library(ggplot2)
Tcell <- readRDS("C:/study/eye/2023.3/newest/Tcell.rds")

###################reset cluster color
"cols <-c(####old version
  '5'='#E0FFFF',
  '3'='#87CEFF',
  '6'='#63B8FF',
  '2'='#00BFFF',
  '10'='#1C86EE',
  '11'='#4F94CD',
  '0'='#2020EE',
  '4'='#00008B',
  
  '1_5'='#F7E4e8',
  '1_2'='#F0b080',
  '1_0'='#FF8347',
  '1_3'='#EE7621',
  '1_1'='#FF4500',
  '1_4'='#FF0000',
  '1_6'='#B22222',
  '12'='#B03060',
  
  '8'='#006400',
  '9'='#76EE00')
"
cols <-c(####new version 0426
  '5'='#E0FFFF',
  '3'='#87CEFF',
  '6'='#63B8FF',
  '2'='#00BFFF',
  '10'='#1C86EE',
  '11'='#4F94CD',
  '0'='#2020EE',
  '4'='#00008B',
  
  '1_5'='#ffffcc',
  '1_2'='#fed976',
  '1_0'='#fd8d3c',
  '1_3'='#fc4e2a',
  '1_1'='#e31a1c',
  '1_4'='#67000d',
  '1_6'='#a50f15',
  '12'='#cc4c02',
  
  '8'='#006400',
  '9'='#76EE00')
Idents(Tcell)<-factor(Idents(Tcell),levels = rownames(as.matrix(cols)))
UMAPPlot(Tcell)+ scale_colour_manual(values = cols)
UMAPPlot(Tcell,split.by="Patient",ncol=3)+ scale_colour_manual(values = cols)
source("C:/study/figure/ShjVlnplot-nolabel.R")
ShjVlnPlot(Tcell,features = c("TRAC","TRBC1", "TRBC2","TRDC","TRGC1","TRGC2", "CD3D", "CD4", "RUNX3","CD8A","CD8B","CD19"),cols=cols)+scale_x_reverse()

ShjVlnPlot(Tcell,features = c("CCR7", "SELL",  "TCF7","IFNG", "TNF","CCL4", "CCL5", "GZMA", "GZMB",  "GZMH", "GZMK","PRF1",  "CD38",  "PDCD1", "CTLA4","FOXP3","IL2RA","CXCR6",'ITGAE',"ITGB1","ZNF683","KLRB1",'CD19'),cols=cols)
DotPlot(Tcell,assay = "RNA",features = c('CD3E', 'CD4','CD8A','CD8B','SELL', 'CCR7', 'CD44', 'CD27','CD28','ITGA1', 'CD69', 'ITGAE', 'GZMK','TRAC', 'TRBC1', 'IFNG','FOXP3','BCL6','CXCR5','PDCD1','RORA','IL17RA','IL17F','RORC','KLRB1', 'TRDC', 'TRGC1','MKI67'))
DotPlot(Tcell,features = c("TRAC","TRBC1", "TRBC2","TRDC","TRGC1","TRGC2", "CD3D", "CD4", "RUNX3","CD8A","CD8B","CD19"),assay = "RNA")
DotPlot(Tcell,features = c("CCR7", "SELL",  "TCF7","IFNG", "TNF","CCL4", "CCL5", "GZMA", "GZMB",  "GZMH", "GZMK","PRF1",  "CD38",  "PDCD1", "CTLA4","FOXP3","IL2RA","CXCR6",'ITGAE',"ITGB1","ZNF683","KLRB1"),assay = "RNA")+theme( axis.text.x = element_text(angle = 45, hjust = 1.2))


####1A CD4/CD8 SUBPLOT

pdf("C:/study/eye/2022.7/figure/Fig1/1a-cd4.pdf",5,4)
FeaturePlot(Tcell,features = c("CD4"),max.cutoff = 'q90',min.cutoff = 'q10',order=T)
dev.off()

pdf("C:/study/eye/2022.7/figure/Fig1/1a-cd8a.pdf",5,4)
FeaturePlot(Tcell,features = c("CD8A"),max.cutoff = 'q90',min.cutoff = 'q10',pt.size = 0.1,order=T)
dev.off()

UMAPPlot(Tcell,label=T)|FeaturePlot(Tcell,features = c("CD4"),max.cutoff = 'q90',min.cutoff = 'q10')
FeaturePlot(Tcell,features = c("CD8A","CD8B"),max.cutoff = 'q90',min.cutoff = 'q10',blend=T)
UMAPPlot(Tcell,split.by="Patient",ncol=3)+scale_color_d3("category20")


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



##Fig2 heatmap
library(dplyr)
markers <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
DoHeatmap(Tcell, features = top10$gene) + NoLegend()
DotPlot(Tcell, features = top10$gene[!duplicated(top10$gene)])+coord_flip()#+ NoLegend()
###
Tcell$orig.ident<-Idents(Tcell)
Tcell <- RenameIdents(Tcell, 
                      '1_0' = 'CD8',
                      '1_1' = 'CD8',
                      '1_2' = 'CD8',
                      '1_3' = 'CD8',
                      '1_4' = 'CD8',
                      '1_5' = 'CD8',
                      '1_6' = 'CD8',
                      '12' = 'CD8',
                      '0' = 'CD4',
                      '2' = 'CD4',
                      '3' = 'CD4',
                      '4' = 'CD4',
                      '5' = 'CD4',
                      '6' = 'CD4',
                      '10' = 'CD4',
                      '11' = 'CD4',
                      '8' = 'gdT',
                      '9' = 'gdT')
Tcell$celltype<-Idents(Tcell)
Idents(Tcell)<-Tcell$orig.ident
DotPlot(Tcell,features = top10$gene[!duplicated(top10$gene)],split.by = "celltype",cols = c('#5e0fa2','#00441b','#DF3F00'))+coord_flip()


#Tcell=ScaleData(Tcell,assay="RNA")
pdf("figures/TcellHeatmap5.pdf",10,12)
DoHeatmap(Tcell,features =markers,assay="RNA", size = 5)
DoHeatmap(subset(Tcell, downsample = 1000),features =markers,assay="RNA", size = 5)
dev.off()

pdf("figures/TcellDotplot5.pdf",10,8)
DotPlot(Tcell,features =markers,assay="RNA")+coord_flip()      ## cols = c("blue", "yellow")
dev.off()
#check whether cluster 2 and 5 are MAIT
istrav7_2=grepl("TRAV7-2",Tcell0616withcongaclusters@meta.data[["CTgene"]])
table<-table(Tcell0616withcongaclusters$CTgene %in% "TRAV1-2")
istrav1_2[i]=TRUE

######################################
#Shj violin plot 2E
source("C:/study/figure/ShjVlnplot-nolabel.R")
source("C:/study/figure/ShjVlnplot.R")
ShjVlnPlot(Tcell,features = c("CCR7", "SELL", "KLF2", "TCF7", "GZMK", "CCL5", "GZMA", "GZMB", "PRF1", "MKI67", "CD38", "IFNG", "TNF", "PDCD1", "CTLA4","FOXP3","IL2RA","CXCR6", "ITGB1","KLRB1",'CD19'),cols = cols)


##改变显示时的顺序
#Tcell$integrated_merge_cluster <- factor(x =Tcell$seurat_clusters, levels = c('10','5','2','9','7','1','4','8','3','6','0'))
##C7到最后
Tcell$integrated_merge_cluster <- factor(x =Tcell$seurat_clusters, levels = c(2,8,15,0,1,3,4,6,7,11,14,5,10,13,9,12,16))


ShjVlnPlot(Tcell , features = c("CCR7", "SELL", "KLF2", "TCF7", "CCL5", "GZMA", "GZMB", "PRF1", "MKI67", "CD38", "IFNG", "TNF", "PDCD1", "CTLA4","FOXP3","IL2RA","CXCR6", "ITGB1","KLRB1","C1QA","C1QB"),assay = "RNA",group.by = 'integrated_merge_cluster')
ShjVlnPlot(Tcell , features = c("CCR7", "SELL", "KLF2", "TCF7", "CCL5", "GZMA", "GZMB", "PRF1", "MKI67", "CD38", "IFNG", "TNF", "PDCD1","CXCR6", "ITGB1","KLRB1","C1QA","C1QB", "CTLA4","FOXP3","IL2RA"),assay = "RNA",group.by = 'integrated_merge_cluster')

#S2A check whether cluster 2 and 5 are MAIT

istrav7_2=grepl("TRAV7-2",Tcell0616withcongaclusters@meta.data[["CTgene"]])
table<-table(Tcell0616withcongaclusters$CTgene %in% "TRAV1-2")
istrav1_2[i]=TRUE

ShjVlnPlot(Tcell,features = c("TRAC","TRBC1","TRBC2","TRDC","TRGC1","TRGC2","CD4","CD8A","CD8B"))#,group.by = 'integrated_merge_cluster'
ShjVlnPlot(Tcell,features = c("TRAV2","TRBV2","TRBV13","TRAV24"),cols=cols)#,group.by = 'integrated_merge_cluster'
DotPlot(Tcell,features = c("TRAV2","TRBV2","TRBV13","TRAV24"),assay = "RNA",scale.min = 0,scale.max = 60)+coord_flip()#,group.by = 'integrated_merge_cluster'
##2E stack vlnplot for common markers


#########S2a-c UMAP before removing C7 and C13
Tcell <- readRDS("C:/study/eye/2023.3/Tcell0322.rds")
#mki67 expression
pdf("C:/study/eye/2022.7/figure/Fig1/S2B-MKI67.pdf",5,4)
FeaturePlot(Tcell,features = c("MKI67"),max.cutoff = 'q90',min.cutoff = 'q10',pt.size = 0.1,order=T)
dev.off()


#0510 check cluster 1-6 and cluster 11 DGE correlated with type I interfore
markers16<-FindMarkers(CD8,ident.1 = "1_6")
markers11<-FindMarkers(CD4,ident.1 = "11")
write.csv(markers16,"marker1-6.csv")
write.csv(markers11,"markers11.csv")


DotPlot(Tcell,assay="RNA",features = c("CCR7", "SELL", "KLF2", "TCF7", "GZMK", "CCL5", "GZMA", "GZMB", "PRF1", "MKI67", "CD38", "IFNG", "TNF", "PDCD1", "CTLA4","FOXP3","IL2RA","CXCR6", "ITGB1","KLRB1"))+
  theme(axis.text.x = element_text(vjust=0.5,angle = 40))


####compare BD vs VKHD
#CD4
RNA_CD4VKHDmarker<-FindMarkers(CD4,ident.1 = "VKH",group.by = "Disease",assay = "RNA")
CD4VKHDmarker<-FindMarkers(CD4,ident.1 = "VKH",group.by = "Disease")
write.csv(RNA_CD4VKHDmarker,"newest/CD4VKHDmarker.csv")
#CD8
RNA_CD8VKHDmarker<-FindMarkers(CD8,ident.1 = "VKH",group.by = "Disease",assay = "RNA")
CD8VKHDmarker<-FindMarkers(CD8,ident.1 = "VKH",group.by = "Disease")
write.csv(RNA_CD8VKHDmarker,"newest/CD8VKHDmarker.csv")

#volcano CD4
library(ggplot2)
library(ggrepel)
RNA_CD4VKHDmarker<-FindMarkers(CD4,ident.1 = "VKH",group.by = "Disease",logfc.threshold = 0.15)
#write.csv(RNA_CD4VKHDmarker,"newest/CD4VKHDmarker-2.csv")
RNA_CD4VKHDmarker <- read.csv("C:/study/eye/2023.3/newest/CD4VKHDmarker-2.csv")
RNA_CD4VKHDmarker$threshold = factor(ifelse(RNA_CD4VKHDmarker$p_val_adj < 0.05 & abs(RNA_CD4VKHDmarker$avg_log2FC) >= 0.25, ifelse(RNA_CD4VKHDmarker$avg_log2FC>= 0.25 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
RNA_CD4VKHDmarker$Gene = RNA_CD4VKHDmarker$X
#RNA_CD4VKHDmarker$p_val_adj<-RNA_CD4VKHDmarker$p_val_adj+1e-310
ggplot(RNA_CD4VKHDmarker,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+ ggtitle("RNA_CD4VKHDmarker")+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
  geom_text_repel(
    data = RNA_CD4VKHDmarker[RNA_CD4VKHDmarker$p_val_adj<0.05&abs(RNA_CD4VKHDmarker$avg_log2FC)>0.25,],
    aes(label = Gene),
    size = 0,#6
    segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
  theme_bw()+#修改图片背景
  theme(
    legend.title = element_blank()#不显示图例标题
  )+
  ylab('-log10 (p-adj)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  geom_vline(xintercept=c(-0.5,0.5),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05

#volcano CD8
RNA_CD8VKHDmarker<-FindMarkers(CD8,ident.1 = "VKH",group.by = "Disease",logfc.threshold = 0.15)
#write.csv(RNA_CD8VKHDmarker,"newest/CD8VKHDmarker-2.csv")
RNA_CD8VKHDmarker <- read.csv("C:/study/eye/2023.3/newest/CD8VKHDmarker-2.csv")
RNA_CD8VKHDmarker$threshold = factor(ifelse(RNA_CD8VKHDmarker$p_val_adj < 0.05 & abs(RNA_CD8VKHDmarker$avg_log2FC) >= 0.25, ifelse(RNA_CD8VKHDmarker$avg_log2FC>= 0.25 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
RNA_CD8VKHDmarker$Gene = RNA_CD8VKHDmarker$X
ggplot(RNA_CD8VKHDmarker,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+ ggtitle("RNA_CD8VKHDmarker")+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
  geom_text_repel(
    data = RNA_CD8VKHDmarker[RNA_CD8VKHDmarker$p_val_adj<0.05&abs(RNA_CD8VKHDmarker$avg_log2FC)>0.25,],
    aes(label = Gene),
    size = 0,#6
    segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
  theme_bw()+#修改图片背景
  theme(
    legend.title = element_blank()#不显示图例标题
  )+
  ylab('-log10 (p-adj)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  geom_vline(xintercept=c(-0.5,0.5),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05

