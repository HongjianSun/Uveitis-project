#0312 SEE KEY clonotype seperately
setwd("C:/study/eye/2023.3")
library(Seurat)
library(scRepertoire)

Tcell <- readRDS("C:/study/eye/source/Tcell.rds")

########################################integrate TCR to Tcell object0412
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

table(Tcell$cloneType,Tcell$subclusterCD8,Tcell$Patient)
a=table(Tcell$cloneType,Idents(Tcell),Tcell$Disease)
a1=as.data.frame(a[,,1])#BD
a2=as.data.frame(a[,,2])#VKH
#draw figure(finished0417)
a1 <- within(a1,{ Var2 <- factor(Var2,levels = rownames(as.matrix(cols)))})
a2 <- within(a2,{ Var2 <- factor(Var2,levels = rownames(as.matrix(cols)))})
ggplot(a1, aes(fill=Var2, y=Freq, x=Var1)) +   geom_bar(position="fill", stat="identity")+ scale_fill_manual(values = cols)+ theme_classic()|
  ggplot(a2, aes(fill=Var2, y=Freq, x=Var1)) +   geom_bar(position="fill", stat="identity")+ scale_fill_manual(values = cols)+ theme_classic()

###reverse 230512
table(Tcell$cloneType,Tcell$subclusterCD8,Tcell$Patient)
a=table(Idents(Tcell),Tcell$cloneType,Tcell$Disease)
a1=as.data.frame(a[,,1])#BD
a2=as.data.frame(a[,,2])#VKH
#draw figure(finished0417)
a1 <- within(a1,{ Var2 <- factor(Var2,levels = rownames(as.matrix(cols)))})
a2 <- within(a2,{ Var2 <- factor(Var2,levels = rownames(as.matrix(cols)))})
ggplot(a1, aes(fill=Var2, y=Freq, x=Var1)) +   geom_bar(position="fill", stat="identity")+ scale_fill_manual(values = brewer.pal(3,"Dark2"))+ theme_classic()|
  ggplot(a2, aes(fill=Var2, y=Freq, x=Var1)) +   geom_bar(position="fill", stat="identity")+ scale_fill_manual(values = brewer.pal(3,"Dark2"))+ theme_classic()



#Tcell <- RenameCells(Tcell,add.cell.id = Tcell$Patient)
#Tcell <- RenameCells(integrated,add.cell.id = integrated$Patient)
#Tcell <- RenameCells(Tcell,new.names =substr(colnames(Tcell),start=1,stop = nchar(colnames(Tcell))-2))
Tcell <- combineExpression(combined, Tcell, cloneCall="aa", proportion = T, cloneTypes=c(Rarely = 0.005, Morderately = 0.02, Largely = 1))
UMAPPlot(Tcell,split.by="cloneType")+ scale_colour_manual(values = cols)

################set new function of clonaloverlay

trace(clonalOverlay, edit=TRUE)#    , contour_var = "ndensity",lwd = 0.5,alpha=0.2
Idents(Tcell)<-Tcell$subclusterCD8
UMAPPlot(Tcell,label=T)
UMAPPlot(CD4,label=T)|UMAPPlot(CD8,label=T)
clonalOverlay(Tcell,reduction = "umap",freq.cutpoint = 0.005, bins =5,facet='Patient')+ scale_colour_manual(values = cols)
clonalOverlay(Tcell,reduction = "umap",freq.cutpoint = 0.005, bins = 5,facet='Disease')+ scale_colour_manual(values = cols)

clonalOverlay(CD4,reduction = "umap",freq.cutpoint = 0.005, bins = 5)+ scale_colour_manual(values = cols)|
clonalOverlay(CD8,reduction = "umap",freq.cutpoint = 0.005, bins =5)+ scale_colour_manual(values = cols)

#################other check of TCR
#4 Visualizing Contigs
quantContig(combined, cloneCall="aa", scale = TRUE)
quantContig_output <- quantContig(combined, cloneCall="aa", 
                                  scale = TRUE, exportTable = TRUE)
quantContig_output

abundanceContig(combined, cloneCall = "aa", scale = FALSE)
abundanceContig(combined, group = "sample", scale = TRUE)
lengthContig(combined, cloneCall="aa", chain = "both",scale=T) 
compareClonotypes(combined, numbers = 20, samples = c("BD-1", "BD-2","BD-3"), 
                  cloneCall="aa", graph = "alluvial")

head(sort(table(Tcell$Patient,Tcell$CTaa),T),n=15)

topclone<-c("CIVRVISPRETSGSRLTF_CASSQRQGAAYEQYF","NA_CASSQRQGAAYEQYF",#BD-1
            "CAVIRNYGQNFVF_CASSFGRDEQFF","CAVIRNYGQNFVF;CAVEERHSGNTPLVF_CASSFGRDEQFF",#BD-2
            "CADRTGNQFYF_CASSYRDRLNTEAFF","CAVNEGGSEKLVF_CASSGGMNYGYTF",#BD-3
            "CALRSNAGNMLTF_CASSPGQGNYGYTF","CAVTRRGDNARLMF_CASSLSGTGELFF",#VKH-1
            "CAASKGNAGNMLTF_CASSGAGTGNSPLHF","CAASTKAAGNKLTF_CASSLRQGGTGELFF",#VKH-2
            "CAASQAGTALIF_CSARDYREGELFF","CAVAHKEYGNKLVF_CASSQGIGGTEAFF")#VKH-3

trace(highlightClonotypes, edit=TRUE)#    , contour_var = "ndensity
Tcell <- highlightClonotypes(Tcell,cloneCall= "aa",sequence = topclone)

DimPlot(Tcell, split.by = "highlight",ncol = 6) +theme(plot.title = element_blank())


###############################################################################################################SEE SEPERATELY

a=as.matrix(table(Tcell$CTaa,Tcell$Patient))
#a[,1]
aa=matrix(nrow = 30, ncol = 6,dimnames = list(c(1:30),c("BD1","BD2","BD3","VKH1","VKH2","VKH3")))
aa[,1]=head(as.matrix(a[order(a[,1],decreasing = T),1]),n=30)/sum(a[,1])#bd1
aa[,2]=head(as.matrix(a[order(a[,2],decreasing = T),2]),n=30)/sum(a[,2])#bd2
aa[,3]=head(as.matrix(a[order(a[,3],decreasing = T),3]),n=30)/sum(a[,3])#bd3
aa[,4]=head(as.matrix(a[order(a[,4],decreasing = T),4]),n=30)/sum(a[,4])#vkh1
aa[,5]=head(as.matrix(a[order(a[,5],decreasing = T),5]),n=30)/sum(a[,5])#vkh2
aa[,6]=head(as.matrix(a[order(a[,6],decreasing = T),6]),n=30)/sum(a[,6])#vkh3
aa
  write.csv(a,"C:/study/eye/2023.3/figure/expanded clones.csv")
  write.csv(aa,"C:/study/eye/2023.3/figure/top clonetypes-0519.csv")
# more than 0.5% expanded clonotypes:   "BD1":5,"BD2":10,"BD3":6,"VKH1":15,"VKH2":12,"VKH3:27"###0323
#######################################################################BD1
#desired_length <- length(files) 
contig_list <- vector(mode = "list", length = 1)
contig_list[[1]]<-read.csv("E:/eye/source/tcr/new/BD-1.csv")
####strip barcode(use for scRNA-seq data)
combined1 <- combineTCR(contig_list,samples="BD-1",cells ="T-AB")
new <- combineExpression(combined1, Tcell, cloneCall="aa", proportion = T, cloneTypes=c(Small = 0.002, Medium = 0.01, Large = 0.05,SuperLarge=1))

new <- highlightClonotypes(new, 
                              cloneCall= "aa", 
                              sequence = rownames(head(sort(table(new@meta.data[["CTaa"]]),T),n=5)))
new$highlight<-factor(new$highlight,levels=rownames(head(sort(table(new@meta.data[["CTaa"]]),T),n=5)))
clonalOverlay(new,reduction = "umap",freq.cutpoint = 0.005, bins = 10)+ scale_colour_manual(values = cols)
DimPlot(new, split.by = "highlight",ncol = 3,label = T) + 
  theme(plot.title = element_blank())+ scale_colour_manual(values = cols)

#########find HLA type
genes<-as.matrix(Tcell@assays[["RNA"]]@data@Dimnames[[1]])
FeaturePlot(Tcell,features = c('HLA-DRB1'),min.cutoff = 'q10',max.cutoff = 'q90')|UMAPPlot(Tcell,group.by="Disease")
DotPlot(Tcell,features = c('HLA-DRB1'),group.by="Disease")|VlnPlot(Tcell,features = c('HLA-DRB1'),group.by="Disease")



#compare in CD4 or CD8
CD4expandmarkers<-FindMarkers(CD4,ident.1 = c("Largely (0.02 < X <= 1)"),group.by = "cloneType",assay = "RNA")
CD8expandmarkers<-FindMarkers(CD8,ident.1 = c("Largely (0.02 < X <= 1)"),group.by = "cloneType",assay = "RNA")
write.csv(CD4expandmarkers,"newest/CD4expandmarkers.csv")
write.csv(CD8expandmarkers,"newest/CD8expandmarkers.csv")


####draw clonoindex plot(0512fix)
clonoindex<-read.csv("clonoproportion.csv")
rownames(clonoindex) <- clonoindex$X
clonoindex$X <- NULL
cols2=as.matrix(cols)
rownames(cols2)=paste0("C",rownames(as.matrix(cols2)))

clonoindex2 <- data.frame(clonotype=factor(rep(c("rare","moderate","large"),ncol(clonoindex)*2),levels=c("rare","moderate","large")),disease=rep(c("VKHD","VKHD","VKHD","BD","BD","BD"),ncol(clonoindex)),
                          cluster=factor(rep(c(colnames(clonoindex)),each=nrow(clonoindex)),levels = rownames(as.matrix(cols2))),
                          index=unlist(clonoindex))
ggplot(clonoindex2, aes(x =clonotype , y = index, color = cluster, group = cluster)) + 
  geom_point() + geom_line(linewidth = 1) + scale_y_log10(breaks = c(0.1,0.2,0.5, 1,2,5,10,20,50,100))+
  facet_grid(~ disease)+ scale_colour_manual(values = cols2)+theme_bw()#+theme(panel.grid = element_blank())


clonalOverlay(CD4,reduction = "umap",freq.cutpoint = 0.005, bins = 10)+ scale_colour_manual(values = cols)|
  clonalOverlay(CD8,reduction = "umap",freq.cutpoint = 0.005, bins =10)+ scale_colour_manual(values = cols)

write.csv(Tcell@meta.data,"Tmetadata.csv")


#volcano CD4
library(ggplot2)
library(ggrepel)
CD4expandmarkers <- read.csv("C:/study/eye/2023.3/newest/CD4expandmarkers.csv")
CD4expandmarkers$threshold = factor(ifelse(CD4expandmarkers$p_val_adj < 0.05 & abs(CD4expandmarkers$avg_log2FC) >= 0.5, ifelse(CD4expandmarkers$avg_log2FC>= 0.5 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
CD4expandmarkers$Gene = CD4expandmarkers$X
ggplot(CD4expandmarkers,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+ ggtitle("CD4expandmarkers")+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
  geom_text_repel(
    data = CD4expandmarkers[CD4expandmarkers$p_val_adj<0.05&abs(CD4expandmarkers$avg_log2FC)>0.5,],
    aes(label = Gene),
    size = 6,
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
CD8expandmarkers <- read.csv("C:/study/eye/2023.3/newest/CD8expandmarkers.csv")
CD8expandmarkers$threshold = factor(ifelse(CD8expandmarkers$p_val_adj < 0.05 & abs(CD8expandmarkers$avg_log2FC) >= 0.5, ifelse(CD8expandmarkers$avg_log2FC>= 0.5 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
CD8expandmarkers$Gene = CD8expandmarkers$X
ggplot(CD8expandmarkers,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+ ggtitle("CD8expandmarkers")+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
  geom_text_repel(
    data = CD8expandmarkers[CD8expandmarkers$p_val_adj<0.05&abs(CD8expandmarkers$avg_log2FC)>0.5,],
    aes(label = Gene),
    size = 5,
    segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
  theme_bw()+#修改图片背景
  theme(
    legend.title = element_blank()#不显示图例标题
  )+
  ylab('-log10 (p-adj)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  geom_vline(xintercept=c(-0.5,0.5),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05

