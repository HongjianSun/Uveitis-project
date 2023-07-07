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

#######################################################################BD2
#desired_length <- length(files) 
contig_list <- vector(mode = "list", length = 1)
contig_list[[1]]<-read.csv("E:/eye/source/tcr/new/BD-2.csv")
####strip barcode(use for scRNA-seq data)
combined1 <- combineTCR(contig_list,samples="BD-2",cells ="T-AB")
new <- combineExpression(combined1, Tcell, cloneCall="aa", proportion = T, cloneTypes=c(Small = 0.002, Medium = 0.01, Large = 0.05,SuperLarge=1))

new <- highlightClonotypes(new, 
                           cloneCall= "aa", 
                           sequence = rownames(head(sort(table(new@meta.data[["CTaa"]]),T),n=10)))
new$highlight<-factor(new$highlight,levels=rownames(head(sort(table(new@meta.data[["CTaa"]]),T),n=10)))
clonalOverlay(new,reduction = "umap",freq.cutpoint = 0.005, bins = 10)+ scale_colour_manual(values = cols)
DimPlot(new, split.by = "highlight",ncol = 3,label = T) + 
  theme(plot.title = element_blank())+ scale_colour_manual(values = cols)

#######################################################################BD3
#desired_length <- length(files) 
contig_list <- vector(mode = "list", length = 1)
contig_list[[1]]<-read.csv("E:/eye/source/tcr/new/BD-3.csv")
####strip barcode(use for scRNA-seq data)
combined1 <- combineTCR(contig_list,samples="BD-3",cells ="T-AB")
new <- combineExpression(combined1, Tcell, cloneCall="aa", proportion = T, cloneTypes=c(Small = 0.002, Medium = 0.01, Large = 0.05,SuperLarge=1))

new <- highlightClonotypes(new, 
                           cloneCall= "aa", 
                           sequence = rownames(head(sort(table(new@meta.data[["CTaa"]]),T),n=6)))
new$highlight<-factor(new$highlight,levels = rownames(head(sort(table(new@meta.data[["CTaa"]]),T),n=6)))
clonalOverlay(new,reduction = "umap",freq.cutpoint = 0.005, bins = 10)+ scale_colour_manual(values = cols)
DimPlot(new, split.by = "highlight",ncol = 3,label = T) + 
  theme(plot.title = element_blank())+ scale_colour_manual(values = cols)

#######################################################################VKH1
#desired_length <- length(files) 
contig_list <- vector(mode = "list", length = 1)
contig_list[[1]]<-read.csv("E:/eye/source/tcr/new/VKH-1.csv")
####strip barcode(use for scRNA-seq data)
combined1 <- combineTCR(contig_list,samples="VKH-1",cells ="T-AB")
new <- combineExpression(combined1, Tcell, cloneCall="aa", proportion = T, cloneTypes=c(Small = 0.002, Medium = 0.01, Large = 0.05,SuperLarge=1))

new <- highlightClonotypes(new, 
                           cloneCall= "aa", 
                           sequence = rownames(head(sort(table(new@meta.data[["CTaa"]]),T),n=15)))
new$highlight<-factor(new$highlight,levels=rownames(head(sort(table(new@meta.data[["CTaa"]]),T),n=15)))
clonalOverlay(new,reduction = "umap",freq.cutpoint = 0.005, bins = 10)+ scale_colour_manual(values = cols)
DimPlot(new, split.by = "highlight",ncol = 3,label = T) + 
  theme(plot.title = element_blank())+ scale_colour_manual(values = cols)

#######################################################################VKH2
#desired_length <- length(files) 
contig_list <- vector(mode = "list", length = 1)
contig_list[[1]]<-read.csv("E:/eye/source/tcr/new/VKH-2.csv")
####strip barcode(use for scRNA-seq data)
combined1 <- combineTCR(contig_list,samples="VKH-2",cells ="T-AB")
new <- combineExpression(combined1, Tcell, cloneCall="aa", proportion = T, cloneTypes=c(Small = 0.002, Medium = 0.01, Large = 0.05,SuperLarge=1))

new <- highlightClonotypes(new, 
                           cloneCall= "aa", 
                           sequence = rownames(head(sort(table(new@meta.data[["CTaa"]]),T),n=12)))
new$highlight<-factor(new$highlight,levels=rownames(head(sort(table(new@meta.data[["CTaa"]]),T),n=12)))
clonalOverlay(new,reduction = "umap",freq.cutpoint = 0.005, bins = 10)+ scale_colour_manual(values = cols)
DimPlot(new, split.by = "highlight",ncol = 3,label = T) + 
  theme(plot.title = element_blank())+ scale_colour_manual(values = cols)

#######################################################################VKH3
#desired_length <- length(files) 
contig_list <- vector(mode = "list", length = 1)
contig_list[[1]]<-read.csv("E:/eye/source/tcr/new/VKH-3.csv")
####strip barcode(use for scRNA-seq data)
combined1 <- combineTCR(contig_list,samples="VKH-3",cells ="T-AB")
new <- combineExpression(combined1, Tcell, cloneCall="aa", proportion = T, cloneTypes=c(Small = 0.002, Medium = 0.01, Large = 0.05,SuperLarge=1))

new <- highlightClonotypes(new, 
                           cloneCall= "aa", 
                           sequence = rownames(head(sort(table(new@meta.data[["CTaa"]]),T),n=27)))
new$highlight<-factor(new$highlight,levels=rownames(head(sort(table(new@meta.data[["CTaa"]]),T),n=27)))
clonalOverlay(new,reduction = "umap",freq.cutpoint = 0.005, bins = 10)+ scale_colour_manual(values = cols)
DimPlot(new, split.by = "highlight",ncol = 5,label = T) + 
  theme(plot.title = element_blank())+ scale_colour_manual(values = cols)


#########find HLA type
genes<-as.matrix(Tcell@assays[["RNA"]]@data@Dimnames[[1]])
FeaturePlot(Tcell,features = c('HLA-DRB1'),min.cutoff = 'q10',max.cutoff = 'q90')|UMAPPlot(Tcell,group.by="Disease")
DotPlot(Tcell,features = c('HLA-DRB1'),group.by="Disease")|VlnPlot(Tcell,features = c('HLA-DRB1'),group.by="Disease")



#########0304 compare clonoexpansion level GEX
C1_3expand<-subset(Tcell,idents ="1_3")#BD-1,BD-2
C1_3expandmarker1<-FindMarkers(C1_3expand,ident.1 = c("Largely (0.02 < X <= 1)"),ident.2 = "Rarely (0 < X <= 0.005)",group.by = "cloneType",assay = "RNA")
C1_3expandmarker2<-FindMarkers(C1_3expand,ident.1 = c("Largely (0.02 < X <= 1)"),ident.2 = "Morderately (0.005 < X <= 0.02)",group.by = "cloneType",assay = "RNA")
C1_3expandmarker3<-FindMarkers(C1_3expand,ident.1 = c("Largely (0.02 < X <= 1)"),group.by = "cloneType",assay = "RNA")

C1_0expand<-subset(Tcell,idents ="1_0")#BD-3
C1_0expandmarker1<-FindMarkers(C1_0expand,ident.1 = c("Largely (0.02 < X <= 1)"),ident.2 = "Rarely (0 < X <= 0.005)",group.by = "cloneType",assay = "RNA")
C1_0expandmarker2<-FindMarkers(C1_0expand,ident.1 = c("Largely (0.02 < X <= 1)"),ident.2 = "Morderately (0.005 < X <= 0.02)",group.by = "cloneType",assay = "RNA")
C1_0expandmarker3<-FindMarkers(C1_0expand,ident.1 = c("Largely (0.02 < X <= 1)"),group.by = "cloneType",assay = "RNA")

C1expand<-subset(Tcell,idents =1)
C1expandmarker1<-FindMarkers(C1expand,ident.1 = c("Largely (0.02 < X <= 1)"),ident.2 = "Rarely (0 < X <= 0.005)",group.by = "cloneType",assay = "RNA")
C1expandmarker2<-FindMarkers(C1expand,ident.1 = c("Largely (0.02 < X <= 1)"),ident.2 = "Morderately (0.005 < X <= 0.02)",group.by = "cloneType",assay = "RNA")
C1expandmarker3<-FindMarkers(C1expand,ident.1 = c("Largely (0.02 < X <= 1)"),group.by = "cloneType",assay = "RNA")


C4expand<-subset(Tcell,idents =4)
C4expandmarker1<-FindMarkers(C4expand,ident.1 = c("Largely (0.02 < X <= 1)"),ident.2 = "Rarely (0 < X <= 0.005)",group.by = "cloneType",assay = "RNA")
C4expandmarker2<-FindMarkers(C4expand,ident.1 = c("Largely (0.02 < X <= 1)"),ident.2 = "Morderately (0.005 < X <= 0.02)",group.by = "cloneType",assay = "RNA")
C4expandmarker3<-FindMarkers(C4expand,ident.1 = c("Largely (0.02 < X <= 1)"),group.by = "cloneType",assay = "RNA")


C0expand<-subset(Tcell,idents =0)
C0expandmarker1<-FindMarkers(C0expand,ident.1 = c("Largely (0.02 < X <= 1)"),ident.2 = "Rarely (0 < X <= 0.005)",group.by = "cloneType",assay = "RNA")
C0expandmarker2<-FindMarkers(C0expand,ident.1 = c("Largely (0.02 < X <= 1)"),ident.2 = "Morderately (0.005 < X <= 0.02)",group.by = "cloneType",assay = "RNA")
C0expandmarker3<-FindMarkers(C0expand,ident.1 = c("Largely (0.02 < X <= 1)"),group.by = "cloneType",assay = "RNA")
#C0expandmarker4<-FindMarkers(C0expand,ident.1 = c("Largely (0.02 < X <= 1)"),ident.2 = c("Rarely (0 < X <= 0.005)","Morderately (0.005 < X <= 0.02)"),group.by = "cloneType",assay = "RNA")

C2expand<-subset(Tcell,idents =2)
C2expandmarker1<-FindMarkers(C2expand,ident.1 = c("Largely (0.02 < X <= 1)"),ident.2 = "Rarely (0 < X <= 0.005)",group.by = "cloneType",assay = "RNA")
C2expandmarker2<-FindMarkers(C2expand,ident.1 = c("Largely (0.02 < X <= 1)"),ident.2 = "Morderately (0.005 < X <= 0.02)",group.by = "cloneType",assay = "RNA")
C2expandmarker3<-FindMarkers(C2expand,ident.1 = c("Largely (0.02 < X <= 1)"),group.by = "cloneType",assay = "RNA")


write.csv(C0expandmarker1,"markers/C0expandemarker large vs medium.csv")
write.csv(C0expandmarker2,"markers/C0expandemarker large vs small.csv")
write.csv(C1expandmarker1,"markers/C1expandemarker large vs medium.csv")
write.csv(C1expandmarker2,"markers/C1expandemarker large vs small.csv")
write.csv(C2expandmarker1,"markers/C2expandemarker large vs medium.csv")
write.csv(C2expandmarker2,"markers/C2expandemarker large vs small.csv")
write.csv(C4expandmarker1,"markers/C4expandemarker large vs medium.csv")
write.csv(C4expandmarker2,"markers/C4expandemarker large vs small.csv")

write.csv(C1_3expandmarker1,"markers/C1_3expandemarker large vs medium.csv")
write.csv(C1_3expandmarker2,"markers/C1_3expandemarker large vs small.csv")
write.csv(C1_0expandmarker1,"markers/C1_0expandemarker large vs medium.csv")
write.csv(C1_0expandmarker2,"markers/C1_0expandemarker large vs small.csv")

write.csv(C0expandmarker3,"markers/C0expandemarker large vs rest.csv")
write.csv(C1expandmarker3,"markers/C1expandemarker large vs rest.csv")
write.csv(C2expandmarker3,"markers/C2expandemarker large vs rest.csv")
write.csv(C4expandmarker3,"markers/C4expandemarker large vs rest.csv")
write.csv(C1_3expandmarker3,"markers/C1_3expandemarker large vs rest.csv")
write.csv(C1_0expandmarker3,"markers/C1_0expandemarker large vs rest.csv")



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
  geom_point() + geom_line() + 
  facet_grid(~ disease)+ scale_colour_manual(values = cols2)


clonalOverlay(CD4,reduction = "umap",freq.cutpoint = 0.005, bins = 10)+ scale_colour_manual(values = cols)|
  clonalOverlay(CD8,reduction = "umap",freq.cutpoint = 0.005, bins =10)+ scale_colour_manual(values = cols)