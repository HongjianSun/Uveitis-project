library(Seurat)
library(dplyr)
setwd("C:/study/eye/2023.3")
#READ HLAresult
files<-list.files(path = "HLA")

desired_length <- length(files) 
contig_list <- vector(mode = "list", length = desired_length)
sample<-paste0("HLA/",files,"/alleles_exp.xls")

integrated <- readRDS("C:/study/eye/2023.3/integrated-new.rds")
integrated <- RenameCells(integrated,add.cell.id = integrated$Patient)
#Tcell <- RenameCells(integrated,add.cell.id = integrated$Patient)
integrated <- RenameCells(integrated,new.names =substr(colnames(integrated),start=1,stop = nchar(colnames(integrated))-2))

for (i in 1:length(contig_list)) {
  
  a<-read.table(sample[i],header = T)
  a <- data.frame(a[,-1], row.names = a[,1])
  colnames(a)<-gsub("\\.","-",colnames(a))
  colnames(a)<-paste0(files[i],"_",colnames(a))
  newrowname<-sapply(strsplit(rownames(a),":"), `[`, 1)
  a$HLA<-newrowname
  aa <- a %>%
    group_by(HLA) %>%
    summarise_all(sum)
  rownames(aa)<-aa$HLA
  #aa<-aa[,-1]
  #contig_list[[i]]<-as.data.frame(t(aa))
  contig_list[[i]]<-as.data.frame(aa)
  #integrated<-AddMetaData(integrated,as.data.frame(t(aa)))
}
merged_df <- Reduce(function(x, y) merge(x, y, by = "HLA", all = TRUE), contig_list)
HLAtype<-merged_df$HLA
merged_df<-merged_df[,-1]
rownames(merged_df)<-HLAtype
merged_df1<-merged_df
merged_df1[is.na(merged_df1)] <- 0
integrated<-AddMetaData(integrated,as.data.frame(t(merged_df1)))

table1<-as.data.frame(table(integrated$B.51,integrated$Patient))
table2<-as.data.frame(table(integrated$DRB1.04,integrated$Patient))
FeaturePlot(integrated,features = c("B.51","DRB1.04"),min.cutoff=-1,max.cutoff=10,split.by = "Patient",order=T)
FeaturePlot(integrated,features = c("B.51","DRB1.04"),split.by = "Patient")
FeaturePlot(integrated,features = c("B.51","DRB1.04"),min.cutoff=1,max.cutoff='q90',split.by = "Patient")
DotPlot(integrated,features = c("B.51","DRB1.04"),group.by="Patient")
VlnPlot(integrated,features = c("B.51","DRB1.04"),group.by="Patient")
DotPlot(Tcell,features = c("B.51","DRB1.04"),group.by = "subclusterCD8")
