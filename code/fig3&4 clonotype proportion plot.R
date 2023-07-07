

####clonotype analysis

#clonotype proportion calculate
proportion=matrix()
for (i in c(1:length(Tcell$file))){
  
  if (Tcell$file[i]=="BD1"){
    proportion[i]=Tcell$Frequency[i]/table(Tcell$file)[1]
    if(is.na(Tcell$Frequency[i])){proportion[i]=1/table(Tcell$file)[1]}
  }
  else if (Tcell$file[i]=="VKH1"){
    proportion[i]=Tcell$Frequency[i]/table(Tcell$file)[2]
    if(is.na(Tcell$Frequency[i])){proportion[i]=1/table(Tcell$file)[2]}
  }
  else if (Tcell$file[i]=="VKH2"){
    proportion[i]=Tcell$Frequency[i]/table(Tcell$file)[3]
    if(is.na(Tcell$Frequency[i])){proportion[i]=1/table(Tcell$file)[3]}
  }
  else {
    print("ERROR!!!!")
    proportion[i]=NA
  }
}

Tcell0616withcongaclusters$cloneproportion<-proportion


##################0701 add conga hits to a clone family
###fix on 0708 (devide 0-10 and 0-3)
table(Tcell0616withcongaclusters$congacluster)
Tcell0616withcongaclusters$cloneproportionfix<-Tcell0616withcongaclusters$cloneproportion

for (i in c(1:11983)){
  if (!is.na(Tcell0616withcongaclusters$congacluster[i])){
    
    if (Tcell0616withcongaclusters$congacluster[i]=="0-10"){
      if (Tcell0616withcongaclusters$file[i]=="VKH1"){
        Tcell0616withcongaclusters$cloneproportionfix[i]=(21)/4557}
      else{Tcell0616withcongaclusters$cloneproportionfix[i]=(5)/6570}
    }
    
    if (Tcell0616withcongaclusters$congacluster[i]=="0-3"){
      if (Tcell0616withcongaclusters$file[i]=="VKH1"){
        Tcell0616withcongaclusters$cloneproportionfix[i]=(9)/4557}
      else{Tcell0616withcongaclusters$cloneproportionfix[i]=(5)/6570}
    }
    
    if (Tcell0616withcongaclusters$congacluster[i]=="4-2"){
      Tcell0616withcongaclusters$cloneproportionfix[i]=(18)/4557}  
    if (Tcell0616withcongaclusters$congacluster[i]=="5-7"){
      Tcell0616withcongaclusters$cloneproportionfix[i]=(15)/856} 
  }
}


######change clonotype class

cloneTypes=c(Small = 0.002, Medium = 0.01, Large = 0.05, 
             SuperLarge=1)
cloneTypes <- c(None = 0, cloneTypes)
for (x in seq_along(cloneTypes)) {
  names(cloneTypes)[x] <- paste0(names(cloneTypes[x]), 
                                 " (", cloneTypes[x - 1], " < X <= ", cloneTypes[x], 
                                 ")")
}
Tcell0616withcongaclusters$cloneType0701<-NA
for (i in 2:length(cloneTypes)) {
  Tcell0616withcongaclusters$cloneType0701 <- ifelse(Tcell0616withcongaclusters$cloneproportion > cloneTypes[i-1] & Tcell0616withcongaclusters$cloneproportion <= cloneTypes[i], names(cloneTypes[i]), 
                                                     Tcell0616withcongaclusters$cloneType0701)
}
table(Tcell0616withcongaclusters$cloneType0701,Tcell0616withcongaclusters$seurat_clusters,Tcell0616withcongaclusters$file)
table(Tcell0616withcongaclusters$cloneType0701,Tcell0616withcongaclusters$seurat_clusters)


################considering ientical+similar clonotype

cloneTypes=c(Small = 0.002, Medium = 0.01, Large = 0.05, 
             SuperLarge=1)
cloneTypes <- c(None = 0, cloneTypes)
for (x in seq_along(cloneTypes)) {
  names(cloneTypes)[x] <- paste0(names(cloneTypes[x]), 
                                 " (", cloneTypes[x - 1], " < X <= ", cloneTypes[x], 
                                 ")")
}
Tcell0616withcongaclusters$cloneType0701fix<-NA
for (i in 2:length(cloneTypes)) {
  Tcell0616withcongaclusters$cloneType0701fix <- ifelse(Tcell0616withcongaclusters$cloneproportionfix > cloneTypes[i-1] & Tcell0616withcongaclusters$cloneproportionfix <= cloneTypes[i], names(cloneTypes[i]), 
                                                        Tcell0616withcongaclusters$cloneType0701fix)
}
table(Tcell0616withcongaclusters$cloneType0701fix,Tcell0616withcongaclusters$seurat_clusters,Tcell0616withcongaclusters$file)
table(Tcell0616withcongaclusters$cloneType0701fix,Tcell0616withcongaclusters$seurat_clusters)

