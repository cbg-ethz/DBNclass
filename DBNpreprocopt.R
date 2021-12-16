#library(stringr)
data2matrix<-function(x,mapping=NULL){
samples<-names(x@gsms)
datatable<-NULL
for(i in samples) {
  x_local<-x@gsms[[i]]@dataTable@table$VALUE
  datatable<-cbind(datatable,x_local)
}
colnames(datatable)<-samples
rownames(datatable)<-x@gsms[[i]]@dataTable@table$ID_REF
  if(!is.null(mapping)) {
    for(i in 1:nrow(datatable)) {
      curaffy<-rownames(datatable)[i]
      if(!is.na(mapping[curaffy,3])) rownames(datatable)[i]<-mapping[curaffy,3]
    }
  }
return(datatable)
}
makeAffyMapping<-function(x){
  newmap<-matrix(nrow=nrow(x@gsms[[1]]@dataTable@table), ncol=3)
  rownames(newmap)<-x@gsms[[1]]@dataTable@table[,1]
  cuttab<-x@gpls[[1]]@dataTable@table

  for(i in 1:nrow(newmap)) {
    curaffy<-rownames(newmap)[i]
    newmap[i,1]<-curaffy
    if(curaffy%in%cuttab[,1]) {
      ind<-which(cuttab[,1]==curaffy)
      colind<-which(grepl("Symbol",colnames(cuttab)))
      newmap[i,2]<-cuttab[ind,colind]
      if(newmap[i,2]=="") {
        newmap[i,2]<-NA
        newmap[i,3]<-NA
      } else {
        newmap[i,3]<-strsplit(newmap[i,2]," ///")[[1]][1]
      }
    } else {
      newmap[i,2]<-NA
      newmap[i,3]<-NA
    }
  }
  colnames(newmap)<-c("affy","symbolnu","symbol")
  return(newmap)
}
remap<-function(x,mapping) {
    for(i in 1:nrow(x)) {
      curaffy<-rownames(x)[i]
      if(!is.na(mapping[curaffy,3])) rownames(x)[i]<-mapping[curaffy,3]
    }
  return(x)
}
remapvector<-function(x,mapping) {
  rownames(mapping)<-mapping[,1]
  for(i in 1:length(x)) {
    curaffy<-x[i]
    if(!is.na(mapping[curaffy,3])) x[i]<-mapping[curaffy,3]
  }
  return(x)
}
procMeta<-function(x,dataset){
  if(dataset=="breast") {
    addcols<-matrix(ncol=3,nrow=nrow(x))
    colnames(addcols)<-c("group","response","dose")
    for(i in 1:nrow(x)) {
      if(grepl(" nonresponder", x[i,2])) {
        addcols[i,2]<-"nonresponder"
      } else if (grepl(" responder", x[i,2])) {
        addcols[i,2]<-"responder"
      } else {
        addcols[i,2]<-"unspecified"
      }
    }

    for(i in 1:nrow(x)) {
      if(grepl("pretreatment", x[i,2])) {
        addcols[i,1]<-"pretreatment"
        addcols[i,3]<-"pre"
      } else if (grepl("Letrozole", x[i,2])) {
        addcols[i,1]<-"posttreatment"
        if(grepl("2.5mg/day", x[i,2])) {
          addcols[i,3]<-"standard"
        } else {
          addcols[i,3]<-"non-standard"
        }
      } else {
        addcols[i,1]<-"unspecified"
        addcols[i,3]<-"unspecified"
      }
    }
  } else if (dataset=="cdk") {

    addcols<-matrix(ncol=5,nrow=nrow(x))
    colnames(addcols)<-c("group","type","time","replicate","dose")
    for(i in 1:nrow(x)) {
      a<-x[i,2]
      if(grepl("Untreated", a)) {
        addcols[i,1]<-"untreated"
        addcols[i,3]<-"0"
        addcols[i,5]<-"0"
      } else if (grepl("treated", a)) {
        addcols[i,1]<-"treated"
        t <- str_match(a, "for\\s*(.*?)\\s*hours")
        addcols[i,3]<-t[,2]
        if(grepl("DMSO",a) | grepl("vehicle",a)) {
          addcols[i,5]<-"DMSO"
        } else {
          addcols[i,5]<-getdose(a)
        }
      } else {
        addcols[i,1]<-"unspecified"
      }
      addcols[i,4] <- sub(".*replicate ", "", a)
    }


    for(i in 1:nrow(x)) {
      if(grepl("PBMCs", x[i,2])) {
        addcols[i,2]<-"PBMC"
      } else if (grepl("DU145", x[i,2])) {
        addcols[i,2]<-"DU145"
      } else if (grepl("HCT116", x[i,2])){
        addcols[i,2]<-"HCT116"
      } else {
        addcols[i,2]<-"unspecified"
      }
    }

  } else if(dataset=="colon") {
    addcols<-matrix(ncol=5,nrow=nrow(x))
    colnames(addcols)<-c("group","time","patient","replicate","ID")
    for(i in 1:nrow(x)) {
      a<-x[i,2]
      if(grepl("Normal", a)) {
        addcols[i,1]<-"normal"
      } else if (grepl("cancer", a)) {
        addcols[i,1]<-"cancer"
      }
      if(grepl("T0",a)) {
        addcols[i,2]<-"T0"
      } else if (grepl("T1",a)){
        addcols[i,2]<-"T1"
      } else if (grepl("T2",a)){
        addcols[i,2]<-"T2"
      } else if (grepl("T3",a)){
        addcols[i,2]<-"T3"
      }
      addcols[i,3] <- sub(" T.*", "", sub(".*patient ", "", a))
      addcols[i,4]<-sub(".* T. ", "", sub(" CL.*", "", a))
      addcols[i,5]<-sub(".*CL", "CL", a)
    }


  }
 x<-cbind(x,addcols)
 return(x)
}
getdose<-function(x){
  if(grepl("3xIC90",x)) {
    return("3xIC90")
  } else if(grepl("IC90",x)) {
    return("IC90")
  } else if(grepl("IC50",x)) {
    return("IC50")
  }
}
plotPCADBN<-function(dbndata,metadata,varpch,varcol,celltype="all") {

  rownames(metadata)<-metadata$sample_id
  if(celltype!="all") {
    loc_samps<-metadata$sample_id[which(metadata$type==celltype)]
    dbndata<-dbndata[,loc_samps]
    metadata<-metadata[loc_samps,]
  }

  var0<-which(apply(dbndata,1,sd)==0)
  if(length(var0)>0) dbndata<-dbndata[-var0,]
  pca_res <- prcomp(t(dbndata), scale. = TRUE)
  plot(pca_res$x[,1],pca_res$x[,2],col=as.integer(factor(metadata[,varcol])),
       pch=as.integer(factor(metadata[,varpch])),
         xlab="PC1",ylab="PC2", main="")
}
expr2dbn<-function(curdata) {
  nsamp<-nrow(curdata)
  curdata1<-curdata[0:(nsamp/2-1)*2+1,]
  curdata2<-curdata[0:(nsamp/2-1)*2+2,]
  colnames(curdata2)<-paste(colnames(curdata),".2",sep="")
  curdata<-cbind(curdata1,curdata2)
  return(curdata)
}
dbn2expr<-function(curdata) {
  nsamp<-nrow(curdata)
  newdata<-matrix(nrow=2*nsamp,ncol=ncol(curdata)/2)
  newdata[0:(nsamp-1)*2+1,]<-curdata[,1:(ncol(curdata)/2)]
  newdata[0:(nsamp-1)*2+2,]<-curdata[,1:(ncol(curdata)/2)+(ncol(curdata)/2)]
  colnames(newdata)<-colnames(curdata)[1:(ncol(curdata)/2)]
  return(newdata)
}
assignmembs<-function(vec1,vec2,prior=NULL) {
  if(is.null(prior)) {
    return(1*(vec2>vec1)+1)
  } else {
    tot<-cbind(vec1,vec2)
    nsamp<-nrow(tot)
    relativeprobabs <-tot # to store the relative probabilities
    for(s in 1:nsamp){ # run through the samples
      clusterscores<-tot[s,]
      maxscorey<-max(clusterscores) # take the maximum
      shifty<-exp(clusterscores-maxscorey)
      rescaley<-shifty/sum(shifty)
      relativeprobabs[s,]<-rescaley # relative probabilities
    }
    temp<-prior*t(relativeprobabs)
    relativeprobabswithtau<-1/colSums(temp)*t(temp) #
    return(1*(relativeprobabswithtau[,2]>relativeprobabswithtau[,1])+1)
  }
}
assignClass<-function(vec,prior=NULL) {
  if(!is.null(prior)) {
    maxscorey<-max(vec) # take the maximum
    shifty<-exp(vec-maxscorey)
    rescaley<-shifty/sum(shifty)
    vectau<-rescaley*prior
    vec<-vectau/sum(vectau)
  }
  return(which.max(vec))
}
DBNaccuracy<-function(estimatedvec,truevec,abs=FALSE) {
  if(abs) {
    truevec[which(truevec=="responder")]<-1
    truevec[which(truevec=="nonresponder")]<-2
    return(length(which(estimatedvec==truevec))/length(truevec))
  } else {
    return(adjustedRandIndex(estimatedvec,truevec))
  }
}
getCVmemb<-function(dbnCVres) {
  return(Reduce("rbind",lapply(dbnCVres,function(x)x$memb)))
}
meanCVerror<-function(dbnCVres) {
  return(Reduce("+",lapply(dbnCVres,function(x)x$prederror))/length(dbnCVres))
}
DBNdatabyslice<-function(DBNdata,metadata,samplecol=5,timecol=4,groupcol=3,stationary=TRUE) {
  rownames(metadata)<-metadata[,1]
  samples<-unique(metadata[,samplecol])
  times<-unique(metadata$time)
  datalist<-list()
  times_numeric<-as.numeric(sub("T","",times))
  groups<-unique(metadata[,groupcol])
  data_init<-NULL
  if(!stationary) initind<-max(times_numeric)+1 else initind<-1

  for(j in groups) {
    datalist[[j]]<-list()
    datalist[[j]][[initind]]<-data.frame()
    if(stationary) {
      datalist[[j]][[2]]<-data.frame()
    } else {
      for(k in 1:max(times_numeric)) {
        datalist[[j]][[k]]<-data.frame()
      }
    }
    for(i in samples) {
      sample_ids<-rownames(metadata[which(metadata[,samplecol]==i & metadata[,groupcol]==j),])
      if(length(sample_ids>0)) {
        Ts<-metadata[sample_ids,timecol]
        T0<-which(Ts=="T0")
        Ts_numeric<-as.numeric(sub("T","",Ts))
        time_order<-order(Ts_numeric)
        time_order<-time_order[-which(Ts=="T0")]

        init_loc<-data.frame(DBNdata[sample_ids[T0],])
        rownames(init_loc)<-sample_ids[T0]
        datalist[[j]][[initind]]<-rbind(datalist[[j]][[initind]],init_loc)
        data_init_mean<-apply(DBNdata[sample_ids[T0],],2,mean)

        data_local<-rbind(data_init_mean,DBNdata[sample_ids[time_order],])
        Ts<-c("T0",Ts[time_order])
        Ts_numeric<-c(0,Ts_numeric[time_order])
        if(stationary){
          for(t in 1:(length(Ts_numeric)-1)) {
            if(Ts_numeric[t+1]==(Ts_numeric[t]+1)) {
              newrow<-matrix(c(data_local[t,],data_local[t+1,]),nrow=1)
              rownames(newrow)<-rownames(data_local)[t+1]
              datalist[[j]][[2]]<-rbind(datalist[[j]][[2]],newrow)
            }
          }
        } else {
          for(t in 1:(length(Ts_numeric)-1)) {
            if(Ts_numeric[t+1]==(Ts_numeric[t]+1)) {
              slicet<-Ts_numeric[t]+1
              newrow<-matrix(c(data_local[t,],data_local[t+1,]),nrow=1)
              rownames(newrow)<-rownames(data_local)[t+1]
              datalist[[j]][[slicet]]<-rbind(datalist[[j]][[slicet]],newrow)
            }
          }
        }

      }
    }
  }
  for(j in groups) {
    for(i in 1:length(datalist[[j]])){
      if(ncol(datalist[[j]][[i]])==2*ncol(DBNdata)) {
        colnames(datalist[[j]][[i]])<-c(colnames(DBNdata),paste(colnames(DBNdata),".2",sep=""))
      }
    }
  }
  return(datalist)
}
DBNdatabysample<-function(DBNdata,metadata,samplecol=5,timecol=4,groupcol=3,stationary=TRUE) {
  rownames(metadata)<-metadata[,1]
  samples<-unique(metadata[,samplecol])
  times<-unique(metadata$time)
  datalist<-list()
  times_numeric<-as.numeric(sub("T","",times))
  times_numeric<-sort(times_numeric)
  groups<-unique(metadata[,groupcol])
  data_init<-NULL

  for(j in groups) {
    datalist[[j]]<-data.frame()
    for(i in samples) {
      sample_ids<-rownames(metadata[which(metadata[,samplecol]==i & metadata[,groupcol]==j),])
      data_local_all<-c()
      if(length(sample_ids>0)) {
        Ts<-metadata[sample_ids,timecol]
        T0<-which(Ts=="T0")
        Ts_numeric<-as.numeric(sub("T","",Ts))
        time_order<-order(Ts_numeric)
        time_order<-time_order[-which(Ts=="T0")]

        #datalist[[j]][[initind]]<-rbind(datalist[[j]][[initind]],data.frame(DBNdata[sample_ids[T0],]))
        data_local_all<-apply(DBNdata[sample_ids[T0],],2,mean)

        data_local<-DBNdata[sample_ids[time_order],]
        Ts<-Ts[time_order]
        Ts_numeric<-Ts_numeric[time_order]
        missingT<-setdiff(times_numeric[-1],Ts_numeric)
        nsmall<-ncol(DBNdata)
        if(length(missingT)>0) {
          if(missingT==1) {
            data_local<-rbind(rep(NA,nsmall),data_local)
          } else if (missingT==2) {
            data_local<-rbind(data_local,data_local[2,])
            data_local[2,]<-rep(NA,nsmall)
          } else if (missingT==3){
            data_local<-rbind(data_local,rep(NA,nsmall))
          }
        }
        for(t in times_numeric[-1]) {
          data_local_all<-c(data_local_all,data_local[t,])
        }
        datalist[[j]]<-rbind(datalist[[j]],data_local_all)
        rownames(datalist[[j]])[nrow(datalist[[j]])]<-i
      }

    }

  }

  basicnames<-colnames(DBNdata[,1:nsmall])
  allnames<-c(basicnames,paste(basicnames,".",rep(times_numeric[-1],each=nsmall),sep=""))
  colnames(datalist[[1]])<-allnames
  colnames(datalist[[2]])<-allnames
  return(datalist)
}
get1DBNsample<-function(dbndata,i) {
  if(i>nrow(dbndata$normal)) {
    return(dbndata$cancer[i-nrow(dbndata$normal),])
  } else {
    return(dbndata$normal[i,])
  }
}
remove1patient<-function(dbndata,sampleids) {
  for (i in 1:length(dbndata)) {
    delrows<-which(rownames(dbndata[[i]])%in%sampleids)
    if(length(delrows)>0) {
      dbndata[[i]]<-dbndata[[i]][-delrows,]
    }
  }
  return(dbndata)
}
DBNdatainit<-function(data,k){
  newdata<-list()
  samppergroup<-nrow(data)/k
  for(i in 1:k) {
    newdata[[i]]<-data[3*(1:samppergroup-1)+i,]
  }
  return(newdata)
}
edgepmDBN<-function(nodelabels,mapping,intlist, dyn,b=0,intsame=1,intpf=1,basepf=2){
  pm<-matrix(basepf,nrow=dyn,ncol=dyn)
  pmbig<-matrix(basepf,nrow=2*dyn,ncol=2*dyn)
  for(i in 1:dyn) {
    for(j in 1:dyn) {
      if(i==j) {
        pm[i,j]<-intsame
      } else {
        gene1<-mapping[nodelabels[i],3]
        gene2<-mapping[nodelabels[j],3]
        if(is.na(gene1) | is.na(gene2)) {
          pm[i,j]<-basepf
        } else {
          com<-paste(gene1,gene2,sep="")
          comr<-paste(gene2,gene1,sep="")
          if(com%in%rownames(intlist)) {
            pm[i,j]<-intpf
          } else if(comr%in%rownames(intlist)){
            pm[i,j]<-intpf
          } else {
            pm[i,j]<-basepf
          }
        }

      }
    }
  }

  pmbig[1:dyn,1:dyn]<-pmbig[1:dyn+dyn,1:dyn+dyn]<-pmbig[1:dyn,1:dyn+dyn]<-pm

  if(b>0){
  }
  return(pmbig)
}
edgepmDBNdirect<-function(nodelabels,intlist, dyn,b=0,intsame=1,intpf=1,basepf=2){
  pm<-matrix(basepf,nrow=dyn,ncol=dyn)
  pmbig<-matrix(basepf,nrow=2*dyn,ncol=2*dyn)
  for(i in 1:dyn) {
    for(j in 1:dyn) {
      if(i==j) {
        pm[i,j]<-intsame
      } else {
        gene1<-nodelabels[i]
        gene2<-nodelabels[j]
        if(is.na(gene1) | is.na(gene2)) {
          pm[i,j]<-basepf
        } else {
          com<-paste(gene1,gene2,sep="")
          comr<-paste(gene2,gene1,sep="")
          if(com%in%rownames(intlist)) {
            pm[i,j]<-intpf
          } else if(comr%in%rownames(intlist)){
            pm[i,j]<-intpf
          } else {
            pm[i,j]<-basepf
          }
        }

      }
    }
  }

  pmbig[1:dyn,1:dyn]<-pmbig[1:dyn+dyn,1:dyn+dyn]<-pmbig[1:dyn,1:dyn+dyn]<-pm

  if(b>0){
  }
  return(pmbig)
}
getInteractionList<-function(amat){
  intlist<-NULL
  labs<-sub("\\..*","",colnames(amat))
  for(i in 1:ncol(amat)) {
    for(j in 1:ncol(amat)) {
      if(amat[i,j]==1) {
        intlocal<-c(labs[i],labs[j])
        intlist<-rbind(intlist,intlocal)
      }
    }
  }
  colnames(intlist)<-c("gene1","gene2")
  rownames(intlist)<-paste(intlist[,1],intlist[,2],sep="")
  if(length(which(duplicated(rownames(intlist))))>0) {
    intlist<-intlist[-which(duplicated(rownames(intlist))),]
  }
  return(intlist)
}
compareInts<-function(int1,int2) {
  ff<-rep(F,nrow(int1))
  for(i in 1:nrow(int1)) {
    gene1<-int1[i,1]
    fg11<-which(int2[,1]==gene1)
    if(length(fg11)>0) {
      fg22<-which(int2[fg11,2]==int1[i,2])
      if(length(fg22)>0) ff[i]<-T
    }
    fg12<-which(int2[,2]==gene1)
    if(length(fg12)>0) {
      fg21<-which(int2[fg12,1]==int1[i,2])
      if(length(fg21)>0) ff[i]<-T
    }
  }
  return(ff)
}
sameInts<-function(int1) {
  return(length(which(int1[,1]==int1[,2])))
}
compareNeighbors<-function(model1,model2) {
  nsmall<-ncol(model1)/2
  res<-list()
  res$neibs<-list()
  res$abs<-vector()
  res$sim<-vector()
  res$rel<-vector()
  nodes<-colnames(model1)[1:nsmall+nsmall]
  nodesclean<-sub("\\..*","",nodes)
  k<-1
  for(i in nodes) {
    res$neibs[[nodesclean[k]]]<-list()
    res$neibs[[nodesclean[k]]][[1]]<-names(which(model1[,i]==1))
    res$neibs[[nodesclean[k]]][[2]]<-names(which(model2[,i]==1))
    n1<-length(res$neibs[[nodesclean[k]]][[1]])
    n2<-length(res$neibs[[nodesclean[k]]][[2]])
    diff12<-setdiff( res$neibs[[nodesclean[k]]][[1]], res$neibs[[nodesclean[k]]][[2]])
    diff21<-setdiff( res$neibs[[nodesclean[k]]][[2]], res$neibs[[nodesclean[k]]][[1]])
    inter<-intersect( res$neibs[[nodesclean[k]]][[2]], res$neibs[[nodesclean[k]]][[1]])
    res$abs[k]<-max(length(diff12),length(diff21))
    res$sim[k]<-length(inter)
    if(length(res$abs[k])==0) res$rel[k]=100 else res$rel[k]<-res$sim[k]/res$abs[k]
    k<-k+1
  }
  names(res$abs)<-names(res$sim)<-names(res$rel)<-nodesclean
  return(res)
}


