#Cross validation for the colon cancer dataset

DBNCV_MAE<-function(i,dbndata,metadata,bl,pm,pl1,levels=c("normal","cancer")){

  K<-length(unique(levels))
  data_classification<-DBNdatabysample(dbndata,metadata)
  nrows<-unlist(lapply(data_classification,nrow))

  totn<-sum(nrows)
  classvector<-rep(levels,nrows)
  trueclass<-classvector[i]
  abs_i<-i-min(which(classvector==trueclass))+1
  patid<-rownames(data_classification[[trueclass]])[i]
  removesamples<-rownames(metadata)[which(metadata$group==trueclass & metadata$patient==patid)]


  data_curs<-DBNdatabyslice(dbndata,metadata,stationary=TRUE)
  data_curns<-DBNdatabyslice(dbndata,metadata,stationary=FALSE)

    data_init<-DBNdatainit(data_curns[[trueclass]][[4]],3)
    for(j in 1:length(data_init)) {
      data_init[[j]]<-data_init[[j]][abs_i,]
    }
    data_curs[[trueclass]]<-remove1patient(data_curs[[trueclass]],removesamples)
    data_curns[[trueclass]]<-remove1patient(data_curns[[trueclass]],removesamples)


  Di<-get1DBNsample(data_classification,i)
  data_classification<-NULL

  data_local<-data_curs[[trueclass]]
  data_local_ns<-data_curns[[trueclass]]

  #fit a stationary model
  pars<-scoreparameters("bge",data_local,
                        dbnpar=list(samestruct=TRUE, slices=4,b=0,datalist=TRUE,stationary=TRUE),
                        DBN=TRUE,edgepmat=pm)
  itfitnormal_s<-iterativeMCMC(pars,plus1it=pl1,hardlim=10,blacklist = bl)
  sampn <- orderMCMC(pars, verbose = FALSE, startspace=itfitnormal_s$endspace,chainout = TRUE,blacklist = bl)
  consn<-modelp(sampn,p=0.9,pdag=FALSE)
  sampn<-NULL


  #second fit non-stationary
  pars_ns<-scoreparameters("bge",data_local_ns,
                           dbnpar=list(samestruct=FALSE, slices=4,b=0,datalist=TRUE,stationary=FALSE),
                           DBN=TRUE,edgepmat=pm)
  itfitnormal_ns<-iterativeMCMC(pars_ns,plus1it=pl1,hardlim=10,blacklist = bl)
  sampn_ns <- orderMCMC(pars_ns, verbose = FALSE, startspace=itfitnormal_ns$endspace,chainout = TRUE,blacklist = bl)
  consn_ns<-modelp(sampn_ns,p=0.9,pdag=FALSE)
  sampn_ns<-NULL


  #fit same DBN, stationary
  data_curs_all<-data_curs$normal
  data_curs_all[[1]]<-rbind(data_curs$normal[[1]],data_curs$cancer[[1]])
  data_curs_all[[2]]<-rbind(data_curs$normal[[2]],data_curs$cancer[[2]])

  parsall<-scoreparameters("bge",data_curs_all,
                           dbnpar=list(samestruct=TRUE, slices=4,b=0,datalist=TRUE,stationary=TRUE),
                           DBN=TRUE,edgepmat=pm)
  itfitall_s<-iterativeMCMC(parsall,plus1it=pl1,hardlim=10,blacklist = bl)
  sampnall <- orderMCMC(parsall, verbose = FALSE, startspace=itfitall_s$endspace,chainout = TRUE,blacklist = bl)
  consall<-modelp(sampnall,p=0.9,pdag=FALSE)
  sampnall<-NULL

  #fit same DBN non stationary
  data_curns_all<-data_curns$normal
  data_curns_all[[1]]<-rbind(data_curns$normal[[1]],data_curns$cancer[[1]])
  data_curns_all[[2]]<-rbind(data_curns$normal[[2]],data_curns$cancer[[2]])
  data_curns_all[[3]]<-rbind(data_curns$normal[[3]],data_curns$cancer[[3]])
  data_curns_all[[4]]<-rbind(data_curns$normal[[4]],data_curns$cancer[[4]])

  parsall_ns<-scoreparameters("bge",data_curns_all,
                           dbnpar=list(samestruct=FALSE, slices=4,b=0,datalist=TRUE,stationary=FALSE),
                           DBN=TRUE,edgepmat=pm)
  itfitall_ns<-iterativeMCMC(parsall_ns,plus1it=pl1,hardlim=10,blacklist = bl)
  sampnall <- orderMCMC(parsall_ns, verbose = FALSE, startspace=itfitall_ns$endspace,chainout = TRUE,blacklist = bl)
  consall_ns<-modelp(sampnall,p=0.9,pdag=FALSE)
  sampnall<-NULL


  res<-list()
  res$pie<-matrix(nrow=8,ncol=3)
  rownames(res$pie)<-c("MAP","sMAP","sMAP_ns","MAP_ns","Cons","sCons","sCons_ns","Cons_ns")


  res$pie[1,]<-DBNCVparmult(itfitall_s$DAG,NULL,pars,Di)
  res$pie[2,]<-DBNCVparmult(itfitnormal_s$DAG,NULL,pars,Di)
  res$pie[3,]<-DBNCVparmult(itfitnormal_ns$DAG,NULL,pars_ns,Di)
  res$pie[4,]<-DBNCVparmult(itfitall_ns$DAG,NULL,pars_ns,Di)

  res$pie[5,]<-DBNCVparmult(itfitall_s$DAG,consall,pars,Di)
  res$pie[6,]<-DBNCVparmult(itfitnormal_s$DAG,consn,pars,Di)
  res$pie[7,]<-DBNCVparmult(itfitnormal_ns$DAG,consn_ns,pars_ns,Di)
  res$pie[8,]<-DBNCVparmult(itfitall_ns$DAG,consall_ns,pars_ns,Di)

return(res)

}
