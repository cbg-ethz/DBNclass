DBNCVcore<-function(dbndata,i,groupvector,bl=NULL,pm=NULL,supervised=TRUE,p=0.5,levels=c("responder","nonresponder"),
                    nruns=1,slices=2){

    K<-length(unique(levels))
    Di<-dbndata[i,]
    dbndata<-dbndata[-i,]
    groupvector<-groupvector[-i]
    res<-list()
    pr<-vector()
    for(j in 1:length(levels)) {
    pr[j]<-length(which(groupvector==levels[j]))/length(which(groupvector%in%levels))
    }
      res$prederror<-matrix(ncol=2,nrow=6)
      rownames(res$prederror)<-c("sMAP","scons", "MAP", "cons","sMAP_sep","scons_sep")
      datatoscore<-matrix(Di,nrow=1)

      set.seed(100)
      membi<-vector()
      res$lls<-vector()

      scoreClassMAP<-vector()
      scoreClassCons<-vector()
      modelMAP<-list()
      modelCons<-list()
      datalocal<-list()
      parlocal<-list()
      for(k in 1:K) {
        datalocal[[k]]<-dbndata[which(groupvector==levels[k]),]
        parlocal[[k]] <- scoreparameters("bge", datalocal[[k]], dbnpar = list(samestruct = TRUE, slices = 2, b = 0),
                                  DBN = TRUE, edgepmat = pm)
        fitresp <- iterativeMCMC(parlocal[[k]], verbose = FALSE, blacklist=bl,hardlimit=9,plus1it=7)
        sampresp <- orderMCMC(parlocal[[k]], verbose = FALSE, startspace=fitresp$endspace, MAP=FALSE,
                              chainout = TRUE, blacklist=bl)

        modelMAP[[k]]<-getDAG(fitresp)
        modelCons[[k]]<-modelp(sampresp,p=p,pdag=FALSE)
        scoreClassMAP[k]<-scoreagainstDBN(parlocal[[k]],modelMAP[[k]],datatoscore)
        scoreClassCons[k]<-scoreagainstDBN(parlocal[[k]],modelCons[[k]],datatoscore)

        fitresp<-NULL
        sampresp<-NULL
      }

      membi[1]<-assignClass(scoreClassMAP,prior=pr)
      membi[2]<-assignClass(scoreClassCons,prior=pr)
      res$prederror[1,]<- DBNCVpar2(modelMAP[[membi[1]]],DBNconsmat=NULL,parlocal[[membi[1]]],Di)
      res$prederror[2,]<-DBNCVpar2(modelMAP[[membi[2]]],DBNconsmat=modelCons[[membi[2]]],parlocal[[membi[2]]],Di)
      res$lls[1]<-scoreClassMAP[membi[1]]
      res$lls[2]<-scoreClassMAP[membi[2]]

      #DBNCVpar(modelMAP[[membi[1]]],DBNconsmat=NULL,datalocal[[membi[1]]],Di,slices = 2, b = 0)
      #DBNCVpar(modelMAP[[membi[2]]],DBNconsmat=modelCons[[membi[2]]],datalocal[[membi[2]]],Di,slices = 2, b = 0)


      for(k in 1:K) {
      parlocal[[k]] <- scoreparameters("bge", datalocal[[k]], dbnpar = list(samestruct = FALSE, slices = 2, b = 0),
                                DBN = TRUE, edgepmat = pm)
      fitresp <- iterativeMCMC(parlocal[[k]], verbose = FALSE, blacklist=bl,hardlimit=9,plus1it=7)
      sampresp <- orderMCMC(parlocal[[k]], verbose = FALSE, startspace=fitresp$endspace,MAP=FALSE,
                            chainout = TRUE, blacklist=bl)
      modelMAP[[k]]<-getDAG(fitresp)
      modelCons[[k]]<-modelp(sampresp,p=p,pdag=FALSE)
      scoreClassMAP[k]<-scoreagainstDBN(parlocal[[k]],modelMAP[[k]],datatoscore)
      scoreClassCons[k]<-scoreagainstDBN(parlocal[[k]],modelCons[[k]],datatoscore)
      fitresp<-NULL
      sampresp<-NULL
      }

      membi[5]<-assignClass(scoreClassMAP,prior=pr)
      membi[6]<-assignClass(scoreClassCons,prior=pr)

      res$prederror[5,]<-DBNCVpar2(modelMAP[[membi[5]]],DBNconsmat=NULL,parlocal[[membi[5]]],Di)
      res$lls[5]<-scoreClassMAP[membi[1]]
      res$prederror[6,]<-DBNCVpar2(modelMAP[[membi[6]]],DBNconsmat=modelCons[[membi[6]]],parlocal[[membi[6]]],Di)
      res$lls[6]<-scoreClassMAP[membi[6]]


      scoreall<-scoreparameters("bge", dbndata, dbnpar = list(samestruct = FALSE, slices = 2, b = 0),
                                DBN = TRUE, edgepmat = pm)
      fitall <- iterativeMCMC(scoreall, verbose = FALSE, blacklist=bl,hardlimit=9,plus1it=7)
      sampall <- orderMCMC(scoreall, verbose = FALSE, startspace=fitall$endspace,MAP=FALSE,
                           chainout = TRUE, blacklist=bl)
      modelCons<-modelp(sampall,p=p,pdag=FALSE)
      modelMAP<-getDAG(fitall)
      sampall<-NULL

      for(k in 1:K) {
        parlocal[[k]]<-scoreparameters("bge", datalocal[[k]], dbnpar = list(samestruct = FALSE, slices = 2, b = 0),
                                  DBN = TRUE, edgepmat = pm)
        scoreClassMAP[k]<-scoreagainstDBN(parlocal[[k]],modelMAP,datatoscore)
        scoreClassCons[k]<-scoreagainstDBN(parlocal[[k]],modelCons,datatoscore)
      }

      membi[3]<-assignClass(scoreClassMAP,prior=pr)
      membi[4]<-assignClass(scoreClassCons,prior=pr)
      res$prederror[3,]<-DBNCVpar2(modelMAP,DBNconsmat=NULL,parlocal[[membi[3]]],Di)
      res$lls[3]<-scoreClassMAP[membi[3]]
      res$prederror[4,]<-DBNCVpar2(modelMAP,DBNconsmat=modelCons,parlocal[[membi[4]]],Di)
      res$lls[4]<-scoreClassMAP[membi[4]]

      names(membi)<-c("sMAP","sConsensus","MAP","Consensus","sMAP_sep","sConsensus_sep")
      res$memb<-membi

    return(res)
}




