args = commandArgs(trailingOnly=TRUE)
nsamp<-as.numeric(args[1])
path<-"DBNsimres/"
DBNsimulation<-function(rep,n,nsamp=100) {
  library(pcalg)
  library(BiDAG)
  library(bnlearn)
  source("dbnsimfns.R")
  source("DBNCVparopt.R")

  path<-"DBNsimres/"

  simDBN<-generateDBN(n,rep+100)
  set.seed(rep+100)
  simdata<-rmvDAG(nsamp,simDBN,errDist="mix",mix=0.1)
  adj<-graph2m(simDBN)
  mean(apply(adj[,1:n+n],2,sum))
  cpadj<-adj
  cpadj[1:n,1:n]<-graph2m(dag2cpdag(m2graph(adj[1:n,1:n])))
  cpadj[1:n+n,1:n+n]<-graph2m(dag2cpdag(m2graph(adj[1:n,1:n])))


  hcfit5<-amat(hc(as.data.frame(simdata), blacklist = makeDBNblacklist(n),maxp=5))
  hcboot<-boot.strength(as.data.frame(simdata), algorithm = "hc", R=100, algorithm.args=list(blacklist = makeDBNblacklist(n),maxp=5))
  hcfitep<-makeAmatBoot(hcboot,n)
  hcfitcons5<-modellistDBN(hcfitep,p=c(0.3,0.5,0.7,0.9,0.99),n*2)


  hcfit3<-amat(hc(as.data.frame(simdata), blacklist = makeDBNblacklist(n),maxp=3))
  hcboot<-boot.strength(as.data.frame(simdata), algorithm = "hc", R=100, algorithm.args=list(blacklist = makeDBNblacklist(n),maxp=3))
  hcfitep<-makeAmatBoot(hcboot,n)
  hcfitcons3<-modellistDBN(hcfitep,p=c(0.3,0.5,0.7,0.9,0.99),n*2)

  hcres<-rbind(compareDBNs(hcfit3,adj,struct = "trans"),
  compareDBNs(hcfit5,adj,struct = "trans"),
  Reduce("rbind",lapply(hcfitcons3,compareDBNs,cpadj,"trans")),
  Reduce("rbind",lapply(hcfitcons5,compareDBNs,cpadj,"trans")))
  hcres<-as.data.frame(hcres)
  hcres$algorithm<-c("hc.3","hc.5",rep("hc.boot.3",5),rep("hc.boot.5",5))
  hcres$FDR<-hcres$FP/(hcres$FP+hcres$TP)
  hcres$threshold<-c(0,0,0.3,0.5,0.7,0.9,0.99,0.3,0.5,0.7,0.9,0.99)
  hcres$maxpar<-c(max(apply(hcfit3[,1:n+n],2,sum)),max(apply(hcfit5[,1:n+n],2,sum)),
                    unlist(lapply(hcfitcons3,function(x)max(apply(x,2,sum)))),
                  unlist(lapply(hcfitcons5,function(x)max(apply(x,2,sum)))))
  hcres$avpar<-c(mean(apply(hcfit3[,1:n+n],2,sum)),mean(apply(hcfit5[,1:n+n],2,sum)),
                 2*unlist(lapply(hcfitcons3,function(x)mean(apply(x,2,sum)))),
                 2*unlist(lapply(hcfitcons5,function(x)mean(apply(x,2,sum)))))

  dbnscore<-scoreparameters("bge",simdata,dbnpar=list(samestruct = FALSE, slices = 2, b = 0, stationary = TRUE, rowids = NULL,
                                                      datalist = NULL),DBN=TRUE)
  dbnfit<-iterativeMCMC(dbnscore,accum=FALSE,alpha=0.5,hardlim=14,cpdag=FALSE,startspace = hcfitcons5[[1]])
  dbnsamp<-orderMCMC(dbnscore,startspace = dbnfit$endspace,MAP=FALSE,chainout=TRUE)
  mcmcep<-edgep(dbnsamp)
  mcmccons<-modellistDBN(mcmcep,p=c(0.3,0.5,0.7,0.9,0.99),n*2)
  mcmcres<-rbind(compareDBNs(dbnfit$DAG,adj,struct="trans"),
                 Reduce("rbind",lapply(mcmccons,compareDBNs,adj,"trans")))
  mcmcres<-as.data.frame(mcmcres)
  mcmcres$algorithm<-c("mcmcMAP",rep("mcmcConsensus",5))
  mcmcres$FDR<-mcmcres$FP/(mcmcres$FP+mcmcres$TP)
  mcmcres$threshold<-c(0,0.3,0.5,0.7,0.9,0.99)
  mcmcres$maxpar<-c(max(apply(dbnfit$DAG[,1:n+n],2,sum)),
                    unlist(lapply(mcmccons,function(x)max(apply(x,2,sum)))))
  mcmcres$avpar<-c(mean(apply(dbnfit$DAG[,1:n+n],2,sum)),
                   2*unlist(lapply(mcmccons,function(x)mean(apply(x,2,sum)))))
  set.seed(rep+100)
  testdata<-rmvDAG(10,simDBN,errDist="mix",mix=0.1)
  mcmcres$MAE<-DBNCV2.sim(getDAG(dbnfit),mcmccons,testdata,dbnscore,n)
  hc3MAE<-DBNCV2.sim(hcfit3,hcfitcons3,testdata,dbnscore,n)
  hc5MAE<-DBNCV2.sim(hcfit5,hcfitcons5,testdata,dbnscore,n)

  hcres$MAE<-c(hc3MAE[1],hc5MAE[1],hc3MAE[2:length(hc3MAE)],
               hc5MAE[2:length(hc5MAE)])


res<-rbind(mcmcres,hcres)
rownames(res)<-c(1:nrow(res))
 saveRDS(res,paste(path,"dbnsim",nsamp,rep,".rds",sep=""))
 return(res)
}
library(parallel)
rep<-c(1:50)
cl <- makeCluster(52)
outputClApply <- parallel::clusterApply(cl, rep, DBNsimulation,
                                 n=100,nsamp=50)
stopCluster(cl)
saveRDS(outputClApply,paste(path,"DBNallsimwparN",nsamp,".rds",sep=""))



