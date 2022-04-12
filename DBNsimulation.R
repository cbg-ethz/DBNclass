#THIS SCRIPT CAN BE USED TO REPLICATE FIGURE 3
#A LINE WITH A PATH TO SOVE THE RESULT MUST BE ADDED IN THE END

##########################################
#If simulations are run via command line #
##########################################
#first argument is the number replicates
#second argument is the number of nodes in one time slice of a DBN
args = commandArgs(trailingOnly=TRUE)
nsamp<-as.numeric(args[1])
n<-as.numeric(args[2])

#define path to save final result
path<-"DBNsimres/"

##########################################

#this function creates a DBN with n nodes in one time slice, 4 time slices and nsamp replicates
#rep - simulation replicate number
#n - number of nodes in one time slice
#nsmal - number of replicates (rowa) in the dataset
DBNsimulation<-function(rep,n,nsamp=100) {
  library(pcalg)
  library(BiDAG)
  library(bnlearn)
  library(dbnR)
  library(bnstruct)

  source("/Users/polinasuter/Downloads/DBNs/DBNclass/dbnsimfns.R")
  source("/Users/polinasuter/Downloads/DBNs/DBNclass/DBNCVparopt.R")

  #define path to save individual results
  path<-"DBNsimres/"


  set.seed(rep+100, kind = "Mersenne-Twister", sample.kind = 'Rounding', normal.kind = "Inversion")

  #simulate transition structure
  simDBN<-genDBN(1.3,n,2,lB=0.4,uB=1.1)

  #simulate data DBN and data
  simdata<-genDataDBN(simDBN,slices=4,ss=nsamp+2)
  rownames(simDBN$dbn)<-colnames(simDBN$dbn)<-c(paste("X",1:n,sep=""),paste("X",1:n,".2",sep=""))

  #define test and train data
  testdata<-simdata$mcmc[1:2+nsamp,]
  simdata$mcmc<-simdata$mcmc[1:nsamp,]
  simdata$gs<-simdata$gs[1:(nsamp*3),]
  colnames(simdata$gs)<-c(paste("X",1:n,sep=""),paste("X",1:n,".2",sep=""))

  adj<-simDBN$dbn
  mean(apply(adj[,1:n+n],2,sum))
  cpadj<-adj
  cpadj[1:n,1:n]<-graph2m(dag2cpdag(m2graph(adj[1:n,1:n])))
  cpadj[1:n+n,1:n+n]<-graph2m(dag2cpdag(m2graph(adj[1:n,1:n])))

  #MMHC from dbnR
  simdata$dbnr<-dbndata2dbnr(simdata$mcmc,4)
  f_dt_sim<-fold_dt(simdata$dbnr, 2)
  row_index<-1:nrow(f_dt_sim)
  delrows<-which(row_index %% 4 == 0)
  f_dt_sim<-f_dt_sim[-delrows,]
  starttime<-Sys.time()
  dbnrfit<-learn_dbn_struc(simdata$dbnr, 2, method="dmmhc",  f_dt=f_dt_sim)
  endtime<-Sys.time()
  runtime<-endtime-starttime
  dbnrfitmat<-amat(dbnrfit)
  dbnrfitmat[1:n,1:n+n]<-dbnrfitmat[1:n+n,1:n]
  dbnrfitmat[1:n+n,1:n]<-0
  dbnrres<-compareDBNs(dbnrfitmat,adj,struct = "trans")
  dbnrres$algorithm<-c("dbnR")
  dbnrres$FDR<-dbnrres$FP/(dbnrres$FP+dbnrres$TP)
  dbnrres$threshold<-0
  dbnrres$maxpar<-max(apply(amat(dbnrfit)[,1:n+n],2,sum))
  dbnrres$avpar<-mean(apply(amat(dbnrfit)[,1:n+n],2,sum))
  dbnrres$time<-as.numeric(runtime)
  dbnrres$units<-units(runtime)

  #hill climbing bnstruct
  dataset <- BNDataset(simdata$gs, rep(FALSE,n), node.sizes=rep(5,n),
                       variables=colnames(simdata$gs), num.time.steps = 2)
  starttime<-Sys.time()
  bnstructfit <- learn.dynamic.network(dataset, num.time.steps = 2, scoring.func="BIC", algo="hc")
  endtime<-Sys.time()
  runtime<-endtime-starttime

  bnstructres<-compareDBNs(bnstructfit@dag,adj,struct = "trans")
  bnstructres$algorithm<-c("bnstruct")
  bnstructres$FDR<-bnstructres$FP/(bnstructres$FP+bnstructres$TP)
  bnstructres$threshold<-0
  bnstructres$maxpar<-max(apply(bnstructfit@dag[,1:n+n],2,sum))
  bnstructres$avpar<-mean(apply(bnstructfit@dag[,1:n+n],2,sum))
  bnstructres$time<-as.numeric(runtime)
  bnstructres$units<-units(runtime)

  #hill climbing adaptation from the bnlearn package
  starttime<-Sys.time()
  hcfit5<-amat(hc(as.data.frame(simdata$gs), blacklist = makeDBNblacklist(n),maxp=5))
  hcboot<-boot.strength(as.data.frame(simdata$gs), algorithm = "hc", R=100, algorithm.args=list(blacklist = makeDBNblacklist(n),maxp=5))
  endtime<-Sys.time()
  runtime5<-endtime-starttime

  hcfitep<-makeAmatBoot(hcboot,n)
  hcfitcons5<-modellistDBN(hcfitep,p=c(0.3,0.5,0.7,0.9,0.99),2*n)


  starttime<-Sys.time()
  hcfit3<-amat(hc(as.data.frame(simdata$gs), blacklist = makeDBNblacklist(n),maxp=3))
  hcboot<-boot.strength(as.data.frame(simdata$gs), algorithm = "hc", R=100, algorithm.args=list(blacklist = makeDBNblacklist(n),maxp=3))
  endtime<-Sys.time()
  runtime3<-endtime-starttime

  hcfitep<-makeAmatBoot(hcboot,n)
  hcfitcons3<-modellistDBN(hcfitep,p=c(0.3,0.5,0.7,0.9,0.99),n*2)

  hcres<-rbind(compareDBNs(hcfit3,adj,struct = "trans"),
               compareDBNs(hcfit5,adj,struct = "trans"),
               Reduce("rbind",lapply(hcfitcons3,compareDBNs,adj,"trans")),
               Reduce("rbind",lapply(hcfitcons5,compareDBNs,adj,"trans")))
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
  hcres$time<-c(as.numeric(runtime3),as.numeric(runtime5),rep(as.numeric(runtime3),5),rep(as.numeric(runtime5),5))
  hcres$units<-c(units(runtime3),units(runtime5),
                 rep(units(runtime3),5),rep(units(runtime5),5))

  dbnscore<-scoreparameters("bge",simdata$mcmc,dbnpar=list(samestruct = TRUE, slices = 4, b = 0, stationary = TRUE, rowids = NULL,
                                                           datalist = NULL),DBN=TRUE)

  starttime<-Sys.time()
  dbnfit<-iterativeMCMC(dbnscore,accum=TRUE,alpha=0.2,plus1it=6,hardlim=10,cpdag=FALSE, alphainit=0.01, scoreout=TRUE)#,addspace=hcfitcons5[[1]])
  dbnsamp<-orderMCMC(dbnscore,scoretable = getSpace(dbnfit),MAP=FALSE,chainout=TRUE)
  endtime<-Sys.time()
  runtime<-endtime-starttime

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

  mcmcres$time<-as.numeric(runtime)
  mcmcres$units<-units(runtime)

  res<-rbind(mcmcres,hcres,dbnrres,bnstructres)
  rownames(res)<-c(1:nrow(res))

  #if(!is.null(path)) {
  #  saveRDS(res,paste(path,"MAR22dbnsimn",nsamp,n,rep,".rds",sep=""))
  #}
  return(res)
}


#this code can be used to run 50 replicates of simulations in parallel
library(parallel)
rep<-c(1:50)
cl <- makeCluster(51)
outputClApply <- parallel::clusterApply(cl, rep, DBNsimulation,
                                        n=n,nsamp=nsamp)
stopCluster(cl)

#ADD LINE TO SAVE THE RESULT!
saveRDS(outputClApply,paste(path,"MAR22all",nsamp,"n",n,".rds",sep=""))

