#code for parallelized LOOCV for GSE5462 dataset

#command line arguments
#arg1 data
#arg2 responsevector
#arg3 name
#arg4 blacklist file
#arg5 penalization matrix file

#arguments specifying the evaluated model
#for example: args<-c("data_19","responsevector","breast","no","pm_19")
makeid<-paste("sup",args[1],sep="")
curdata<-readRDS(paste("data/",args[3],"/",args[1],".rds",sep=""))#
responsevector<-readRDS(paste("data/",args[3],"/",args[2],".rds",sep=""))#
if (args[4]!="no") bl<-readRDS(paste("data/",args[3],"/",args[4],".rds",sep="")) else bl<-NULL
if (args[5]!="no") pm<-readRDS(paste("data/",args[3],"/",args[5],".rds",sep="")) else pm<-NULL

DBNCV<-function(ind,dbndata, groupvector, bl=NULL, pm=NULL,
                dataset="breast",id="supervised", nruns=1){
  library(BiDAG)
  library(clue)
  library(pcalg)
  library(mclust)
  source("DBNCVopt.R")
  source("DBNCVparopt.R")
  source("DBNpreprocopt.R")

  res<-DBNCVcore(dbndata,i=ind,groupvector=groupvector, bl=bl, pm=pm, supervised=TRUE, p=0.5,
                 levels=c("responder","nonresponder"),nruns=1,slices=2)
  file<-paste(id,ind,".rds",sep="")
  path<-paste(dataset,"/",file,sep="") #path for saving the result of each CV run
  saveRDS(res,path)
  return(res)
}

library(parallel)
cl <- makeCluster(54)
rep<-c(1:52)
outputClApply <- parallel::clusterApply(cl, rep, DBNCV,
                                        dbndata=curdata,
                                        groupvector=responsevector,
                                        dataset="breast",id="sup",
                                        bl=bl,pm=pm,
                                        nruns=1)
stopCluster(cl)
