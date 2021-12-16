#code for parallelized LOOCV for GSE37182 dataset

#command line arguments
#arg1 data
#arg2 metadata
#arg3 name
#arg4 blacklist file
#arg5 penalization matrix file

#should be rub for all sets of parameters, below are optimal
#args<-c("data_122","colon_metadata","colon","bl_122i","no")
args = commandArgs(trailingOnly=TRUE)
makeid<-paste("MAE",args[1],args[4],args[5],sep="")
curdata<-readRDS(paste(path,args[1],".rds",sep=""))#
metadata<-readRDS(paste(path,args[2],".rds",sep=""))#
if (args[4]!="no") bl<-readRDS(paste(path,args[4],".rds",sep="")) else bl<-NULL
if (args[5]!="no") pm<-readRDS(paste(path,args[5],".rds",sep="")) else pm<-NULL

DBNCV_MAEcore<-function(ind,dbndata, metadata,  bl=NULL, pm=NULL,
                dataset="colon",id="MAE"){
  library(bnclustOmics)
  library(BiDAG)
  library(clue)
  library(pcalg)
  library(mclust)
  source("DBNCVcolonMAEopt.R")
  source("DBNCVparopt.R")
  source("DBNpreprocopt.R")

  res<-DBNCV_MAE(i=ind,dbndata,metadata=metadata,
                  bl=bl, pm=pm, pl1=4)
  file<-paste(ind,".rds",sep="")
  saveRDS(res,file)
  return(res)
}

library(parallel)
cl <- makeCluster(31)
rep<-c(1:29)
outputClApply <- parallel::clusterApply(cl, rep, DBNCV_MAEcore,
                                        dbndata=curdata,
                                        metadata=metadata,
                                        dataset="colon",id=makeid,
                                        bl=bl,pm=pm)
stopCluster(cl)

