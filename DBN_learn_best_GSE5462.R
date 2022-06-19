#This is script for learning the model with the highest predictive accuracy
#for the dataset GSE5462
#The best model was defined by running parallelized cross-validation for all
#cases of structural and parametric assumptions ("CV_GSE5462_parallel.R")
#the results of CV are in the folder "biological_data_results"
#Both MAP and consensus models are learned
#Non-STRING edges are blacklisted; "bl_125.rds"
#Intra-edges are allowed
#Initial and transition structures are learned separately, i.e. samestruct = FALSE


#data for the dataset GSE5462
#original data can be found here
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5462
#selected features (dbndata) include differentially expressed genes and their transcription factors
#dbndata<-readRDS("data/data_125.rds")

#STRING based blacklisting
curbl<-readRDS("data/bl_125.rds")
#threshold for consensus model
p<-0.9

library(BiDAG)
#score object, edgepmat = NULL denotes uniform prior over structures
scoreall<-scoreparameters("bge", dbndata, dbnpar = list(samestruct = FALSE, slices = 2, b = 0),
                          DBN = TRUE, edgepmat = NULL)
#find MAP dag
fitall <- iterativeMCMC(scoreall, verbose = TRUE, blacklist=curbl,hardlimit=9,plus1it=7)
#obtain a sample from the posterior distribution
sampall <- orderMCMC(scoreall, verbose = FALSE, startspace=fitall$endspace, MAP=FALSE,
                     chainout = TRUE, blacklist=curbl)
#estimate consensus model by performing model averaging
consall<-modelp(sampall,p=p,pdag=FALSE)

compareDBNs(consall,getDAG(fitall))



