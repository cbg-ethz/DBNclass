#This is script for learning the model with the highest predictive accuracy
#for the dataset GSE37182
#The best model was defined by running parallelized cross-validation for all
#cases of structural and parametric assumptions ("CV_GSE37182_parallel.R")
#the results of CV are in the folder "biological_data_results"
#Structures are learned separately for cancer and normal biopsies
#Both MAP and consensus models are learned
#Intra-edges are not allowed (bl_122i.rds)


#search space for sampling identified using iterativeMCMC for both groups
commonspace<-readRDS("data/comcolon122.rds")
#data for the dataset data for the dataset GSE37182: differentially expressed genes and their transcription factors
datanorm<-readRDS("data/datanorm122.rds")
datacancer<-readRDS("data/datacancer122.rds")
#curbl - blacklist prohibiting intra-edges
curbl<-readRDS("data/bl_122i.rds")
#threshold for consensus model
p<-0.9

library(BiDAG)
pars_ns<-scoreparameters("bge",datanorm,
                         dbnpar=list(samestruct=FALSE, slices=4,b=0,datalist=TRUE,stationary=FALSE),
                         DBN=TRUE,edgepmat=NULL)
sampn_ns <- orderMCMC(pars_ns, startspace=commonspace,chainout = TRUE,
                      verbose=TRUE,blacklist=curbl)
consn_ns<-modelp(sampn_ns,p=0.9,pdag=FALSE)

