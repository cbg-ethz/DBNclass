#This is script for learning the model with the highest predictive accuracy
#for the dataset GSE37182
#The best model was defined by running parallelized cross-validation for all
#cases of structural and parametric assumptions ("CV_GSE37182_parallel.R")
#the results of CV are in the folder "biological_data_results"
#Structures are learned separately for cancer and normal biopsies
#Both MAP and consensus models are learned
#Intra-edges are not allowed (bl_122i.rds)


#search space for sampling identified using iterativeMCMC for both groups
commonspace<-readRDS("/Users/polinasuter/Downloads/DBNs/DBNclass/data/comcolon122.rds")
#data for the dataset data for the dataset GSE37182: differentially expressed genes and their transcription factors
datanorm<-readRDS("/Users/polinasuter/Downloads/DBNs/DBNclass/data/datanorm122.rds")
datacancer<-readRDS("/Users/polinasuter/Downloads/DBNs/DBNclass/data/datacancer122.rds")
#curbl - blacklist prohibiting intra-edges
curbl<-readRDS("/Users/polinasuter/Downloads/DBNs/DBNclass/data/bl_122i.rds")
#threshold for consensus model
p<-0.9

library(BiDAG)
pars_ns_normal<-scoreparameters("bge",datanorm,
                         dbnpar=list(samestruct=FALSE, slices=4,b=0,datalist=TRUE,stationary=FALSE),
                         DBN=TRUE,edgepmat=NULL)
sampn_ns_normal<-orderMCMC(pars_ns_normal, startspace=commonspace,chainout=TRUE,
                      verbose=TRUE,blacklist=curbl)
consn_ns_normal<-modelp(sampn_ns_normal,p=0.9,pdag=FALSE)

pars_ns_cancer<-scoreparameters("bge",datacancer,
                         dbnpar=list(samestruct=FALSE, slices=4,b=0,datalist=TRUE,stationary=FALSE),
                         DBN=TRUE,edgepmat=NULL)
sampn_ns_cancer <-orderMCMC(pars_ns_cancer, startspace=commonspace,chainout = TRUE,
                     verbose=TRUE,blacklist=curbl)
consn_ns_cancer<-modelp(sampn_ns_cancer,p=0.9,pdag=FALSE)

compareDBNs(consn_ns_normal,consn_ns_cancer)


