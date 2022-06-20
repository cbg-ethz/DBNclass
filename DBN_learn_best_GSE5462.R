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
dbndata<-readRDS("data/data_125.rds")

#STRING based blacklisting
curbl<-readRDS("data/bl_125.rds")
#threshold for consensus model
p<-0.9

library(BiDAG)
set.seed(100)
#score object, edgepmat = NULL denotes uniform prior over structures
scoreall<-scoreparameters("bge", dbndata, dbnpar = list(samestruct = FALSE, slices = 2, b = 0),
                          DBN = TRUE, edgepmat = NULL)
#find MAP dag
fitall <- iterativeMCMC(scoreall, verbose = TRUE, blacklist=curbl,hardlimit=9,plus1it=7)
#obtain a sample from the posterior distribution
sampall <- orderMCMC(scoreall, verbose = FALSE, startspace=fitall$endspace,chainout = TRUE, blacklist=curbl)
#estimate consensus model by performing model averaging
breast_net<-modelp(sampall,p=p,pdag=FALSE)

#for visualization
cellcycle_genes<-readRDS("visual/cell_cycle_genes.rds")
cellcycle_genes_plot<-c(cellcycle_genes,paste(cellcycle_genes,".",sep=""))
response_genes<-readRDS("visual/response_genes.rds")
response_genes_plot<-c(response_genes,paste(response_genes,".",sep=""))
data_125_map<-readRDS("visual/gene125symbols.rds")

colnames(breast_net)<-rownames(breast_net)<-data_125_map
genes_subset<-unique(c(cellcycle_genes,response_genes,
                       paste(cellcycle_genes,".2",sep=""),
                       paste(response_genes,".2",sep="")))
breast_subnet<-breast_net[genes_subset,genes_subset]

source("visual/plothelp.R")

plotGSE5462(breast_subnet,"trans",highlight1 = cellcycle_genes_plot,
           highlight2 = c(response_genes_plot),ts=28)



