#this file shows major steps of gene selection for the GSE5462 dataset
#please note that the data can be found on the GEO server and is not allowed to be published on gthub
library(limma)
library(OmnipathR)

##get TF-targets interactions list from the Omnipath database
intTF <- import_dorothea_interactions(
  resources=c("DoRothEA"),
  dorothea_levels = c('A'),
  organism=org
)

#first construct the design matrix
colon_group<-colon_metadata$group
colon_group<-colon_metadata$time

cond<-factor(colon_group)
#set t=0 as a reference
cond <- relevel(cond, ref="T0")
design.trt=model.matrix(~cond)
colnames(design.trt)

#code for DGE analysys of previously normalized data
corfit <- duplicateCorrelation(colon_normalized, design.trt,ndups=1)
fit_colon <- lmFit(colon_rma, design.trt, cor = corfit$consensus.correlation)
efit_colon<-eBayes(fit_colon)

#get differentially expressed genes at t=1
de_colon_T1c <- topTable(efit_colon, coef = "condT1",n=10000)
length(which(de_colon_T1c[,"adj.P.Val"]<0.1)) #0

#get differentially expressed genes at t=2
de_colon_T2c <- topTable(efit_colon, coef = "condT2",n=10000)
length(which(de_colon_T2c[,"adj.P.Val"]<pval)) #1
head(de_colon_T2)
de_colon<-union(de_colon,de_colon_T2c$ID[which(de_colon_T2c$adj.P.Val<pval)])

#get differentially expressed genes at t=3
de_colon_T3c <- topTable(efit_colon, coef = "condT3",n=10000)
length(which(de_colon_T3c[,"adj.P.Val"]<pval)) #17
de_colon_T3c$ID[which(de_colon_T3c$adj.P.Val<pval)]

#get the union
de_colon<-union(de_colon,de_colon_T3c$ID[which(de_colon_T3c$adj.P.Val<pval)])

#extend the list with TF
respTF<-unique(intTF$source_genesymbol[which(intTF$target_genesymbol%in%de_colon)])

#final set of genes
GSE37182_genes<-unique(c(respTF,de_colon))
