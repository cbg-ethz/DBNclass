#this file shows major steps of gene selection for the GSE5462 dataset
#please note that the data can be found on the GEO server and is not allowed to be published on gthub
library(limma)
library(OmnipathR)

#get TF-targets interactions list from the Omnipath database
intTF <- import_dorothea_interactions(
  resources=c("DoRothEA"),
  dorothea_levels = c('A'),
  organism=org
)


#code for DGE analysys of previously normalized data
#needs to be performed for resonders vs non-responders and
#for pre-treatment vs post-treatment
corfit <- duplicateCorrelation(breast_rma, design.trt,ndups=1)
fit_breast <- lmFit(breast_rma, design.trt, cor = corfit$consensus.correlation)
efit_breast<-eBayes(fit_breast)
de_breast <- topTable(efit_breast, coef = "condpost",n=22283)
de_breast<-de_breast[which(de_breast[,"adj.P.Val"]<0.05),]

#determine breast_DE_time
#determine breast_DE_response

#get gene names
#breast_DE_response_symbol
#breast_DE_time_symbol

breast_DE_genes_symbol<-unique(breast_DE_response_symbol,breast_DE_time_symbol)

#extend the list with TF
respTF<-unique(intTF$source_genesymbol[which(intTF$target_genesymbol%in%breast_DE_genes_symbol)])

#final set of genes
GSE5462_genes<-unique(c(respTF,breast_DE_genes_symbol))
