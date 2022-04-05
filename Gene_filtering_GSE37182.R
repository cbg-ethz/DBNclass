#this file shows major steps of gene selection for the GSE5462 dataset
#please note that the data can be found on the GEO server and is not allowed to be published on gthub
library(limma)
library(OmnipathR)

#get TF-targets interactions list
intTF <- import_dorothea_interactions(
  resources=c("DoRothEA"),
  dorothea_levels = c('A'),
  organism=org
)


#code for DGE analysys of previously normalized data
#needs to be performed for resonders vs non-responders and
#for pre-treatment vs post-treatment
corfit <- duplicateCorrelation(colon_rma, design.trt,ndups=1)
fit_colon <- lmFit(colon_rma, design.trt, cor = corfit$consensus.correlation)
efit_colon<-eBayes(fit_colon)
de_colon <- topTable(efit_colon, coef = "condpost",n=22283)
de_colon<-de_colon[which(de_colon[,"adj.P.Val"]<0.05),]


#get gene names
#de_colon_symbol

colon_DE_genes_symbol<-unique(colon_DE_response_symbol,colon_DE_time_symbol)

#extend the list with TF
respTF<-unique(intTF$source_genesymbol[which(intTF$target_genesymbol%in%colon_DE_genes_symbol)])

#final set of genes
GSE37182_genes<-unique(c(respTF,colon_DE_genes_symbol))