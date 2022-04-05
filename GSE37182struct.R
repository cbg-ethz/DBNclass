#This script shows how to create different structural constraints for DBN models
#for GSE37182 dataset

source("DBNpreprocopt.R")

#this function can be used to construct a STRING-based penalization matrix (prior()
#de_breast_expanded - set of genes included in the analalysis
#stringintDBN = list of STRING interactions, two columens gene1 and gene2
#number of genes
pm_122<-edgepmDBNdirect(de_colon_expanded,stringintDBN,122)

#this function can be used to construct a STRING-based blacklist
#de_breast_expanded = set of genes included in the analalysis
#stringintDBN - list of STRING interactions, two columens gene1 and gene2
#number of genes
bl_122<-edgepmDBNdirect(de_colon_expanded,stringintDBN,122,intsame=0,intpf=0,basepf=1)

#prohibiting intra-edges, 0 means edge is allowed, 1 means edge is blacklisted
bl_122i<-matrix(1,nrow=2*122,ncol=2*122)
bl_122i[1:122,123:ncol(bl_122i)]<-0


