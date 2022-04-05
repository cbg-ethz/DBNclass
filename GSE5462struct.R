#This script shows how to create different structural constraints for DBN models
#for GSE5462 dataset

source("DBNpreprocopt.R")

#this function can be used to construct a STRING-based penalization matrix (prior()
#de_breast_expanded - set of genes included in the analalysis
#stringintDBN = list of STRING interactions, two columens gene1 and gene2
#number of genes
pm_125<-edgepmDBNdirect(de_breast_expanded,stringintDBN,125)

#this function can be used to construct a STRING-based blacklist
#de_breast_expanded = set of genes included in the analalysis
#stringintDBN - list of STRING interactions, two columens gene1 and gene2
#number of genes
bl_125<-edgepmDBNdirect(de_breast_expanded,stringintDBN,125,intsame=0,intpf=0,basepf=1)

#prohibiting intra-edges, 0 means edge is allowed, 1 means edge is blacklisted
bl_125i<-matrix(1,nrow=2*125,ncol=2*125)
bl_125i[1:125,126:ncol(bl_125i)]<-0


