#This script shows how to create different structural constraints for DBN models
#_58 and _122 are models for the GSE37182 dataset


#STRING-based prior
pm_58<-edgepmDBNdirect(de_colon,stringintDBN,58)
pm_122<-edgepmDBNdirect(de_colon_expanded,stringintDBN,122)

#STRING-based blacklisting
bl_58<-edgepmDBNdirect(de_colon,stringintDBN,58,intsame=0,intpf=0,basepf=1)
bl_122<-edgepmDBNdirect(de_colon_expanded,stringintDBN,122,intsame=0,intpf=0,basepf=1)

#prohibiting intra-edges, 0 means edge is allowed, 1 means edge is blacklisted
bl_58i<-matrix(1,nrow=2*58,ncol=2*58)
bl_58i[1:58,59:116]<-0

bl_122i<-matrix(1,nrow=2*122,ncol=2*122)
bl_122i[1:122,123:ncol(bl_122i)]<-0


