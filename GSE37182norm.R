#first download the data
data_colon <- getGEO("...colon_GSE37182_family.soft.gz")
colon_map<-makeAffyMapping(data_colon)
table_colon<-data2matrix(data_colon,mapping=colon_map)
table_colon_bc<-medianNormalization(table_colon, noLogTransform=TRUE)

