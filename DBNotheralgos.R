#This scropt can be used to replicate the results from Table 1

#random forest
library(ranger)
library(naivebayes)

levs<-c("nonresponder","responder")
curdatarf<-data_125
data_rf<-as.data.frame(curdatarf[assessvec,])
data_rf$resp<-responsevector[assessvec]
data_rf$resp<-factor(data_rf$resp)
colnames(data_rf)<-make.names(colnames(data_rf))
rf_cv_glob<-vector()
nb_cv_glob<-vector()
for(j in 1:100) {
  rf_cv<-vector()
  nb_cv<-vector()
  for(i in 1:nrow(data_rf)) {
    Di<-i
    data_rf.train <- data_rf[-Di, ]
    data_rf.test <- data_rf[Di, ]
    rg.data_rf <- ranger(resp ~ ., data = data_rf.train)
    pred.data_rf <- predict(rg.data_rf, data = data_rf.test[-length(data_rf.test)])
    rf_cv[i]<-as.character(pred.data_rf$predictions)
    model <- naive_bayes(resp ~ ., data = data_rf.train)
    p <- predict(model, data_rf.test[-length(data_rf.test)], type = 'prob')
    nb_cv[i]<-colnames(p)[which.max(p)]
  }
  rf_cv_glob[j]<-1-length(which(rf_cv!=data_rf$resp))/52
  nb_cv_glob[j]<-1-length(which(nb_cv!=data_rf$resp))/52
}

CVclassRF<-mean(rf_cv_glob)
CVclassNB<-mean(nb_cv_glob)




