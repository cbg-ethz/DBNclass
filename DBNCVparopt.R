#The function estimates MAP parameters of DBN model

#@param incidence adjacency matrix of a DBN
#@param scorepar  object of class score parameters (package BiDAG)
MAPparametersDBN<-function(incidence,scorepar) {
  res<-list()
  res$mus<-list()
  res$betas<-list()
  scorepar$split<-TRUE
  incidence<-BiDAG:::DBNbacktransform_l(incidence,scorepar)

  if(scorepar$stationary) {
    res$betas$init<-matrix(0,nrow=scorepar$firstslice$n,ncol=scorepar$firstslice$n)
    res$mus$init<-scorepar$firstslice$muN
    mainnodes<-scorepar$firstslice$mainnodes
    for (j in mainnodes)  {
      parentnodes <- which(incidence$init[,j]==1)
      if(length(parentnodes)>0) {
        res$betas$init[parentnodes,j]<-MAPparametersDBNcore(j,parentnodes,scorepar$firstslice)
      }
    }
    res$betas$trans<-matrix(0,nrow=scorepar$otherslices$n,ncol=scorepar$otherslices$n)
    res$mus$trans<-vector()
    MAPparametersDBNmu(j,parentnodes, scorepar$otherslices)

    mainnodes<-scorepar$otherslices$mainnodes
    for (j in mainnodes)  {
      parentnodes <- which(incidence$trans[,j]==1)
      res$mus$trans[j]<-MAPparametersDBNmu(j,parentnodes, scorepar$otherslices)
      if(length(parentnodes)>0) {
        res$betas$trans[parentnodes,j]<-MAPparametersDBNcore(j,parentnodes,scorepar$otherslices)
      }
    }
    res$betas$trans<-BiDAG:::DBNtransform(res$betas$trans,scorepar)
    res$mus$trans<-res$mus$trans[c(1:scorepar$n+scorepar$nsmall,1:scorepar$nsmall)]

  } else {
    for(i in 2:scorepar$slices-1) {
      res$betas[[i+1]]<-matrix(0,nrow=scorepar$paramsets[[i]]$n,ncol=scorepar$paramsets[[i]]$n)
      res$mus[[i+1]]<-vector()
      mainnodes<-scorepar$paramsets[[i]]$mainnodes
      for (j in mainnodes)  {
        parentnodes <- which(incidence$trans[,j]==1)
        res$mus[[i+1]][j]<-MAPparametersDBNmu(j,parentnodes, scorepar$paramsets[[i]])
        if(length(parentnodes)>0) {
          res$betas[[i+1]][parentnodes,j]<-MAPparametersDBNcore(j,parentnodes,scorepar$paramsets[[i]])
        }
      }
      res$betas[[i+1]]<-BiDAG:::DBNtransform(res$betas[[i+1]],scorepar)
      res$mus[[i+1]]<-res$mus[[i+1]][c(1:scorepar$n+scorepar$nsmall,1:scorepar$nsmall)]

    }
    slices<-scorepar$slices
    res$betas[[1]]<-matrix(0,nrow=scorepar$paramsets[[slices]]$n,ncol=scorepar$paramsets[[slices]]$n)
    res$mus[[1]]<-vector()
    mainnodes<-scorepar$paramsets[[slices]]$mainnodes
    for (j in mainnodes)  {
      parentnodes <- which(incidence$init[,j]==1)
      res$mus[[1]][j]<-MAPparametersDBNmu(j,parentnodes, scorepar$paramsets[[slices]])
      if(length(parentnodes)>0) {
        res$betas[[1]][parentnodes,j]<-MAPparametersDBNcore(j,parentnodes,scorepar$paramsets[[slices]])
      }
    }
  }
  return(res)
}
#function for computing MAE for one test sample, only works for 2 slices
DBNCVpar2<-function(DBNmaxmat,DBNconsmat=NULL,scorepar,Di,returnDi=FALSE){
  res<-vector()
  #estimate parameters
  if(!is.null(DBNconsmat)) {
    ord<-tryCatch(orderdagDBN(DBNconsmat))
    if(is.character(ord)) {
      DBNmat<-1*(DBNconsmat&DBNmaxmat)
      ord<-orderdagDBN(DBNmat)
    } else {
      DBNmat<-DBNconsmat
    }
  } else {
    DBNmat<-DBNmaxmat
    ord<-orderdagDBN(DBNmat)
  }

  mappar<-MAPparametersDBN(DBNmat,scorepar)
  b<-scorepar$bgn
  slices<-2
  nsmall<-scorepar$nsmall
  DBNmat[1:(b+nsmall),1:(b+nsmall)]<-0


  #take Di at point 0 and propagate to the next time slice
  staticnodes<-c()
  if(b>0) staticnodes<-Di[1:b]
  firstslice<-Di[1:(b+nsmall)]
  estDi<-c()
  ord2<-ord[-which(ord%in%c(1:(b+nsmall)))]
  prevslice<-firstslice
  curslice<-vector(length=nsmall) #create empty vector for the new slice
  curslice2<-c(prevslice,curslice) #append to the previous slice with estimated values
  for(i in ord2) {
    curslice2[i]<-mappar$mus$trans[i]
    parentnodes<-which(mappar$betas$trans[,i]!=0)
    if(length(parentnodes)>0) {
      coefs<-mappar$betas$trans[parentnodes,i]
      coefval<-unlist(curslice2[parentnodes])
      curslice2[i]<-curslice2[i]+sum(as.numeric(coefval*coefs))
    }
  }
  estDi<-c(estDi,curslice2[1:nsmall+b+nsmall])
  prevslice<-c(staticnodes,curslice2[1:nsmall+b+nsmall])

  res<-vector()
  res[1]<-sum(abs(Di[(b+nsmall+1):length(Di)]-estDi[1:length(estDi)]))
  res[2]<-sum((Di[(b+nsmall+1):length(Di)]-estDi[1:length(estDi)])^2)

  names(res)<-c("abs","sq")
}
#function for computing MAE for one test sample, works for multiple slices
DBNCVparmult2<-function(DBNmaxmat,DBNconsmat=NULL,scorepar,Di){
  res<-vector()
  #estimate parameters
  if(!is.null(DBNconsmat)) {
    ord<-tryCatch(orderdagDBN(DBNconsmat))
    if(is.character(ord)) {
      DBNmat<-1*(DBNconsmat&DBNmaxmat)
      ord<-orderdagDBN(DBNmat)
    } else {
      DBNmat<-DBNconsmat
    }
  } else {
    DBNmat<-DBNmaxmat
    ord<-orderdagDBN(DBNmat)
  }

  mappar<-MAPparametersDBN(DBNmat,scorepar)
  b<-scorepar$bgn
  slices<-2
  nsmall<-scorepar$nsmall
  DBNmat[1:(b+nsmall),1:(b+nsmall)]<-0


  #take Di at point 0 and propagate to the next time slice
  staticnodes<-c()
  if(b>0) staticnodes<-Di[1:b]
  firstslice<-Di[1:(b+nsmall)]
  estDi<-c()
  ord2<-ord[-which(ord%in%c(1:(b+nsmall)))]
  prevslice<-firstslice
  curslice<-vector(length=nsmall) #create empty vector for the new slice
  curslice2<-c(prevslice,curslice) #append to the previous slice with estimated values
  for(i in ord2) {
    curslice2[i]<-mappar$mus$trans[i]
    parentnodes<-which(mappar$betas$trans[,i]!=0)
    if(length(parentnodes)>0) {
      coefs<-mappar$betas$trans[parentnodes,i]
      coefval<-unlist(curslice2[parentnodes])
      curslice2[i]<-curslice2[i]+sum(as.numeric(coefval*coefs))
    }
  }
  estDi<-c(estDi,curslice2[1:nsmall+b+nsmall])
  prevslice<-c(staticnodes,curslice2[1:nsmall+b+nsmall])

  res<-vector()
  res[1]<-sum(abs(Di[(b+nsmall+1):length(Di)]-estDi[1:length(estDi)]))
  res[2]<-sum((Di[(b+nsmall+1):length(Di)]-estDi[1:length(estDi)])^2)

  names(res)<-c("abs","sq")
  return(res)
}
DBNCVest<-function(DBNmat,DBNparmus, DBNparbetas, Di, ord, b=0){
  res<-vector()
  #estimate parameters

  nsmall<-(ncol(DBNmat)-b)/2
  DBNmat[1:(b+nsmall),1:(b+nsmall)]<-0


  #take Di at point 0 and propagate to all time slices in an order defined by the DAG
  staticnodes<-c()
  if(b>0) staticnodes<-Di[1:b]
  ord2<-ord[-which(ord%in%c(1:(b+nsmall)))]
  curslice<-rep(0,nsmall) #create empty vector for the new slice
  curslice2<-unlist(c(Di,curslice))#append to the previous slice with estimated values
  for(i in ord2) {
    curslice2[i]<-DBNparmus[i]
    parentnodes<-which(DBNparbetas[,i]!=0)
    if(length(parentnodes)>0) {
      coefs<-DBNparbetas[parentnodes,i]
      coefval<-unlist(curslice2[parentnodes])
      curslice2[i]<-curslice2[i]+sum(as.numeric(coefval*coefs))
    }
  }
  return(curslice2[1:nsmall+b+nsmall])
}
DBNCVparmult<-function(DBNmaxmat,DBNconsmat=NULL,scorepar,Di){
  b<-scorepar$bgn
  slices<-scorepar$slices
  nsmall<-scorepar$nsmall
  stationary<-scorepar$stationary

  res<-vector()
  if(!is.null(DBNconsmat)) {
    ord<-tryCatch(orderdagDBN(DBNconsmat))
    if(is.character(ord)) {
      DBNmat<-1*(DBNconsmat&DBNmaxmat)
      ord<-orderdagDBN(DBNmat)
    } else {
      DBNmat<-DBNconsmat
    }
  } else {
    DBNmat<-DBNmaxmat
    ord<-orderdagDBN(DBNmat)
  }

  DBNpar<-MAPparametersDBN(DBNmat,scorepar)

  if(stationary) {
    newDi<-Di[1:nsmall]
    for(i in 2:slices) {
      newDi<-DBNCVest(DBNmat,DBNpar$mus$trans, DBNpar$betas$trans, newDi,ord,b=b)
      res[[i-1]]<-sum(abs(newDi-Di[1:nsmall+nsmall*(i-1)]))
    }
  } else {
    newDi<-Di[1:nsmall]
    for(i in 2:slices) {
      newDi<-DBNCVest(DBNmat,DBNpar$mus[[i]], DBNpar$betas[[i]], newDi, ord, b=b)
      res[[i-1]]<-sum(abs(newDi-Di[1:nsmall+nsmall*(i-1)]))
    }
  }
  return(res)
}
DBNCVpar<-function(DBNmaxmat,DBNconsmat=NULL,dbndata2,Di,slices=2,b=0){
 betamatrix<-DBNmaxmat
 betamatrix[1:ncol(DBNmaxmat),1:ncol(DBNmaxmat)]<-0
 res<-vector()
  #estimate parameters
  if(!is.null(DBNconsmat)) {
    ord<-tryCatch(orderdagDBN(DBNconsmat))
    if(is.character(ord)) {
      DBNmat<-1*(DBNconsmat&DBNmaxmat)
      ord<-orderdagDBN(DBNmat)
    } else {
      DBNmat<-DBNconsmat
    }
  } else {
    DBNmat<-DBNmaxmat
    ord<-orderdagDBN(DBNmat)
  }

  nsmall<-(ncol(DBNmat)-b)/2
  DBNmat[1:(b+nsmall),1:(b+nsmall)]<-0
  colnames(dbndata2)<-make.names(colnames(dbndata2), allow_=F)
  colnames(DBNmat)<-rownames(DBNmat)<-make.names(colnames(DBNmat),allow_=F)
  DBNparam<-DBNparamfit(DBNmat,dbndata2)
  rownames(betamatrix)<-rownames(DBNmat)
  #take Di at point 0 and propagate to all time slices in an order defined by the DAG
  nodes<-colnames(dbndata2)
  staticnodes<-c()
  if(b>0) staticnodes<-Di[1:b]
  firstslice<-Di[1:(b+nsmall)]
  fnames<-nodes[1:(b+nsmall)]
  names(firstslice)<-fnames
  snames<-nodes[1:nsmall+b+nsmall]
  estDi<-c()
  ord2<-ord[-which(ord%in%c(1:(b+nsmall)))]
  prevslice<-firstslice
  muvec<-vector()
  for(s in 2:slices) {
    curslice<-vector(length=nsmall) #create empty vector for the new slice
    names(curslice)<-snames #name it
    curslice2<-c(prevslice,curslice) #append to the previous slice with estimated values
    for(i in ord2) {
        curmodel<-DBNparam[[nodes[i]]]
        curslice2[i]<-curmodel$coefficients[1]
        muvec[i]<-curmodel$coefficients[1]
        #print(i)
        #print(curmodel$coefficients)
        if(length(curmodel$coefficients)>1) {
        coefs<-curmodel$coefficients[-1]
        coefval<-unlist(curslice2[names(curmodel$coefficients)[-1]])
        betamatrix[names(curmodel$coefficients)[-1],i]<-unlist(curslice2[names(curmodel$coefficients)[-1]])
        curslice2[i]<-curslice2[i]+sum(as.numeric(coefval*coefs))
        }
    }
    estDi<-c(estDi,curslice2[1:nsmall+b+nsmall])
    prevslice<-c(staticnodes,curslice2[1:nsmall+b+nsmall])
    names(prevslice)<-fnames
  }

  res<-vector()
  res[1]<-sum(abs(Di[(b+nsmall+1):length(Di)]-estDi[1:length(estDi)]))
  res[2]<-sum((Di[(b+nsmall+1):length(Di)]-estDi[1:length(estDi)])^2)

  names(res)<-c("abs","sq")
  return(res)
}
orderdagDBN<-function(adj) {
  n<-ncol(adj)
  allnodes<-c(1:n)
  curnodes<-c(1)
  order<-c()
  cntr<-1
  while(length(curnodes)<n & cntr<n) {
    npar<-apply(adj,2,sum)
    curnodes<-which(npar==0)
    order<-c(order,setdiff(curnodes,order))
    adj[curnodes,]<-0
    cntr<-cntr+1
  }

  if(sum(adj)==0) return(order)
  else warning("not a DAG")

}
MAPparametersDBNcore<-function(j,parentnodes, param) {
  lp<-length(parentnodes) # number of parents

  Sigma <- param$SigmaN
  A <- Sigma[j,j]
  if(lp==0){# no parents
    return(numeric(0))
  } else {
    D <- as.matrix(Sigma[parentnodes,parentnodes])
    choltemp<-chol(D)
    B <- Sigma[j,parentnodes]
    C <- backsolve(choltemp,B,transpose=TRUE)
    E <- backsolve(choltemp,C) #computing betas
    myE<-B%*%solve(D) #same as E this is how it is done in formulas
    return(myE)
  }
}
MAPparametersDBNmu<-function(j,parentnodes, param) {
  lp<-length(parentnodes) # number of parents

  m<-param$muN[j]
  Sigma <- param$SigmaN
  A <- Sigma[j,j]
  if(lp>0){
    D <- as.matrix(Sigma[parentnodes,parentnodes])
    choltemp<-chol(D)
    B <- Sigma[j,parentnodes]
    C <- backsolve(choltemp,B,transpose=TRUE)
    myE<-B%*%solve(D) #same as E this is how it is done in formulas
    m<-m-B%*%solve(D)%*%param$muN[parentnodes]
  }
  return(m)
}
DBNparamfit<-function(dag,data) {

  res<-list()
  data<-as.data.frame(data)
  ordery<-orderdagDBN(dag)
  nodes<-colnames(data)
  for(i in ordery) {
    pari<-which(dag[,i]==1)
    vars<-colnames(dag)[pari]

    if(length(pari)>0) {
      parents<-paste(colnames(data)[pari], collapse="+")
      fm <- as.formula(paste(colnames(data)[i], "~", parents))
      res[[i]]<-lm(fm, data)
    } else {
      res[[i]]<-lm(data[,i] ~ 1)
    }
  }
  names(res)<-nodes
  res
}
