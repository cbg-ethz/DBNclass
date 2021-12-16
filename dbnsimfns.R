makeDBNblacklist<-function(n) {
  bl<-matrix(nrow=n*n*2,ncol=2)
  colnames(bl)<-c("from","to")
  for(i in 1:n) {
    bl[1:n+(i-1)*n,"from"]<-paste("X",i,".2",sep="")
    bl[1:n+(i-1)*n,"to"]<-paste("X",1:n,sep="")
  }
  k<-max(1:n+(i-1)*n)+1
  for(i in 1:n) {
    for(j in 1:n) {
      bl[k,"from"]<-paste("X",i,sep="")
      bl[k,"to"]<-paste("X",j,sep="")
      k<-k+1
    }
  }
  return(bl)
}
modellistDBN<-function(pedges, p=c(0.5,0.6,0.7,0.8,0.9,0.99),n) {
  res<-list()
  n<-ncol(pedges)
  k<-1
  for(i in p) {
    res[[k]]<-matrix(0,nrow=n,ncol=n)
    ones<-which(pedges>i)
    res[[k]][ones]<-1
    res[[k]][,1:(n/2)]<-0
    k<-k+1
  }
  return(res)
}
makeAmatBoot<-function(bootresult,n) {
  res<-matrix(0,nrow=2*n,ncol=2*n)
  colnames(res)<-rownames(res)<-c(paste("X",1:n,sep=""),paste("X",1:n,".2",sep=""))
  for(i in 1:nrow(bootresult)) {
    res[bootresult$from[i],bootresult$to[i]]<-bootresult$strength[i]
  }
  return(res)
}
DBNCV2.sim<-function(mapmodel,consmodels,testdata,localscore,n){
  if(!is.null(consmodels)) {
    resmat<-matrix(nrow=1+length(consmodels),ncol=nrow(testdata))
    for(i in 1:nrow(testdata)) {
      resmat[1,i]<-DBNCVpar2(mapmodel, NULL,localscore,testdata[i,])[1]/n
      for(j in 1:length(consmodels)) {
        resmat[j+1,i]<-DBNCVpar2(mapmodel, consmodels[[j]],localscore,testdata[i,])[1]/n
      }
    }
    #return(resmat)
    return(apply(resmat,1,mean))
  } else {
    resvec<-vector()
    for(i in 1:nrow(testdata)) {
      resvec[i]<-DBNCVpar2(mapmodel, NULL,localscore,testdata[i,])[1]/n
    }
    return(mean(resvec))
  }
}
DBNCVmult.sim<-function(mapmodel,consmodels,testdata,localscore,n){
  if(!is.null(consmodels)) {
    resmat<-matrix(nrow=1+length(consmodels),ncol=nrow(testdata))
    for(i in 1:nrow(testdata)) {
      resmat[1,i]<-DBNCVparmult2(mapmodel, NULL,localscore,testdata[i,])[1]/n
      for(j in 1:length(consmodels)) {
        resmat[j+1,i]<-DBNCVparmult2(mapmodel, consmodels[[j]],localscore,testdata[i,])[1]/n
      }
    }
    #return(resmat)
    return(apply(resmat,1,mean))
  } else {
    resvec<-vector()
    for(i in 1:nrow(testdata)) {
      resvec[i]<-DBNCVparmult2(mapmodel, NULL,localscore,testdata[i,])[1]/n
    }
    return(mean(resvec))
  }
}
generateDBN<-function(n=40,seed=100){
  set.seed(seed)
  dag1<-randomDAG(n=n,prob=2/n,lB=0.4, uB=2, V=paste("X",1:n,sep=""))
  set.seed(seed)
  dag2<-randomDAG(n=n,prob=2/n,lB=0.4, uB=2, V=paste("X",1:n,".2",sep=""))
  interedges<-generateInter(n)
  return(merge2DBN(dag1,dag2,interedges))
}
generateInter<-function(n) {
  edges2self<-sample.int(n,n/2)
  otheredges.target<-c(1:n)[-edges2self]
  otheredges.source<-sample.int(n,length(otheredges.target))
  edgelist<-list()
  for(i in 1:n) {
    edgelist[[i]]<-integer(0)
    if(i%in%edges2self) edgelist[[i]]<-i+n
    if(i%in%otheredges.target) {
      whichi<-which(otheredges.target==i)
      edgelist[[i]]<-c(edgelist[[i]],otheredges.source[whichi]+n)
    }
  }
  names(edgelist)<-paste("X",1:n,sep="")
  weights<-runif(n,min=0.4,max=2)
  weightlist<-list()
  for(i in 1:n) {
    weightlist[[i]]<-list()
    weightlist[[i]]$weight<-weights[i]
    names(weightlist)[[i]]<-paste("X",i,"|","X",edgelist[[i]]-n,".2",sep="")
  }
  res<-list()
  res$edgeslist<-edgelist
  res$weightlist<-weightlist
  return(res)
}
merge2DBN<-function(dag1,dag2,interedges) {
  nsmall<-length(dag1@nodes)
  commonDBN<-dag1
  commonDBN@nodes<-c(dag1@nodes,dag2@nodes)
  transpars<-dag2@edgeL
  for(i in 1:length(transpars)) {
    if(length(transpars[[i]])>0) transpars[[i]]$edges<-transpars[[i]]$edges+nsmall
    #take union with daginter
    commonDBN@edgeL[[i]]$edges<-c(commonDBN@edgeL[[i]]$edges,interedges$edgeslist[[i]])
  }
  commonDBN@edgeL<-c(commonDBN@edgeL,transpars)
  commonDBN@edgeData@data<-c(dag1@edgeData@data,interedges$weightlist,dag2@edgeData@data)
  return(commonDBN)
}
genDBN<-function(d,n,slices=2,lB=0.1, uB=1, wmpct=0, wmmin=1,wmmax=3, shdpct=0, wm=NULL){
  intstr<-genDAG(d,n)
  dbnmat<-matrix(0,nrow=n*slices,ncol=n*slices)
  wmmat<-matrix(0,nrow=n*slices,ncol=n*slices)
  for(i in 1:slices) {
    dbnmat[1:n+(i-1)*n,1:n+(i-1)*n]<-intstr$dag
    wmmat[1:n+(i-1)*n,1:n+(i-1)*n]<-intstr$wm
  }

  edges2self<-sample.int(n,n/2)
  transwself<-runif(n/2,min=lB,max=uB)
  otheredges.target<-c(1:n)[-edges2self]
  otheredges.source<-sample.int(n,length(otheredges.target))
  transwother<-runif(length(otheredges.target),min=lB,max=uB)

  trans<-matrix(0,nrow=n,ncol=n)
  wm<-matrix(0,nrow=n,ncol=n)

  for(i in 1:length(edges2self)) {
    trans[edges2self[i],edges2self[i]]<-1
    wm[edges2self[i],edges2self[i]]<-transwself[i]
  }
  for(i in 1:length(otheredges.target)) {
    trans[otheredges.source[i],otheredges.target[i]]<-1
    wm[otheredges.source[i],otheredges.target[i]]<- transwother[i]

  }

  for(i in 1:(slices-1)) {
    dbnmat[1:n+(i-1)*n,1:n+i*n]<-trans
    wmmat[1:n+(i-1)*n,1:n+i*n]<-wm
  }

  res<-list()
  res$dbn<-dbnmat
  res$wm<-wmmat
  return(res)
}
genDataDBN<-function(DBN,slices=2,ss=100) {
  dbn<-DBN$dbn
  wm<-DBN$wm
  n<-ncol(dbn)
  nsmall<-n/2
  npar<-apply(dbn,2,sum)

  means<-rep(0,n)
  orderx<-orderdag(dbn)
  sigmas<-rep(0.5,n)

  #first generate weight matrix
  edges<-which(dbn!=0)
  nedges<-length(edges)
  #define order of a dag
  ordery<-rev(orderx)


  #generate data
  datas<-matrix(nrow=ss, ncol=n)
  for(i in ordery) {
    if(npar[i]==0) {
      datas[,i]<-rnorm(ss,mean=means[i])
    } else {
      pari<-as.vector(which(dbn[,i]!=0))
      datas[,i]<-wm[pari,i] %*% t(datas[,pari])+rnorm(ss,mean=means[i],sd=sigmas[i])
    }
  }


  ordery<-ordery[-which(ordery%in%c(1:nsmall))]
  dbn[,1:nsmall]<-0
  wm[,1:nsmall]<-0

  datal<-list()
  datal[[1]]<-datas
  for(k in 2:(slices-1)) {
    datal[[k]]<-cbind(datal[[k-1]][,1:nsmall+nsmall],matrix(0,nrow=ss,ncol=nsmall))
    for(i in ordery) {
      if(npar[i]==0) {
        datal[[k]][,i]<-rnorm(ss,mean=means[i])
      } else {
        pari<-as.vector(which(dbn[,i]!=0))
        datal[[k]][,i]<-wm[pari,i] %*% t(datal[[k]][,pari])+rnorm(ss,mean=means[i],sd=sigmas[i])
      }
    }
  }
  if(slices>2) {
    for(k in 2:(slices-1)) {
      datas<-cbind(datas,datal[[k]][,1:nsmall+nsmall])
    }
  }

  res<-list()
  res$gs<-Reduce("rbind",datal)
  res$mcmc<-datas
  return(res)
}
genDAG<-function(d, n, lB=0.1, uB=1, wmpct=0, wmmin=1,wmmax=3, shdpct=0, wm=NULL) {
  res<-list()
  if(shdpct>0) wmpct<-0
  if(is.null(wm)) {
    dag<-graph2m(pcalg::randomDAG(n, d*2/n, V=as.character(1:n)))
    edges<-which(dag==1)
    nedges<-sum(dag)
    wm<-dag
    edgeWeights<-runif(nedges, lB, uB)
    wm[edges]<-edgeWeights
  } else {
    dag<-matrix(0,nrow=n,ncol=n)
    edges<-which(abs(wm)>lB)
    dag[edges]<-1
    nedges<-length(edges)

    if(wmpct>0) {
      nme<-ceiling(wmpct*nedges/100)
      newedges<-sample(edges,nme)
      nweights<-wm[edges]
      nweights[newedges]<-sapply(nweights, function(x)x^(sample(c(1,-1),1)*runif(1,min=wmmin,max=wmmax)))
      wm[edges]<-sapply(nweights,max,lB)
    } else if(shdpct>0) {
      shdy<-max(ceiling(shdpct*nedges/100),2)
      newdag<-modifydag(dag,shdy) #get new dag with a defined shd with the start one
      newedges<-which(newdag==1) #edges of the new dag
      commonedges<-which((1*(dag & newdag))==1) #edges which are common in old/new dags
      newedges<-setdiff(newedges,commonedges)
      newn<-length(newedges)
      newwm<-matrix(0,nrow=n,ncol=n)
      newwm[commonedges]<-wm[commonedges]
      newwm[newedges]<-runif(newn, lB, uB)
      res$dag<-newdag
      res$wm<-newwm
      return(res)
    }
  }
  res$dag<-dag
  res$wm<-wm
  return(res)
}
orderdag<-function(adj) {
  n<-ncol(adj)
  allnodes<-c(1:n)
  curnodes<-c(1)
  order<-c()
  cntr<-1
  while(length(curnodes)<n & cntr<n) {
    npar<-apply(adj,2,sum)
    curnodes<-which(npar==0)
    order<-c(setdiff(curnodes,order),order)
    adj[curnodes,]<-0
    cntr<-cntr+1
  }

  if(sum(adj)==0) return(order)
  else stop("not a DAG")

}

