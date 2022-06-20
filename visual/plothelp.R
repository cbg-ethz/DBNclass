plotGSE5462<-function(DBN,struct=c("init","trans"),b=0,highlight1,highlight2,ts,...){

  colnames(DBN)<-sub("\\.2","\\.",colnames(DBN))
  cols2<-c("#fd8d3c","#5ab4ac")
  dyn<-(ncol(DBN)-b)/2

  old.par<-par(no.readonly = TRUE)
  oldgraphpar<-graph.par()
  on.exit(par(old.par))
  on.exit(graph.par(oldgraphpar),add=TRUE)

  a<-d<-1.2
  c<-12

  if(!is.matrix(DBN)) {
    DBN<-graph2m(DBN)
  }

  nodelabs<-colnames(DBN)
  statcol<-"lightgrey"
  dyn1col<-"#f7f4f9"
  dyn2col<-"#d4b9da"


  if(struct=="init") {

    nodelabs<-nodelabs[1:(dyn+b)]

    if(b>0){
      legadj<-matrix(0,nrow=2,ncol=2)
      colnames(legadj)<-c("stat","1")
      legadj[1,2]<-1
      legendG<-m2graph(legadj)
      staticnames<-nodelabs[1:b]
      legcol<-c(statcol,dyn1col)
      legw<-rep(a,2)
      legh<-rep(d,2)
      legf<-rep(c,2)
      names(legcol)<-names(legw)<- names(legh)<-names(legf)<-c("stat","1")
    }

    dynamicnames<-nodelabs[1:dyn+b]

    adj<-DBN[1:(dyn+b),1:(dyn+b)]
    arcslist<-BiDAG:::adjacency2edgel(adj,nodes=nodelabs)
    graph.obj = new("graphNEL", nodes = nodelabs, edgeL = arcslist,
                    edgemode = 'directed')
    subGList<-list()
    sg<-list()

    graph.par(list(nodes=list(col=dyn1col, lty="solid", lwd=2, ...),graph=list(...,cex.main=1.5)))
    if(b!=0) {
      sg1 = subGraph(dynamicnames, graph.obj)
      sg2 = subGraph(staticnames, graph.obj)
      sgL = list(list(graph=sg1, cluster = TRUE),
                 list(graph=sg2, cluster = TRUE))
      graph.obj <- Rgraphviz::layoutGraph(graph.obj, subGList= sgL)
      graph::nodeRenderInfo(graph.obj)[["fill"]][staticnames] = statcol
      graph::nodeRenderInfo(graph.obj)[["fill"]][dynamicnames] = dyn1col

    } else {
      graph.obj <- Rgraphviz::layoutGraph(graph.obj)
    }

    if(b>0) {
      layout(matrix(c(1,1,1,1,1,3,3,
                      1,1,1,1,1,2,2,
                      1,1,1,1,1,3,3), nrow = 3, ncol = 7, byrow = TRUE))

      graph::nodeRenderInfo(legendG)[["fill"]]["stat"] = statcol
      graph::edgeRenderInfo(legendG)[["lwd"]]["stat~1"] = 0
      Rgraphviz::renderGraph(graph.obj )

      graph::plot(legendG,attrs=list(graph=list(rankdir="TB"),edge=list(lwd=0)),
                  nodeAttrs=list(fillcolor=legcol,height=legh,width=legw,fontsize=legf),
                  main="nodes (t):",cex.main=1.5)
    } else {
      graph.par(list(nodes=list(col=dyn1col, lty="solid", lwd=2, ...),graph=list(main="nodes (t):")))
      Rgraphviz::renderGraph(graph.obj )
    }
  }

  if(struct=="trans") {
    if(b>0) {
      legadj<-matrix(0,nrow=3,ncol=3)
      legadj[1,2]<-1
      legadj[2,3]<-1
      colnames(legadj)<-c("stat","i","i+1")
      legendG<-m2graph(legadj)
      legcol<-c(statcol,dyn1col,dyn2col)
      names(legcol)<-c("stat","i","i+1")
      legw<-rep(a,3)
      legh<-rep(d,3)
      legf<-rep(c,3)
      names(legw)<- names(legh)<-names(legf)<-c("stat","i","i+1")

      colvector = c(rep(statcol,b), rep(dyn1col,dyn),rep(dyn2col,dyn))
      shapevec = c(rep("triangle",b),rep("circle",dyn),rep("box",dyn))
      wvec=c(rep(7,b),rep(6,dyn),rep(5,dyn))
      names(colvector)<-names(shapevec)<-names(wvec)<-nodelabs
    } else {
      legadj<-matrix(0,nrow=2,ncol=2)
      legadj[1,2]<-1
      legw<-rep(a,2)
      legh<-rep(d,2)
      legf<-rep(c,2)
      colnames(legadj)<-c("i","i+1")
      legendG<-m2graph(legadj)
      legcol<-c(dyn1col,dyn2col)
      names(legcol)<-names(legw)<- names(legh)<-names(legf)<-c("i","i+1")
      colvector = c(rep(dyn1col,dyn),rep(dyn2col,dyn))
      names(colvector)<-nodelabs
      shapevec = c(rep("circle",dyn),rep("box",dyn))
      names(colvector)<-names(shapevec)<-nodelabs
    }

    adjt<-BiDAG:::DBNcut(DBN[1:(b+2*dyn),1:(b+2*dyn)],dyn,b)
    graph.obj<-m2graph(adjt)

    dyn1names<-nodelabs[1:dyn+b]
    dyn2names<-nodelabs[1:dyn+b+dyn]
    sgDyn1 = subGraph(dyn1names, graph.obj)
    sgDyn2 = subGraph(dyn2names, graph.obj)

    if(b>0) {
      staticnames<-nodelabs[1:b]
      sgStat = subGraph(staticnames, graph.obj)
      sgL = list(list(graph=sgStat, cluster = TRUE, attrs = c(rankdir="TB",rank="sink")),
                 list(graph=sgDyn1, cluster = TRUE, attrs = c(rank="same")),
                 list(graph=sgDyn2, cluster = TRUE, attrs = c(rank="same")))
    } else {
      sgL = list(list(graph=sgDyn1, cluster = TRUE, attrs = c(rank="same")),
                 list(graph=sgDyn2, cluster = TRUE, attrs = c(rank="same")))
    }

    graph.par(list(nodes=list(col=dyn1col, lty="solid", lwd=2, textsize=ts),graph=list()))
    graph.obj <- Rgraphviz::layoutGraph(graph.obj, subGList= sgL)

    if(b>0) graph::nodeRenderInfo(graph.obj)[["fill"]][staticnames] = statcol
    graph::nodeRenderInfo(graph.obj)[["fill"]][dyn1names] = dyn1col
    graph::nodeRenderInfo(graph.obj)[["col"]][dyn1names] = dyn1col

    graph::nodeRenderInfo(graph.obj)[["fill"]][dyn2names] = dyn2col
    graph::nodeRenderInfo(graph.obj)[["col"]][dyn2names] = dyn2col


    layout(matrix(c(1,1,1,1,1,3,3,
                    1,1,1,1,1,2,2,
                    1,1,1,1,1,3,3), nrow = 3, ncol = 7, byrow = TRUE))

    graph::nodeRenderInfo(graph.obj)[["col"]][highlight1] = cols2[1]
    graph::nodeRenderInfo(graph.obj)[["col"]][highlight2] = cols2[2]


    alledges12<-matrix(colnames(DBN)[which(DBN[1:dyn,]>0,arr.ind=TRUE)],ncol=2)
    alledgesf12 <- apply(alledges12, 1, paste, collapse = "~")
    graph::edgeRenderInfo(graph.obj)[["lwd"]]<-2
    graph::edgeRenderInfo(graph.obj)[["col"]]<-"#2c7bb6"

    for(i in 1:nrow(alledges12)) {
      graph::edgeRenderInfo(graph.obj)[["lwd"]][alledgesf12[i]] = 2
      graph::edgeRenderInfo(graph.obj)[["col"]][alledgesf12[i]] = "#d7191c"
    }

    nodelabs<-nodelabs[1:dyn]
    nodelabs2<-paste(nodelabs,".",sep="")
    alledges12<-matrix(c(nodelabs,nodelabs2),ncol=2)
    alledgesf12 <- apply(alledges12, 1, paste, collapse = "~")
    for(i in 1:nrow(alledges12)) {
      graph::edgeRenderInfo(graph.obj)[["lwd"]][alledgesf12[i]] = 2
      graph::edgeRenderInfo(graph.obj)[["col"]][alledgesf12[i]] = "grey"
    }

    #print(edgeRenderInfo(graph.plot)[["col"]])


    Rgraphviz::renderGraph(graph.obj)

    graph::plot(legendG,attrs=list(graph=list(rankdir="TB")),
                nodeAttrs=list(fillcolor=legcol,height=legh,width=legw,fontsize=legf),
                main="nodes (t):",cex.main=1.5)
  }
}
