#ref: Gomez S., Jensen P. & Arenas A. (2009). Analysis of community structure in networks of correlated data. PHYS REV E, 80, 016114.
#Newman M.E. (2006). Finding community structure in networks using the eigenvectors of matrices. PHYS REV E, 74, 036104.

newman_eigenvector_mod<-function(adj, plot.dendrogram=F){
  #function to test if a graph is connected; idea by Ed Scheinerman, circa 2006
  #source: http://www.ams.jhu.edu/~ers/matgraph/; routine: matgraph/@graph/isconnected.m
  is.connect<-function(adj1){
    adj1<-as.matrix(adj1)
    if(nrow(adj1)<2){return(FALSE)}else{
      if(length(which(colSums(adj1)==0))>0){return(FALSE)}else{
        x<-c(1,rep(0, nrow(adj1)-1))
        while(1){
          y<-x
          x<-adj1%*%x + x
          x1<-rep(0,length(x))
          x1[x>0]<-1
          x<-x1
          if(all(x==y)){break}
        }
        if(sum(x)<length(x)){return(FALSE)}else{return(TRUE)}
      } 
    }
  }
  
  if(!is.matrix(adj)){stop("The input must be a matrix")}
  if(!isSymmetric(adj)){stop("The input matrix must be symmetric")}
  
  A0<-adj
  n<-nrow(A0)
  
  A0.pos<-A0; A0.pos[A0.pos<0]=0
  A0.neg<-A0; A0.neg[A0.neg>0]=0
  A0.neg<-(-A0.neg)

  k.pos<-colSums(A0.pos)
  m.pos<-sum(k.pos)/2
  k.neg<-colSums(A0.neg)
  m.neg<-sum(k.neg)/2
  
  
  if(m.pos==0){tmp1<-0}else{tmp1<-k.pos %o% k.pos/2/m.pos}  #<<<<<<<<
  if(m.neg==0){tmp2<-0}else{tmp2<-k.neg %o% k.neg/2/m.neg}  #<<<<<<<<
  B0<- A0 - tmp1 + tmp2   #<<<<<<<<<<<
  
  modules<-rep(1, n)
  Queue<-vector("list", n)    #record the temporal divisions
  Queue[[1]]<-modules
  a<-1
  r<-2   # index of loops
  
  #main function
  while(a<=max(modules)){  # while there is always a divisible subgraph
    
    # compute modularity matrix
    i.remain<-which(modules==a) # nodes in current (sub)graph to partition
    n1<-length(i.remain)
    
    A1<-A0[i.remain,i.remain]
    
    #k1.pos<-colSums(A0.pos[i.remain, i.remain])
    #k1.neg<-colSums(A0.neg[i.remain, i.remain])
    #tmp1.pos<- k.pos[i.remain]%o%k.pos[i.remain]/2/m.pos
    #tmp1.neg<- k.neg[i.remain]%o%k.neg[i.remain]/2/m.neg
    #tmp2.pos<- k1.pos-k.pos[i.remain]*sum(k1.pos[i.remain])/2/m.pos
    #tmp2.neg<- (k1.neg-k.neg[i.remain]*sum(k1.neg[i.remain])/2/m.neg)  
    #if(m.pos==0){tmp1.pos<-0; tmp2.pos<-0}
    #if(m.neg==0){tmp1.neg<-0; tmp2.neg<-0}
    
    #B1<-A1 - tmp1.pos + tmp1.neg
    #B0<- tmp2.pos - tmp2.neg        
    #B<-B1
    #diag(B)<-diag(B)-B0
    
    B<-B0[i.remain,i.remain]  #first part of Eq. 51 in Newman 2006
    diag(B)<-diag(B)-colSums(B)  #minus second part of Eq. 51 in Newman 2006
    
    E<-eigen(B)
    
    #if indivisible, terminate and check the next queue
    if(max(E$values)<=1e-5){
      a<-a+1
      next
    } 
    
    #if delta_Q < 0, terminate and check the next one
    i.max<-which.max(E$values)
    v1<-E$vectors[,i.max]
    i.pos<-which(v1>0)
    i.neg<-which(v1<0)
    if(sum(B[i.pos,i.pos])+sum(B[i.neg,i.neg])<=0){
      a<-a+1
      next
    }
    
    
    A1[A1>0]<-1
    A1[A1<0]<-0
    #if either of the subgraphs is disconnected, terminate and check the next queue 
    #if(vertex_connectivity(induced_subgraph(G1.pos, i.pos))<=0 | vertex_connectivity(induced_subgraph(G1.pos, i.neg))<=0){ 
    if(!is.connect(A1[i.pos,i.pos]) | !is.connect(A1[i.neg,i.neg])){
      a<-a+1
      next
    }
    modules[modules>=a]<-modules[modules>=a]+1
    modules[i.remain[i.pos]]<-modules[i.remain[i.pos]]-1
    Queue[[r]]<-modules
    r<-r+1
  }
  Queue<-Queue[1:(r-1)]
  M<-list(membership=modules, membership.temp=Queue)
  
  #plot dendrogram
  if(plot.dendrogram){
    names.node<-colnames(adj)
    if(is.null(names.node)){names.node<-1:n}
    n.step<-length(Queue)
    if(n.step==0){stop("Not able to generate dendrogram because of no division")}
    
    M.temp<-rep(NA, n.step)
    Queue1<-Queue
    ii<-order(modules)
    X<-vector("list", n.step)
    for(i in 1:n.step){
      M.temp[i]<-newman_modularity_mod(adj, Queue[[i]])
      Queue1[[i]]<-Queue[[i]][ii]
      i.module<-unique(Queue1[[i]])
      n.module<-length(i.module)
      x<-rep(NA, n.module)
      for(j in 1:n.module){
        x[j]<-mean(which(Queue1[[i]]==i.module[j]))
      }
      X[[i]]<-x
    }
    Y<-(-M.temp)
    plot(1:n, rep(min(Y)-0.1, n), type="p", pch=20, col="red", ylim=c(min(Y)-0.1, max(Y)), xaxt="n",
         xlab="nodes", ylab="-modularity")
    axis(1, at=1:n, labels=names.node[ii])
    
    for(i in 1:(n.step-1)){
      queue1<-Queue1[[i]]
      i.module<-unique(queue1)
      queue2<-Queue1[[i+1]]
      a<-0
      for(j in 1:length(X[[i]])){
        NN<-length(unique(queue2[queue1==i.module[j]]))
        for(r in 1:NN){
          a<-a+1
          lines(c(X[[i]][j], X[[i+1]][a]), c(Y[i], Y[i+1]), col="black", type="b")
        }
      }
    }
    queue1<-Queue1[[n.step]]
    i.module<-unique(queue1)
    
    a<-0
    for(i in 1:length(X[[n.step]])){
      for(j in 1:length(which(queue1==i.module[i]))){
        lines(c(X[[n.step]][i], j+a), c(Y[n.step], Y[n.step]-0.1), col="black", type="b")
      }
      a<-a+j
    }
    
  }
  
  return(M)
}
