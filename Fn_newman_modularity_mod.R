#ref: Gomez S., Jensen P. & Arenas A. (2009). Analysis of community structure in networks of correlated data. PHYS REV E, 80, 016114.
#Newman M.E. (2006). Finding community structure in networks using the eigenvectors of matrices. PHYS REV E, 74, 036104.

#membership:  a vector

newman_modularity_mod<-function(adj, membership, decomp=F){
  if(!is.matrix(adj)){stop("The input must be a matrix")}
  if(!isSymmetric(unname(adj))){stop("The input matrix must be symmetric")}
  
  n<-nrow(adj)
  A0<-adj
  k<-colSums(A0)
  m<-sum(k)/2
  
  n.m<-length(unique(membership))
  
  delta<-matrix(0, n, n)
  for(i in 1:n.m){
    tmp<-which(membership==i)
    delta[tmp,tmp]<-1
  }
  
  A0.pos<-A0; A0.pos[A0.pos<0]=0
  A0.neg<-A0; A0.neg[A0.neg>0]=0
  A0.neg<-(-A0.neg)
  
  k.pos<-colSums(A0.pos)
  m.pos<-sum(k.pos)/2
  k.neg<-colSums(A0.neg)
  m.neg<-sum(k.neg)/2
  
  if(m.pos==0){x1<-0 }else{ x1<-k.pos%o%k.pos/2/m.pos}
  if(m.neg==0){x2<-0 }else{ x2<-k.neg%o%k.neg/2/m.neg}
  Q<-(A0-x1+x2)*delta
  
  if(decomp==F){return(sum(Q)/2/(m.pos+m.neg))}else{
    Q.decomp.mod<-rep(NA, n.m)
    for(i in 1:n.m){
      tmp<-which(membership==i)
      Q.decomp.mod[i]<-sum(Q[tmp,tmp])/2/(m.pos+m.neg)
    }
    Q.decomp.node<-(rowSums(Q)+colSums(Q))/4/(m.pos+m.neg)
    Q.decomp.node.rescale<-(Q.decomp.node-min(Q.decomp.node))/diff(range(Q.decomp.node))
    return(list(Q=sum(Q)/2/(m.pos+m.neg), Q.decomp.mod=Q.decomp.mod, 
                Q.decomp.node=Q.decomp.node, Q.decomp.node.rescale=Q.decomp.node.rescale))
  }
  
}