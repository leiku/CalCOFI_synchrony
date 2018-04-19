#' wavelet multiple linear regression models
#' 
#' @param trans0f the response variable, a time x frequencies wavelet transform (or concatenated set of wavelet transforms, with dims  (time x locations) x frequencies)
#' @param transnf a list of predictor variables containing wavelet transforms having the same dimensions as the response variable
#' @param n the number of predictor variables (i.e., the number of elements in the list transnf)
#' @param freqs a vector of indices (rows) corresponding to focal frequencies 
#'
#' @return \code{wmodel} returns a frequencies x predictor variables matrix of complex coefficients
#' @note Developed by Lawrence Sheppard and Daniel Reuman; R code by Jonathan Walter.
#' @author Jonathan Walter, \email{jonathan.walter@ku.edu}; Lawrence Sheppard, \email{lwsheppard@ku.edu}; Daniel Reuman, \email{reuman@ku.edu}

#' @examples
#' #Generate some artificial data, fit models, measure goodness of fit
#' library(biwavelet) #package for computing wavelet transforms
#' 
#' time=100
#' p1 = 10
#' p2 = 7
#'
#' sig1<-rep(sin(seq(0,2*pi,length.out=p1)[-p1]),length.out=time)+rnorm(time,0,0.1)
#' sig2<-rep(sin(seq(0,2*pi,length.out=p2)[-p2]),length.out=time)+rnorm(time,0,0.1)
#' sig3<-sig1+sig2+rnorm(time,0,0.1) #sig3 is the sum of artificial signals sig1 and sig2
#' 
#' par(mfrow=c(3,1)) #plot these 
#' plot(sig1, type="l")
#' plot(sig2, type="l")
#' plot(sig3, type="l")
#' 
#' trans0f<-wt(d=cbind(1:time,sig3), dj=1/20, do.sig=F)$wave
#' transnf<-list(wt(d=cbind(1:time,sig1), dj=1/20, do.sig=F)$wave,
#'              wt(d=cbind(1:time,sig2), dj=1/20, do.sig=F)$wave)
#'              
#' X1<-wmodel(trans0f,transnf,n=2,freqs=1:nrow(trans0f)) #compute coefficients
#' matX1_1<-matrix(X1[,1],nrow(trans0f),time) # reshape into matrices of coefficients
#' matX1_2<-matrix(X1[,2],nrow(trans0f),time)
#' model<-transnf[[1]]*matX1_1 + transnf[[2]]*matX1_2 #calculate predicted wavelet transform
#' model.resids<-trans0f-model #calculate model residuals
#' d2<-mean(t(model.resids)*Conj(t(model.resids)),na.rm=T) #calculate squared residuals to assess model fit; smaller is better
#' print(d2)
#' 
#' par(mfrow=c(2,1)) #plot the output: these should be close but not identical.
#' image(Re(trans0f))
#' image(Re(model))
#' @export
#' 

wmodel.test<-function(trans0f, transnf, n, freqs){ ## wrapper function to loop through frequencies and remove NA's (e.g., scallopping)
  if(n==1){stop("Error: n should be higher than 1")}
  
  pownf<-vector("list", n)
  for(nn in 1:n){
    pownf[[nn]]<-colMeans(transnf[[nn]]*Conj(transnf[[nn]]), na.rm=T)
  }
  
  nm=expand.grid(1:n,1:n)
  nm<-nm[nm[,1]!=nm[,2],]
  nm<-nm[nm[,1]<nm[,2],]
  A<-array(1, dim=c(n,n,length(freqs)))
  for(ii in 1:nrow(nm)){
    nn<-nm[ii,1]
    mm<-nm[ii,2]
    A[nn,mm,]<-colMeans(transnf[[mm]]*Conj(transnf[[nn]]),na.rm=T)/pownf[[nn]]
    A[mm,nn,]<-Conj(A[nn,mm,])*pownf[[nn]]/pownf[[mm]]
  }
  
  B<-matrix(NA, n, length(freqs))
  for(nn in 1:n){
    B[nn,]<-colMeans(trans0f*Conj(transnf[[nn]]),na.rm=T)/pownf[[nn]]
  }
  
  X=matrix(NA, length(freqs), n)
  for(j in freqs){
    X[j,]<-Conj(solve(Conj(A[,,j]),Conj(B[,j])))
  }
  return(X)
}

