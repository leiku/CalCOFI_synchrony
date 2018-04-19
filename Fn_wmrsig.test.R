#' Surrogate testing for wavelet multiple linear regression models
#' Lei made some modification to allow for 2 combined surrogate
#' 
#' Function to conduct surrogate testing of individual covariates within wavelet multiple regression models. Squared residuals are calculated for the actual model, and ranked against model fits (d2) using surrogate data, as a function of frequency
#' 
#' @param indata A list of data matrices (location X time)
#' @param s Index in indata of the covariate to surrogate (can be 2 numbers now)
#' @param r Index in indata of the response variable
#' @param n Number of covariates in the model
#' @param n.surrog Number of surrogates to test, defaults to 100
#' @param surr.test If \code{TRUE}, perform surrogate testing. Default is \code{FALSE}
#' @param surr.type what type of surrogates to use? Defaults to fft, choose aaft for the amplitude adjusted Fourier transform
#' @param tsranges A vector giving the min and max of a range of timescales of interest for sig testing, or else a matrix with rows giving the same.
#'
#' @return \code{wmrsig} returns a list consisting of:
#' \item{d2f}{squared residuals by frequency from model using real data}
#' \item{coefs}{The coefficients from the model}
#' \item{d2f.surr}{If surr.test equals \code{TRUE}, squared residuals by frequency from model using surrogate data for one predictor. Each column is a different surrogate}
#' \item{pval}{If surr.test equals \code{TRUE}, the p-value, as determined by ranking, as a function of frequency, the actual d2f against the surrogated d2f values}
#' \item{timescales}{a vector of timescales}
#' \item{sigtest}{p-values for timescale bands of interest, from \code{surrogtest}}
#' \item{pred.wt}{matrix of the predicted wavelet transforms of the model, in timescale X (location X time)}
#' 
#' @note Use AAFT surrogates if the variable you are testing cannot be normalized. Developed by Lawrence Sheppard and Daniel Reuman; R code by Thomas Anderson
#' @author Thomas Anderson, \email{anderstl@gmail.edu}; Lawrence Sheppard, \email{lwsheppard@ku.edu}; Daniel Reuman, \email{reuman@ku.edu}
#' @examples
#'x1<-matrix(rnorm(200,0,0.1),10,20)
#'xx1<-x1-rowMeans(x1) #center the variable
#'x2<-matrix(rnorm(200,0,0.1),10,20)
#'xx2<-x2-rowMeans(x2) #center the variable
#'x3<-matrix(rnorm(200,0,0.1),10,20)
#'xx3<-x3-rowMeans(x3) #center the variable
#'dat.list<-list(x1=xx1,x2=xx2,x3=xx3)
#'res<-wmrsig(indata = dat.list,r=1,s=2,n=2,n.surrog=100,surr.test = T)
#'print(res)
#' @export

wmrsig.twos1<-function(indata, s, r, n, n.surrog=100, surr.test=FALSE, surr.type="fft", tsranges=c(0,Inf)){
  d1<-sapply(X=indata,FUN=function(x){return(dim(x)[1])}) 
  d2<-sapply(X=indata,FUN=function(x){return(dim(x)[2])}) 
  if (!all(d1[2:length(d1)]==d1[1]) || !all(d2[2:length(d2)]==d2[1]))
  {
    stop("Error: all matrices must be same dimensions")
  }
  
  #flatten input data
  wt.covs<-list()
  for(i in 1:length(indata)){
    wt.covs[[i]]<-warray(indata[[i]],times=1:ncol(indata[[i]]),scale.min=2,scale.max.input=NULL,sigma=1.05,f0=1)$wave.array
  }
  
  timescales<-wt(indata[[i]][1,],times=1:ncol(indata[[i]]))$timescales
  
  covs.wide<-lapply(wt.covs, function(x) apply(x, 3, c)) #for each row, arranged like: site1year1, site2year1, ... site1year2, site2year2, ...
  
  resp<-covs.wide[[r]] #define response
  covs.list<-covs.wide[-r] #remove response from covariate list, thus defining covariates
  
  #make coefficients
  X1<-wmodel.test(trans0f = resp,transnf = covs.list,n=n,freqs = 1:ncol(resp)) #calculate coefficients
  
  model<-vector("list", ncol(X1))
  for(i in 1:ncol(X1)){
    model[[i]]<-t(t(covs.list[[i]])*X1[,i])
  }
  
  model1<-Reduce('+',model) #calculate predicted wavelet transform by summing predicted values
  model.resids<-resp-model1 #calculate model residuals
  d2f<-colMeans(model.resids*Conj(model.resids),na.rm=T) #calculate squared residuals by frequency to assess model fit; smaller is better
  d2f<-d2f[is.finite(d2f)]
 
  
  #To conduct surrogate testing
  if(surr.test){ 
    if(surr.type=="fft"){surrog<-ffts(rbind(indata[[s[1]]],indata[[s[2]]]), nsurrogs=n.surrog)} #make surrogates  <<---------
    if(surr.type=="aaft"){surrog<-aafts(rbind(indata[[s[1]]],indata[[s[2]]]), nsurrogs=n.surrog)} # <<-------------
    
    d2f.surr<-matrix(ncol=n.surrog,nrow=ncol(resp))
    for(j in 1:length(surrog)){
      n.surr<-nrow(surrog[[j]])  #<<---------------
      surrog.tmp<-list(surrog[[j]][1:(n.surr/2),], surrog[[j]][(n.surr/2+1):n.surr,])  #<<---------
      for(i.s in 1:2){    # <<--------------
        wts.surr<-warray(surrog.tmp[[i.s]],times=1:ncol(indata[[1]]),scale.min=2,scale.max.input=NULL,sigma=1.05,f0=1)$wave.array #take wavelet transform of surrogate data  <<<------------
        
        covs.list[[s[i.s]-1]]<-apply(wts.surr, 3, c) #replace covariate with surrogate data  <<---------
      }
      
      X1s<-wmodel.test(trans0f = resp,transnf = covs.list,n=n,freqs = 1:ncol(resp)) #make coefficients
      model<-vector("list", ncol(X1s))
      for(i in 1:ncol(X1s)){
        model[[i]]<-t(t(covs.list[[i]])*X1s[,i])
      }
      model2<-Reduce('+',model) #calculate predicted wavelet transform
      model.resids<-resp-model2 #calculate model residuals
      d2f.surr[,j]<-colMeans(model.resids*Conj(model.resids),na.rm=T) #calculate squared residuals by frequency to assess model fit; smaller is better
    }
    d2f.surr<-d2f.surr[is.finite(rowMeans(d2f.surr)),]
    pval<-c()
    for(i in 1:nrow(d2f.surr)){
      pval[i]<-sum(Re(d2f[i])>Re(d2f.surr[i,]),na.rm=T)/n.surrog
    }
    
    if(is.null(nrow(tsranges))){tsranges=matrix(tsranges, nrow=1, ncol=2, byrow=T)}
    sigtest<-surrogtest(Mod(d2f),Mod(t(d2f.surr)),timescales=timescales,tsranges=tsranges)
    
    return(list(d2f=d2f,d2f.surrog=d2f.surr,pval=pval,coefs=X1,timescales=timescales, sigtest=sigtest,pred.wt=model1))
  }
  else{
    return(list(d2f=d2f,coefs=X1,timescales=timescales,pred.wt=model1))
  }
}



wmrsig.twos<-function(indata, s, r, n, n.surrog=100, surr.test=FALSE, surr.type="fft", tsranges=c(0,Inf)){
  d1<-sapply(X=indata,FUN=function(x){return(dim(x)[1])}) 
  d2<-sapply(X=indata,FUN=function(x){return(dim(x)[2])}) 
  if (!all(d1[2:length(d1)]==d1[1]) || !all(d2[2:length(d2)]==d2[1]))
  {
    stop("Error: all matrices must be same dimensions")
  }
  
  #flatten input data
  wt.covs<-list()
  for(i in 1:length(indata)){
    wt.covs[[i]]<-warray(indata[[i]],times=1:ncol(indata[[i]]),scale.min=2,scale.max.input=NULL,sigma=1.05,f0=1)$wave.array
  }
  
  timescales<-wt(indata[[i]][1,],times=1:ncol(indata[[i]]))$timescales
  
  covs.wide<-list()
  for(j in 1:length(wt.covs)){
    xx<-c()
    for(i in 1:dim(wt.covs[[1]])[1]){
      xx<-cbind(xx,t(wt.covs[[j]][i,,]))
    }
    covs.wide[[j]]<-xx
  }
  
  resp<-covs.wide[[r]] #define response
  covs.list<-covs.wide[-r] #remove response from covariate list, thus defining covariates
  
  #make coefficients
  X1<-wmodel.test(trans0f = t(resp),transnf = lapply(covs.list,t),n=n,freqs = 1:nrow(resp)) #calculate coefficients
  coefs<-list()
  for(i in 1:ncol(X1)){
    coefs[[i]]<-matrix(X1[,i],nrow(resp),ncol(resp)) # reshape into matrices of coefficients
  }
  model<-list()
  for(i in 1:length(coefs)){
    model[[i]]<-covs.list[[i]]*coefs[[i]] #multiple coefficients by the data to get predicted values
  }
  model1<-Reduce('+',model) #calculate predicted wavelet transform by summing predicted values
  model.resids<-resp-model1 #calculate model residuals
  d2f<-colMeans(t(model.resids)*Conj(t(model.resids)),na.rm=T) #calculate squared residuals by frequency to assess model fit; smaller is better
  d2f<-d2f[is.finite(d2f)]
  #To conduct surrogate testing
  if(surr.test){ 
    if(surr.type=="fft"){surrog<-ffts(rbind(indata[[s[1]]],indata[[s[2]]]), nsurrogs=n.surrog)} #make surrogates  <<---------
    if(surr.type=="aaft"){surrog<-aafts(rbind(indata[[s[1]]],indata[[s[2]]]), nsurrogs=n.surrog)} # <<-------------
    
    d2f.surr<-matrix(ncol=n.surrog,nrow=nrow(resp))
    for(j in 1:length(surrog)){
      n.surr<-nrow(surrog[[j]])  #<<---------------
      surrog.tmp<-list(surrog[[j]][1:(n.surr/2),], surrog[[j]][(n.surr/2+1):n.surr,])  #<<---------
      for(i.s in 1:2){    # <<--------------
        wts.surr<-warray(surrog.tmp[[i.s]],times=1:ncol(indata[[1]]),scale.min=2,scale.max.input=NULL,sigma=1.05,f0=1)$wave.array #take wavelet transform of surrogate data  <<<------------
        wts.surr.wide<-c()
        for(i in 1:dim(wts.surr)[1]){
          wts.surr.wide<-cbind(wts.surr.wide,t(wts.surr[i,,])) #convert surrogate data to a freq x (location x time)
        }
        covs.list[[s[i.s]-1]]<-wts.surr.wide #replace covariate with surrogate data  <<---------
      }
      
      X1s<-wmodel.test(trans0f = t(resp),transnf = lapply(covs.list,t),n=n,freqs = 1:nrow(resp)) #make coefficients
      coefs<-list()
      for(i in 1:ncol(X1s)){
        coefs[[i]]<-matrix(X1s[,i],nrow(resp),ncol(resp)) # reshape into matrices of coefficients
      }
      model<-list()
      for(i in 1:length(coefs)){
        model[[i]]<-covs.list[[i]]*coefs[[i]] #multiple coefficients by the surrogate data to get predicted values
      }
      model2<-Reduce('+',model) #calculate predicted wavelet transform
      model.resids<-resp-model2 #calculate model residuals
      d2f.surr[,j]<-colMeans(t(model.resids)*Conj(t(model.resids)),na.rm=T) #calculate squared residuals by frequency to assess model fit; smaller is better
    }
    d2f.surr<-d2f.surr[is.finite(rowMeans(d2f.surr)),]
    pval<-c()
    for(i in 1:nrow(d2f.surr)){
      pval[i]<-sum(Re(d2f[i])>Re(d2f.surr[i,]),na.rm=T)/n.surrog
    }
    
    if(is.null(nrow(tsranges))){tsranges=matrix(tsranges, nrow=1, ncol=2, byrow=T)}
    sigtest<-surrogtest(Mod(d2f),Mod(t(d2f.surr)),timescales=timescales,tsranges=tsranges)
    
    return(list(d2f=d2f,d2f.surrog=d2f.surr,pval=pval,coefs=X1,timescales=timescales, sigtest=sigtest,pred.wt=model1))
  }
  else{
    return(list(d2f=d2f,coefs=X1,timescales=timescales,pred.wt=model1))
  }
}

