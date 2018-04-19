rm(list=ls())
library(Reumannplatz)
source('D:/Git/CalCOFI_synchrony/Fn_wmodel.test.R')

time=100
p1 = 10
p2 = 7

sig1<-rep(sin(seq(0,2*pi,length.out=p1)[-p1]),length.out=time)+rnorm(time,0,0.1)
sig2<-rep(sin(seq(0,2*pi,length.out=p2)[-p2]),length.out=time)+rnorm(time,0,0.1)
sig1<-sig1-mean(sig1)
sig2<-sig2-mean(sig2)
sig3<-sig1+sig2+rnorm(time,0,0.1) #sig3 is the sum of artificial signals sig1 and sig2
sig3<-sig3-mean(sig3)

op<-par(mfrow=c(3,1), oma=c(1,1,1,1), mar=c(3,3,1,1),mgp=c(2.2,0.5,0)) 
plot(sig1, type="l")
plot(sig2, type="l")
plot(sig3, type="l")
par(op)

trans0f<-wt(sig3, 1:time)$wave
transnf<-list(wt(sig1, 1:time)$wave,
              wt(sig2, 1:time)$wave)
             
X1<-wmodel(trans0f,transnf,n=2,freqs=1:ncol(trans0f)) #compute coefficients
matX1_1<-matrix(X1[,1],ncol(trans0f),time) # reshape into matrices of coefficients
matX1_2<-matrix(X1[,2],ncol(trans0f),time)
model<-transnf[[1]]*t(matX1_1) + transnf[[2]]*t(matX1_2) #calculate predicted wavelet transform
model.resids<-trans0f-model #calculate model residuals
d1<-mean(t(model.resids)*Conj(t(model.resids)),na.rm=T) #calculate squared residuals to assess model fit; smaller is better
print(d1)

X2<-wmodel.test(trans0f,transnf,n=2,freqs=1:ncol(trans0f)) #compute coefficients
matX1_1<-matrix(X2[,1],ncol(trans0f),time) # reshape into matrices of coefficients
matX1_2<-matrix(X2[,2],ncol(trans0f),time)
model2<-transnf[[1]]*t(matX1_1) + transnf[[2]]*t(matX1_2) #calculate predicted wavelet transform
model.resids2<-trans0f-model2 #calculate model residuals
d2<-mean(t(model.resids2)*Conj(t(model.resids2)),na.rm=T) #calculate squared residuals to assess model fit; smaller is better
print(d2)



op<-par(mfrow=c(2,2), oma=c(1,1,1,1), mar=c(3,3,1,1),mgp=c(2.2,0.5,0))  #plot the output: these should be close but not identical.
image(Re(trans0f))
image(Re(model))
image(Re(trans0f))
image(Re(model2))
par(op)
