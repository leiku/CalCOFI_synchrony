rm(list=ls())
library(Reumannplatz)
library(ecodist)

######### get data ready ############


ds.two<-readRDS("Data_TwoLayers_clean.RDS")
i.season<-2 #spring
ds.two<-ds.two[[i.season]]

dists<-readRDS("res_distance_from_shore.RDS")
i1<-which(dists<240) #near-shore
i2<-which(dists>=240) #off-shore

#offshore - all subsequent analyses are far-shore
C.deep<-ds.two[[2]][[3]][i2,]   #Chla.deep
C.shal<-ds.two[[1]][[3]][i2,]   #Chla.shallow
T.deep<-ds.two[[2]][[1]][i2,]   #temperature.deep
N.deep<-ds.two[[2]][[2]][i2,]   #nutrite.deep

X<-list(C.deep=C.deep, C.shal=C.shal, T.deep=T.deep, N.deep=N.deep)

######## wavelet multiple linear regression ###########
source('Fn_wmrsig.test.R')
tsranges2<-matrix(c(2,Inf, 2,4, 4,Inf), 3,2, byrow=T)  # all, short (2-4), long (4-12)
ans.sC<-wmrsig(X, r=1, n=3, s=2, n.surrog = 10000, surr.test = T, surr.type="fft", tsranges = tsranges2) #surrogate C
ans.sTN<-wmrsig.twos(X, r=1, n=3, s=c(3,4), n.surrog = 100, surr.test = T, surr.type="fft", tsranges = tsranges2) #surrogate T&N
ans.oC<-wmrsig(X[1:2], r=1, n=1, s=2, n.surrog = 10000, surr.test = T, surr.type="fft", tsranges = tsranges2)  #only contain C
ans.oTN<-wmrsig.twos(X[c(1,3,4)], r=1, n=2, s=c(2,3), n.surrog = 10000, surr.test = T, surr.type="fft", tsranges = tsranges2)  #only contain T&N
ans<-list(sC=ans.sC, sTN=ans.sTN, oC=ans.oC, oTN=ans.oTN)
saveRDS(ans, "res_wmrsig_three.RDS")

#plot
ans<-readRDS("res_wmrsig_three.RDS")

tiff("Fig_wmrsig_three.tif", width=5, height=8, units="in",res=300,compression = "lzw")
op<-par(mfrow=c(2,1), oma=c(1,1,1,1), mar=c(3,3,1,1),mgp=c(2,0.5,0))
names.panel<-c("both","solely")
for(i in 1:2){
  ans1<-ans[[2*i-1]]  #sC
  ans2<-ans[[2*i]]  #sTN
  plot(ans1$timescales, ans1$pval, typ="l", col="blue", ylim=c(0,1),xlab="timescale", ylab="p-value", cex.lab=1.3)
  lines(ans2$timescales, ans2$pval, col="red")
  lines(c(0,20),c(0.05, 0.05), lty="dashed")
  mtext(paste0("(",letters[i],") ", names.panel[i]), side=3, adj=0.02, line=-1.2, cex=0.9)
  if(i==2){legend("topright", legend=c("surr Chl.shallow", "surr T.deep"), lty="solid", col=c("blue","red"), horiz=F, cex=1)}
}
par(op)
dev.off()

################# matrix regression model selection ##############
#test the shading hypothesis by comparing the matrix model Cdeep~Cshal+Tdeep+Ndeep with the matrix models
#Cdeep~Tdeep+Ndeep and Cdeep~Cshal using model selection

Y<-lapply(X, synmat, method="pearson")   #list of correlation matrix

modelnames<-list(2, c(3,4), c(2,3,4))
ans3<-lno.weights(mats=Y[1:4], model.names=modelnames, n=3, nrand=1000, maxruns=1000) 
ans4<-summed.weights(varnames=names(Y[1:4]), weights=ans3)
#D: There is possibly a problem with summed.weights, definitely a problem with the example in summed.weights, possibly a problem with the way it is called here. Lei pls fix.

ans.mrm<-list(weights=ans3, importance=ans4)
saveRDS(ans.mrm, "res_mrm_three.RDS")

ans.mrm<-readRDS("res_mrm_three.RDS")

#test Cdeep~Tdeep+Ndeep vs Cdeep~Tdeep and Cdeep~Ndeep, likewise Cshal~Tshal+Nshal vs Cshal~Tshal and Cshal~Nshal
#D: Lei to do
ans3<-lno.weights(mats=Y[c(1,3,4)], n=3, nrand=1000, maxruns=1000) 
#D: add a call to summed.weights once fixed
#D: need to get the variables for shallow up above, then similar call for shallow

################ matrix regression tests ##############
Y<-lapply(X, synmat, method="pearson")   #list of correlation matrix
model1<-matregtest(resp=Y[[1]],preds=Y[2:4],drop=2:3,numperm=10000)
model2<-matregtest(resp=Y[[1]],preds=Y[2:4],drop=1,numperm=10000)


############### TEST LIMITING FACTORS #########
######## wavelet multiple linear regression #######
tsranges2<-matrix(c(2,12, 2,4, 4,12), 3,2, byrow=T)  # all, short (2-4), long (4-12)
ans.sT.TN<-wmrsig(X[c(1,3,4)], r=1, n=2, s=2, n.surrog = 10000, surr.test=T, tsranges=tsranges2)  #surrogate C
ans.sN.TN<-wmrsig(X[c(1,3,4)], r=1, n=2, s=3, n.surrog = 10000, surr.test=T, tsranges=tsranges2)  #surrogate T
ans.oT<-wmrsig(X[c(1,3)], r=1, n=1, s=2, n.surrog = 10000, surr.test = T, tsranges = tsranges2)  #only contain T
ans.oN<-wmrsig(X[c(1,4)], r=1, n=1, s=2, n.surrog = 10000, surr.test = T, tsranges = tsranges2)  #only contain N

ans<-list(sT=ans.sT.TN, sN=ans.sN.TN, oT=ans.oT, oN=ans.oN)
saveRDS(ans, "res_wmrsig_two_TN.RDS")

tiff("Fig_wmrsig_two_TN.tif", width=8, height=8, units="in",res=300,compression = "lzw")
op<-par(mfcol=c(2,1), oma=c(1,1,1,1), mar=c(3,3,1,1),mgp=c(2,0.5,0))
names.panel<-c("both","solely")
ans<-readRDS('res_wmrsig_two_TN.RDS')
for(i in 1:2){
  ans1<-ans[[2*i-1]]  #sC
  ans2<-ans[[2*i]]  #sT
  plot(ans1$timescales, ans1$pval, typ="l", col="blue", ylim=c(0,1),xlab="timescale", ylab="p-value", cex.lab=1.3)
  lines(ans2$timescales, ans2$pval, col="red")
  lines(c(0,20),c(0.05, 0.05), lty="dashed")
  mtext(paste0("(",letters[i],") ", names.panel[i]), side=3, adj=0.02, line=-1.2, cex=0.9)
  if(i==2){legend("topright", legend=c("surr T.deep", "surr N.deep"), lty="solid", col=c("blue","red"), horiz=F, cex=1.2)}
}
par(op)
dev.off()

