rm(list=ls())
library(Reumannplatz)

######### get data ready ############
DS.two<-readRDS("Data_TwoLayers.RDS")
for(i.season in 1:4){
  for(f in 1:3){
    tmp<-CleanData(DS.two[[i.season]][[1]][[f]], normalize=T, each.ts=T, rescale=T, do.plot=F)
    DS.two[[i.season]][[1]][[f]]<-tmp$cleandat
    tmp<-CleanData(DS.two[[i.season]][[2]][[f]], normalize=T, detrend=T, rescale=T, do.plot=F)
    DS.two[[i.season]][[2]][[f]]<-tmp$cleandat
  }
}
saveRDS(DS.two,"Data_TwoLayers_clean.RDS")


########### Wavelet multiple linear regression ###########
ds.two<-readRDS("Data_TwoLayers_clean.RDS")
i.season<-2 #spring
ds.two<-ds.two[[i.season]]

dists<-readRDS("res_distance_from_shore.RDS")
i1<-which(dists<240) #near-shore
i2<-which(dists>=240) #off-shore

#offshore
C.deep<-ds.two[[2]][[3]][i2,]   #Chla.deep
C.shal<-ds.two[[1]][[3]][i2,]   #Chla.shallow
T.deep<-ds.two[[2]][[1]][i2,]   #temperature.deep
N.deep<-ds.two[[2]][[2]][i2,]   #nutrite.deep

X<-list(C.deep=C.deep, C.shal=C.shal, T.deep=T.deep, N.deep=N.deep)

###########################################################
################## three predictors ########################
#########################################################

######## wavelet multiple linear regression ###########
source('Fn_wmrsig.test.R')
tsranges2<-matrix(c(2,12, 2,4, 4,12), 3,2, byrow=T)  # all, short (2-4), long (4-12)
ans.sC<-wmrsig(X, r=1, n=3, s=2, n.surrog = 10000, surr.test = T, surr.type="fft", tsranges = tsranges2) #surrogate C
ans.sTN<-wmrsig.twos(X, r=1, n=3, s=c(3,4), n.surrog = 10000, surr.test = T, surr.type="fft", tsranges = tsranges2) #surrogate T&N
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


################# matrix regression ##############
Y<-lapply(X, synmat, method="pearson")   #list of correlation matrix

modelnames<-list(2, c(3,4), c(2,3,4))
modelnames<-lapply(modelnames, as.integer)
ans3<-lno.weights(mats=Y[1:4], model.names=modelnames, n=3, nrand=1000, maxruns=1000) 
ans4<-summed.weights(varnames=names(Y[1:4]), weights=ans3)

ans.mrm<-list(weights=ans3, importance=ans4)
saveRDS(ans.mrm, "res_mrm_three.RDS")

ans.mrm<-readRDS("res_mrm_three.RDS")

###########################################################
################## two predictors ########################
#########################################################

########### Wavelet multiple linear regression ###########
tsranges2<-matrix(c(2,12, 2,4, 4,12), 3,2, byrow=T)  # all, short (2-4), long (4-12)
ans.sC.CT<-wmrsig(X[1:3], r=1, n=2, s=2, n.surrog = 10000, surr.test=T, tsranges=tsranges2)  #surrogate C
ans.sT.CT<-wmrsig(X[1:3], r=1, n=2, s=3, n.surrog = 10000, surr.test=T, tsranges=tsranges2)  #surrogate T
#seperately
ans.oC<-wmrsig(X[1:2], r=1, n=1, s=2, n.surrog = 10000, surr.test = T, tsranges = tsranges2)  #only contain C
ans.oT<-wmrsig(X[c(1,3)], r=1, n=1, s=2, n.surrog = 10000, surr.test = T, tsranges = tsranges2)  #only contain T
ans<-list(sC.CT=ans.sC.CT, sT.CT=ans.sT.CT, oC=ans.oC, oT=ans.oT)
saveRDS(ans, "res_wmrsig_Chla.shallow_T.deep.RDS")

ans.sC.CN<-wmrsig(X[c(1,2,4)], r=1, n=2, s=2, n.surrog = 10000, surr.test=T, tsranges=tsranges2)  #surrogate C
ans.sN.CN<-wmrsig(X[c(1,2,4)], r=1, n=2, s=3, n.surrog = 10000, surr.test=T, tsranges=tsranges2)  #surrogate N
#seperately
ans.oC<-wmrsig(X[1:2], r=1, n=1, s=2, n.surrog = 10000, surr.test = T, tsranges = tsranges2)  #only contain C
ans.oN<-wmrsig(X[c(1,4)], r=1, n=1, s=2, n.surrog = 10000, surr.test = T, tsranges = tsranges2)  #only contain N
ans<-list(sC.CN=ans.sC.CN, sN.CN=ans.sN.CN, oC=ans.oC, oN=ans.oN)
saveRDS(ans, "res_wmrsig_Chla.shallow_N.deep.RDS")

#synexp.Csurf<-modelsyncexp(X[[1]],X[[2]], 1:28, tsrange=c(2,4))
#synexp.Tdeep<-modelsyncexp(X[[1]],X[[3]], 1:28, tsrange=c(2,4))

#plot

tiff("Fig_wmrsig_two.tif", width=8, height=8, units="in",res=300,compression = "lzw")
op<-par(mfcol=c(2,2), oma=c(1,1,1,1), mar=c(3,3,1,1),mgp=c(2,0.5,0))
names.panel<-c("both","solely")
ans<-readRDS('res_wmrsig_Chla.shallow_T.deep.RDS')
for(i in 1:2){
  ans1<-ans[[2*i-1]]  #sC
  ans2<-ans[[2*i]]  #sT
  plot(ans1$timescales, ans1$pval, typ="l", col="blue", ylim=c(0,1),xlab="timescale", ylab="p-value", cex.lab=1.3)
  lines(ans2$timescales, ans2$pval, col="red")
  lines(c(0,20),c(0.05, 0.05), lty="dashed")
  mtext(paste0("(",letters[i*2-1],") ", names.panel[i]), side=3, adj=0.02, line=-1.2, cex=0.9)
  if(i==2){legend("topright", legend=c("surr Chl.shallow", "surr T.deep"), lty="solid", col=c("blue","red"), horiz=F, cex=1.2)}
}
ans<-readRDS('res_wmrsig_Chla.shallow_N.deep.RDS')
for(i in 1:2){
  ans1<-ans[[2*i-1]]  #sC
  ans2<-ans[[2*i]]  #sN
  plot(ans1$timescales, ans1$pval, typ="l", col="blue", ylim=c(0,1),xlab="timescale", ylab="p-value", cex.lab=1.3)
  lines(ans2$timescales, ans2$pval, col="red")
  lines(c(0,20),c(0.05, 0.05), lty="dashed")
  mtext(paste0("(",letters[i*2],") ", names.panel[i]), side=3, adj=0.02, line=-1.2, cex=0.9)
  if(i==2){legend("topright", legend=c("surr Chl.shallow", "surr N.deep"), lty="solid", col=c("blue","red"), horiz=F, cex=1.2)}
}
par(op)
dev.off()


################# matrix regression ##############

#generate one leave-n-out score for a single model
#ans1<-lno.score(mats=Y[1:3], resp=1, pred=2:3, n=3, maxruns=1000)

#Rank all models for a given response by lno score
#ans2<-lno.ranking(mats=Y[1:3], n=3, maxruns=1000, rank.mod=T)

#Generate lno weights (frequency of occurrence of a given top model) using resampling procedures
#ans3<-lno.weights(mats=Y[1:3], model.names=modelnames, n=3, nrand=1000, maxruns=1000) 

#calculate importance
#ans4<-summed.weights(varnames=names(Y[1:3]), weights=ans3)


