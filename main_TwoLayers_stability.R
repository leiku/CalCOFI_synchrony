
rm(list=ls())
ds.two<-readRDS("Data_TwoLayers.RDS")  #original or detrend?
i.season<-2 #spring
ds.two<-ds.two[[i.season]]

c.top<-ds.two[[1]][[3]]
c.bottom<-ds.two[[2]][[3]]

############# Fig5: CoVar #########
cv<-matrix(NA, nrow(c.top), 3)
row.names(cv)<-row.names(c.top)
colnames(cv)<-c("top","bottom","all")
for(i in 1:nrow(c.top)){
  cv[i,1]<-sd(c.top[i,])/mean(c.top[i,])
  cv[i,2]<-sd(c.bottom[i,])/mean(c.bottom[i,])
  cv[i,3]<-sd(c.top[i,]+c.bottom[i,])/mean(c.top[i,]+c.bottom[i,])
}

## surrogate
cv.rand<-matrix(NA, nrow(c.top), 1000)
for(i in 1:nrow(c.top)){
  x<-c.top[i,]
  y<-c.bottom[i,]
  for(j in 1:1000){
    y1<-sample(y, length(y), replace=F)
    cv.rand[i,j]<-sd(x+y1)/mean(x+y1)
  }
}
saveRDS(list(cv=cv, cv.rand=cv.rand), "res_CoVar_rand.RDS")

#### plot
CV<-readRDS("res_CoVar_rand.RDS")
cv<-CV[[1]]
cv.rand<-CV[[2]]
cv.rand.ci<-t(apply(cv.rand, 1, quantile, c(0.025, 0.5, 0.975)))

cols<-rep("purple", nrow(cv))
cols[cv[,3]<=cv.rand.ci[,1]]<-"forestgreen"
cols[cv[,3]>=cv.rand.ci[,3]]<-"red"

tiff("Fig5_cv_distance_shore_original.tif", width=8, height=5, units="in",res=300,compression = "lzw")
op<-par(oma=c(1,1,1,1), mar=c(3,3,1,1),mgp=c(2,0.5,0))

plot(dists, cv[,3], ylim=range(cbind(cv[,3], cv.rand.ci)), pch=16, col=cols, xlab="distance from shore (km)", ylab="coefficient of variation")
arrows(dists, cv.rand.ci[,1], dists, cv.rand.ci[,3], length=0.05, angle=90, code=3)
lines(c(240, 240),c(0,1.5), lty="dashed")
par(op)
dev.off()


