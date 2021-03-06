rm(list=ls())
library(Reumannplatz)
library(igraph)


###### Fig2: Correlations between shallow and deep Chla  ##################
ds.two<-readRDS("Data_TwoLayers_detrend.RDS")
i.season<-2 #spring
ds.two<-ds.two[[i.season]]

c.top<-ds.two[[1]][[3]]
c.bottom<-ds.two[[2]][[3]]
#get correlation between top and bottom values of Chla
rr<-rep(NA, nrow(c.top))
for(i in 1:nrow(c.top)){
  z<-cor.test(c.top[i,], c.bottom[i,])
  rr[i]<-z$estimate
}
inds1<-which(rr<0)
inds2<-which(rr>=0)

##distance from shore
network_layout.site<-function(input.site){
  n.site<-length(input.site)
  site.coordinate<-read.csv(file="D:/KU/Work/CalCOFI/Locations/CalCOFIStaPos113.csv")
  Site1<-site.coordinate$Line*10000+site.coordinate$Station
  Site1[Site1==818046.9]<-820047
  tmp<-match(input.site,Site1)
  
  longitude=-site.coordinate$Dlongitude[tmp]
  latitude<-site.coordinate$Station.Dlatitude[tmp]
  layout.position<-cbind(longitude,latitude)
  return(layout.position)
}
locations<-network_layout.site(row.names(c.top))

library(geosphere)
site.coastal<-c(934026.4, 900027.7, 868032.5, 833039.4, 817043.5, 800050.5)
loc.coastal<-network_layout.site(site.coastal)
Line1<-as.numeric(row.names(c.top))%/%1000*0.1
Line2<-site.coastal%/%1000*0.1
Line2<-c(Line2,76.7)
loc.coastal<-rbind(loc.coastal,c(-120.75,35.17))
dists<-rep(NA,nrow(locations))
for(i in 1:nrow(locations)){
  j<-which(abs(Line2-Line1[i])<0.5)
  dists[i]<-distm(locations[i,], loc.coastal[j,])
}
dists<-dists*0.001  # m -> km
saveRDS(dists,"res_distance_from_shore.RDS")

#### timeseries
tiff("Fig2_correlation_ts_detrend.tif", width=7, height=9, units="in",res=300,compression = "lzw")
#layout(matrix(c(1,2,3,3,4,4,5,6), nrow=2, ncol=4), widths=c(2,4,1,2), heights=c(1,1)); layout.show(6)
layout(matrix(c(1,3,5,2,3,5,2,4,5), nrow=3, ncol=3), widths=c(4,3,1), heights=c(2,5,3)); layout.show(5)
op<-par(oma=c(0.5,0.5,0.5,0.5), mar=c(3.5,3.2,0.5,1.2),mgp=c(2.2,0.5,0))

ii<-38
plot(1:28+1983, c.top[ii,], xlab=NA, ylab=NA, pch=20, col="darkgreen", cex.lab=1.2, yaxt="n")
lines(1:28+1983, c.top[ii,],  col="darkgreen")
axis(2,col.axis="darkgreen",col.lab="darkgreen")

par(usr=c(par("usr")[1:2], range(c.bottom[ii,])))
points(1:28+1983, c.bottom[ii,], pch=20, col="purple")
lines(1:28+1983, c.bottom[ii,],  col="purple")
axis(4,col.axis="purple",col.lab="purple")
mtext(paste0("(a)"), side=3, line=-1.2, adj=0.05, cex=1)

ii<-8
plot(1:28+1983, c.top[ii,], xlab=NA, ylab=NA, pch=20, col="darkgreen", cex.lab=1.2, yaxt="n")
lines(1:28+1983, c.top[ii,],  col="darkgreen")
axis(2,col.axis="darkgreen",col.lab="darkgreen")

par(usr=c(par("usr")[1:2], range(c.bottom[ii,])))
points(1:28+1983, c.bottom[ii,], pch=20, col="purple")
lines(1:28+1983, c.bottom[ii,],  col="purple")
axis(4,col.axis="purple",col.lab="purple")
legend("topright",c("shallow","deep"),col=c("darkgreen","purple"),lty="solid",pch=20, horiz=T,cex=1.2)
mtext(paste0("(b)"), side=3, line=-1.2, adj=0.05, cex=1)


#### map

#cols1<-heat.colors(300)
#cols2<-topo.colors(600)

pal1=colorRampPalette(c("yellow", "red"), space="rgb")
pal2=colorRampPalette(c( "cyan", "blue"), space="rgb")
cols1<-pal1(100)
cols2<-pal2(100)
#pie(rep(1, 100), col = cols1, labels=NA, border=NA)
#pie(rep(1, 100), col = cols2, labels=NA, border=NA)
colors<-rep(NA, nrow(locations))
colors[inds1]<-cols1[ceiling(-100*rr[inds1])]
colors[inds2]<-cols2[ceiling(100*rr[inds2])]

plot(locations[,1], locations[,2], col=colors, pch=20, xlab="longitude", ylab="latitude", cex=4, cex.lab=1.2)
points(locations[38,1], locations[38,2], col="darkgrey", cex=5, pch=1)
points(locations[8,1], locations[8,2], col="darkgrey", cex=5, pch=1)
text(locations[38,1], locations[38,2]+0.3, "a", col="black",cex=1.1)
text(locations[8,1], locations[8,2]+0.3, "b", col="black",cex=1.1)
mtext("(c)", side=3, line=-1.2, adj=0.05, cex=1)

#color bar
op<-par(mar = c(3.5, 0.5, 0.5, 3))  
plot(c(0.8,1), c(-1,1), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
axis(4, seq(-1,1,len=11), las=1, labels=seq(-1,1,len=11))
scale<-length(cols2)-1
for (i in 1:(length(cols2)-1)) {
  y = (i-1)/scale + 0
  rect(0,y,10,y+1/scale, col=cols2[i], border=NA)
}

scale<-length(cols1)-1
for (i in 1:(length(cols1)-1)) {
  y = (i-1)/scale - 1
  rect(0,y,10,y+1/scale, col=cols1[101-i], border=NA)
}


#### correlations vs distance from shore
op<-par(mar=c(3.5,3,0.5,0.5))
plot(dists, rr, xlab="distance from shore (km)", ylab="correlation", pch=16, cex.lab=1.2)
#points(dists[pp<0.05], rr[pp<0.05], pch=16)
lines(c(-100,800),c(0,0),lty="dashed")
lines(c(240, 240),c(-1,1), lty="dashed")
mtext("(d)", side=3, line=-1.2, adj=0.95)
text(100, -0.5, "near-shore", col="red",cex=1.2)
text(500, 0.5, "off-shore", col="red",cex=1.2)

par(op)
dev.off()

