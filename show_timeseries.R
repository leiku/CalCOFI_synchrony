rm(list=ls())
library(Reumannplatz)
library(igraph)

###### fig.test: show all the timeseries, and the index of sites #########
ds.two<-readRDS("Data_TwoLayers.RDS")
i.season<-2 #spring
ds.two<-ds.two[[i.season]]

c.top<-ds.two[[1]][[3]]
c.bottom<-ds.two[[2]][[3]]

network_layout.site<-function(input.site){
  n.site<-length(input.site)
  site.coordinate<-read.csv(file="CalCOFIStaPos113.csv")
  Site1<-site.coordinate$Line*10000+site.coordinate$Station
  Site1[Site1==818046.9]<-820047
  tmp<-match(input.site,Site1)
  
  longitude=-site.coordinate$Dlongitude[tmp]
  latitude<-site.coordinate$Station.Dlatitude[tmp]
  layout.position<-cbind(longitude,latitude)
  return(layout.position)
}

positions<-network_layout.site(as.numeric(rownames(c.top)))
tiff("fig_test_locations.tif", width=7, height=7, units="in",res=300,compression = "lzw")
plot(positions[,1], positions[,2], pch=1, cex=3, xlab="longitude", ylab="latitude")
text(positions[,1], positions[,2], 1:55, cex=1)
dev.off()

pdf("fig_test_timeseries_original.pdf", width=12, height=17)
op<-par(mfrow=c(11,5), oma=c(3,3,3,3), mar=c(1.5,1.5,1.5,1.5),mgp=c(2.2,0.5,0))
for(ii in 1:55){
  plot(1:28+1983, c.top[ii,], xlab=NA, ylab=NA, pch=20, col="darkgreen", cex.lab=1.2, yaxt="n")
  lines(1:28+1983, c.top[ii,],  col="darkgreen")
  axis(2,col.axis="darkgreen",col.lab="darkgreen")
  
  par(usr=c(par("usr")[1:2], range(c.bottom[ii,])))
  points(1:28+1983, c.bottom[ii,], pch=20, col="purple")
  lines(1:28+1983, c.bottom[ii,],  col="purple")
  axis(4,col.axis="purple",col.lab="purple")
  
  mtext(paste0("(",ii,")"), side=3, line=-1.2, adj=0.05, cex=0.9, col="red")
}
par(fig = c(0, 1, 0, 1), oma=c(0.5,0.5,1,0), mar = c(2.5, 2.5, 0, 0), mgp=c(1,0,0), new = TRUE)
plot(NA, xlim=c(0,1),ylim=c(0,1), xlab="year",ylab="concentration of Chl-a",
     type = "n", bty = "n", xaxt = "n", yaxt = "n", cex.lab=1.2, font.lab=1.5)
legend("top", c("shallow","deep"), col=c("darkgreen","purple"),lty="solid",pch=20, horiz=T,cex=1.2)
dev.off()

## show timeseries of site 11 as well as 6 and 10
ds<-readRDS("Data_all_seasons_interpolate.RDS")
i.season<-2

x<-ds[[i.season]][[3]]
plot(NA, xlim=range(1:28+1983),ylim=c(0,1.7),xlab="year",ylab="Chl-a")
for(i in 1:9){
  lines(1:28+1983, x[[i]][11,],col=rainbow(9)[i])
}
legend("top", names(x), pch=20,col=rainbow(9), horiz=T)

pdf("fig_test_timeseries_original_site11.pdf", width=10, height=8)
op<-par(mfrow=c(3,3), oma=c(3,3,3,3), mar=c(1.5,1.5,1.5,1.5),mgp=c(2.2,0.5,0))
for(i in 1:9){
  plot(1:28+1983, x[[i]][11,], ylim=range(x[[i]][c(6,10,11),]), xlab=NA, ylab=NA, pch=20, col="black", cex.lab=1.2)
  lines(1:28+1983, x[[i]][11,],  col="black")
  
  points(1:28+1983, x[[i]][6,], pch=20, col="red")
  lines(1:28+1983, x[[i]][6,],  col="red")
  points(1:28+1983, x[[i]][10,], pch=20, col="blue")
  lines(1:28+1983, x[[i]][10,],  col="blue")
  
  mtext(paste0("(",i,") ",names(x)[i],"m"), side=3, line=-1.2, adj=0.05, cex=0.9, col="red")
}
par(fig = c(0, 1, 0, 1), oma=c(0.5,0.5,1,0), mar = c(2.5, 2.5, 0, 0), mgp=c(1,0,0), new = TRUE)
plot(NA, xlim=c(0,1),ylim=c(0,1), xlab="year",ylab="concentration of Chl-a",
     type = "n", bty = "n", xaxt = "n", yaxt = "n", cex.lab=1.2, font.lab=1.5)
legend("top", c("11","6","10"), col=c("black","red","blue"),lty="solid",pch=20, horiz=T,cex=1.2)
dev.off()




#################### show correlations in map ################
library(igraph)

##locations
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

plot.correlation.map <- function(fig.name, panel.names, RR, Locations){
  tiff(fig.name, width=9, height=8, units="in",res=300,compression = "lzw")
  #layout(matrix(c(1,2,3,3,4,4,5,6), nrow=2, ncol=4), widths=c(2,4,1,2), heights=c(1,1)); layout.show(6)
  layout(matrix(c(1,3,2,4,5,5), nrow=2, ncol=3), widths=c(4,4,1), heights=c(1,1)); layout.show(5)
  op<-par(oma=c(0.5,0.5,0.5,0.5), mar=c(3.5,3.2,0.5,1.2),mgp=c(2.2,0.5,0))
  
  pal1=colorRampPalette(c("yellow", "red"), space="rgb")
  pal2=colorRampPalette(c( "cyan", "blue"), space="rgb")
  cols1<-pal1(100)
  cols2<-pal2(100)
  #pie(rep(1, 100), col = cols1, labels=NA, border=NA)
  #pie(rep(1, 100), col = cols2, labels=NA, border=NA)
  
  for(i in 1:4){
    inds1<-which(RR[[i]][[1]]<0)
    inds2<-which(RR[[i]][[1]]>=0)
    colors<-rep(NA, nrow(Locations))
    colors[inds1]<-cols1[ceiling(-100*RR[[i]][[1]][inds1])]
    colors[inds2]<-cols2[ceiling(100*RR[[i]][[1]][inds2])]
    plot(Locations[,1], Locations[,2], col=colors, pch=20, xlab="longitude", ylab="latitude", cex=4, cex.lab=1.2)
    inds<-which(RR[[i]][[2]]>0.05)  #indicating insignificant correlations
    points(Locations[inds,1], Locations[inds,2], col="white", pch=20, cex=1)
    mtext(paste0("(",letters[i],") ", panel.names[i]), side=3, line=-1.5, adj=0.05)
  }
  
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
  
  par(op)
  dev.off()
}


#get correlations

ds.two<-readRDS("Data_TwoLayers_detrend.RDS")
i.season<-2 #spring
ds.two<-ds.two[[i.season]]

locations<-network_layout.site(row.names(ds.two[[1]][[1]]))

####################### C T ##########
RR<-vector("list",7)
names(RR)<-c("CsCd","TsTd","CsTs","CdTd","NsNd","CsNs","CdNd")
i1<-c(3, 1, 3, 3, 2, 3, 3) # T N C
j1<-c(1, 1, 1, 2, 1, 1 ,2) # shallow deep
i2<-c(3, 1, 1, 1, 2, 2, 2)
j2<-c(2, 2, 1, 2, 2, 1, 2)

for(a in 1:7){
  rr<-rep(NA, nrow(locations))
  names(rr)<-row.names(ds.two[[1]][[1]])
  pp<-rr
  for(i in 1:nrow(locations)){
    x<-ds.two[[j1[a]]][[i1[a]]][i,]
    y<-ds.two[[j2[a]]][[i2[a]]][i,]
    rr[i]<-cor(x,y)
    y1<-aafts(y, 1000)
    z1<-rep(NA,1000)
    for(j in 1:1000){
      z1[j]<-cor(x,as.vector(y1[[j]]))
    }
    if(rr[i]>0){
      pp[i]<-sum(rr[i]<z1)/1000
    }else{pp[i]<-sum(rr[i]>z1)/1000}
    
  }
  
  RR[[a]]<-list(r=rr,pvalue=pp)
  show(a)
}
saveRDS(RR, "res_correlation_significance_detrend.RDS")

plot.correlation.map("fig_cor_map_C_N1.tif",names(RR[c(1,5:7)]), RR[c(1,5:7)], locations)
plot.correlation.map("fig_cor_map_C_T1.tif",names(RR[1:4]), RR[1:4], locations)




################## timeseries for shallow near-shore chla, T, and N ###################
ds.two<-readRDS("Data_TwoLayers_detrend.RDS")
i.season<-2 #spring
ds.two<-ds.two[[i.season]]

dists<-readRDS("res_distance_from_shore.RDS")
i1<-which(dists<240) #near-shore
i2<-which(dists>=240) #off-shore

RR<-readRDS("res_correlation_significance_detrend.RDS")

name.layer<-c("shallow","deep")
name.factor<-c("temperature","nitrate")
cols<-list(shallow=c("darkgreen","purple"), deep=c("blue","red"))

RR1<-RR[c(3,4,6,7)]
a<-1
pdf("fig_test_timeseries_chla_TN_detrend.pdf", width=12, height=17)
for(jj in 1:2){ #factor
  for(j in 1:2){ #shallow or deep
    X<-ds.two[[j]] 
    op<-par(mfrow=c(7,4), oma=c(3,3,3,3), mar=c(1.5,1.5,1.5,1.5),mgp=c(2.2,0.5,0))
    for(i in 1:28){
      plot(1:28+1983, X[[3]][i1[i],], xlab=NA, ylab=NA, pch=20, col=cols[[j]][1], cex.lab=1.2, yaxt="n")
      lines(1:28+1983, X[[3]][i1[i],],  col=cols[[j]][1])
      axis(2,col.axis=cols[[j]][1],col.lab=cols[[j]][1])
      
      par(usr=c(par("usr")[1:2], range(X[[jj]][i1[i],])))
      points(1:28+1983, X[[jj]][i1[i],], pch=20, col=cols[[j]][2])
      lines(1:28+1983, X[[jj]][i1[i],],  col=cols[[j]][2])
      axis(4,col.axis=cols[[j]][2],col.lab=cols[[j]][2])
      
      mtext(paste0("(",i,") site ",i1[i]), side=3, line=-1.2, adj=0.05, cex=0.9, col="black")
      rr<-RR1[[a]][[1]][i1[i]]
      pp<-RR1[[a]][[2]][i1[i]]
      mtext(paste0("r=",round(rr,3)), side=3, line=-1.2, adj=0.95, cex=0.9, col="black")
      mtext(paste0("p=",round(pp,3)), side=3, line=-2.5, adj=0.95, cex=0.9, col="black")
    }
    par(fig = c(0, 1, 0, 1), oma=c(0.5,0.5,0.5,0.5), mar = c(2, 1, 1, 1), mgp=c(1,0,0), new = TRUE)
    plot(NA, xlim=c(0,1),ylim=c(0,1), xlab="year",ylab=NA,
         type = "n", bty = "n", xaxt = "n", yaxt = "n", cex.lab=1.2, font.lab=1.5)
    mtext("concentration of Chl-a", side=2, line=-0.8, adj=0.5, cex=1, col=cols[[j]][1])
    mtext(name.factor[jj], side=4, line=-0.8, adj=0.5, cex=1, col=cols[[j]][2])
    mtext(name.layer[j], side=3, line=-1, adj=0.5, cex=1, col="black")
    par(op)
    a<-a+1
  }
}
dev.off()



name.factor<-c("temperature","nitrate","chl-a")
cols<-list(Temp=c("purple","red"), N=c("purple","red"), chla=c("darkgreen","blue"))

RR2<-RR[c(2,5,1)]
pdf("fig_test_timeseries_chla_TN_shallow_vs_deep_detrend.pdf", width=12, height=17)
for(j in 1:3){ #T N chla
  op<-par(mfrow=c(7,4), oma=c(3,3,3,3), mar=c(1.5,1.5,1.5,1.5),mgp=c(2.2,0.5,0))
  for(i in 1:28){
    x1<-ds.two[[1]][[j]][i1[i],]  #shallow
    x2<-ds.two[[2]][[j]][i1[i],]  #deep
    
    plot(1:28+1983, x1, xlab=NA, ylab=NA, pch=20, col=cols[[j]][1], cex.lab=1.2, yaxt="n")
    lines(1:28+1983, x1,  col=cols[[j]][1])
    axis(2,col.axis=cols[[j]][1],col.lab=cols[[j]][1])
    
    par(usr=c(par("usr")[1:2], range(x2)))
    points(1:28+1983, x2, pch=20, col=cols[[j]][2])
    lines(1:28+1983, x2,  col=cols[[j]][2])
    axis(4,col.axis=cols[[j]][2],col.lab=cols[[j]][2])
    
    mtext(paste0("(",i,") site ",i1[i]), side=3, line=-1.2, adj=0.05, cex=0.9, col="black")
    rr<-RR2[[j]][[1]][i1[i]]
    pp<-RR2[[j]][[2]][i1[i]]
    mtext(paste0("r=",round(rr,3)), side=3, line=-1.2, adj=0.95, cex=0.9, col="black")
    mtext(paste0("p=",round(pp,3)), side=3, line=-2.5, adj=0.95, cex=0.9, col="black")
  }
  par(fig = c(0, 1, 0, 1), oma=c(0.5,0.5,0.5,0.5), mar = c(2, 1, 1, 1), mgp=c(1,0,0), new = TRUE)
  plot(NA, xlim=c(0,1),ylim=c(0,1), xlab="year",ylab=NA,
       type = "n", bty = "n", xaxt = "n", yaxt = "n", cex.lab=1.2, font.lab=1.5)
  mtext("shallow", side=2, line=-0.8, adj=0.5, cex=1, col=cols[[j]][1])
  mtext("deep", side=4, line=-0.8, adj=0.5, cex=1, col=cols[[j]][2])
  mtext(name.factor[j], side=3, line=-1, adj=0.5, cex=1, col="black")
  par(op)
}
dev.off()




