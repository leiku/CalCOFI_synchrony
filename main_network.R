rm(list=ls())
library(Reumannplatz)
library(igraph)

######### get clean data ############
ds<-readRDS("Data_all_seasons_interpolate.RDS")

ds.clean<-ds
for(s in 1:4){
  for(f in 1:3){  #factor 
    for(d in 1:9){ #depth
      tmp1<-ds[[s]][[f]][[d]]
      #y<-CleanData.mod(tmp1, each.ts=T)
      y<-CleanData(tmp1, normalize=F, detrend=T, rescale=F, do.plot=F )
      ds.clean[[s]][[f]][[d]]<-y$cleandat
    }
    show(f)
  }
}

#saveRDS(ds.clean, "Data_all_seasons_clean.RDS")   # interpolate, normalize (Box-Cox), detrend, rescale
saveRDS(ds.clean, "Data_all_seasons_detrend.RDS")    # detrend

######## synchrony network ###########

#get position
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

#plot
plot.cluster.Calcofi<-function(figure.name, members0, node.contribution, locations, name.depth){
  tiff(figure.name, 
       width=9, height=8, units="in",res=300,compression = "lzw")
  layout(matrix(c(1,4,4,7,2,5,5,8,3,6,6,9,10,10,11,11), nrow=4, ncol=4), widths=c(4,4,4,1), heights=c(4,2,2,4)); layout.show(10)
  op<-par(oma=c(3,3,1,1), mar=c(1,1,1,1),mgp=c(2,0.5,0))
  n.m<-length(unique(members0))
  if(n.m>3){stop("The number of modules is larger than 2. Please change the color setting")}
  inds1<-which(members0==1)
  inds2<-which(members0==2)
  
  pal1=colorRampPalette(c("yellow", "red"), space="rgb")
  pal2=colorRampPalette(c( "cyan", "blue"), space="rgb")
  cols1<-pal1(100)
  cols2<-pal2(100)
  #pie(rep(1, 100), col = cols1, labels=NA, border=NA)
  
  colors<-rep(NA, length(members0))
  colors[inds1]<-cols1[ceiling(100*node.contribution[inds1])+1]
  colors[inds2]<-cols2[ceiling(100*node.contribution[inds2])+1]
  if(n.m>2){inds3<-which(members0==3); cols3<-topo.colors(1000); colors[inds3]<-cols3[300-round(300*node.contribution[inds3])+1]}
  
  for(d in 1:9){
    a<-55*(d-1)+1; b<-55*d
    plot(locations[,1], locations[,2], col=colors[a:b], pch=20, xlab=NA, ylab=NA, cex=3)
    mtext(paste0("(", letters[d],")"), side=3, line=-1.2, adj=0.05, cex=0.8)
    mtext(paste0("(",name.depth[d]," m)"), side=1, line=-1.2, adj=0.95, cex=0.8)
  }
  #color bar
  op<-par(mar = c(2.5, 1, 2.5, 2))  
  plot(c(0.8,1), c(0,1), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
  axis(4, seq(0,1,len=6), las=1, labels=seq(0,1,len=6))
  scale<-length(cols1)-1
  for (i in 1:(length(cols2)-1)) {
    y = (i-1)/scale + 0
    rect(0,y,10,y+1/scale, col=cols1[i], border=NA)
  }
  mtext("cluster 1", side=3, line=0, cex=0.9)
  
  plot(c(0.8,1), c(0,1), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
  axis(4, seq(0,1,len=6), las=1, labels=seq(0,1,len=6))
  scale<-length(cols2)-1
  for (i in 1:(length(cols2)-1)) {
    y = (i-1)/scale + 0
    rect(0,y,10,y+1/scale, col=cols2[i], border=NA)
  }
  mtext("cluster 2", side=3, line=0, cex=0.9)
  
  par(fig = c(0, 1, 0, 1), oma=c(1,1,0,0), mar = c(3, 3, 0, 0), new = TRUE)
  plot(NA, xlim=c(0,1),ylim=c(0,1), xlab="Longitude",ylab="Latitude",
       type = "n", bty = "n", xaxt = "n", yaxt = "n", cex.lab=1.2, font.lab=2)
  par(op)
  dev.off()
}

#################### for original data
source('SynNetModulesSimple.mod.pearson.R')
source('Fn_newman_eigenvector_mod.R')
source('Fn_newman_modularity_mod.R')

ds<-readRDS("Data_all_seasons_detrend.RDS")    # detrend  (should use detrend data or original data?)
i.season<-2 #spring

ds1<-ds[[i.season]][[3]]
ds2<-ds1[[1]]
for(d in 2:9){
  ds2<-rbind(ds2, ds1[[d]])
}

name.season<-names(ds)
name.depth<-names(ds[[1]][[1]])[1:9]
Sites<-as.integer(row.names(ds[[1]][[1]][[1]]))
locations<-network_layout.site(Sites)

#plot 
#ans<-SynNetModules(ds2, 1:28, coords=NA, method="pearson", nsurrogs = NA, do.plot=F) # errors: the output synmat is a function
#decomp<-ModularityDecomp(ans$synmat, ans$modules$membership)   # errors: it is not for negative networks

ans<-SynNetModulesSimple.mod.pearson(ds2, 1:28)
members<-ans$modules$membership
decomp<-newman_modularity_mod(ans$synmat, members, decomp=T)
plot.cluster.Calcofi(paste0("Fig1_cluster_membership_negative_",name.season[i.season],"_pearson_detrend.tif"), 
                     members, decomp$Q.decomp.node.rescale, locations=locations, name.depth=name.depth)



