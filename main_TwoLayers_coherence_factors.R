rm(list=ls())
library(Reumannplatz)

ds.two<-readRDS("Data_TwoLayers_clean.RDS")
i.season<-2 #spring
ds.two<-ds.two[[i.season]]


######## sig test: shallow vs. deep #########
dists<-readRDS("res_distance_from_shore.RDS")
i1<-which(dists<240) #near-shore
i2<-which(dists>=240) #off-shore
i.dist<-list(i1,i2)

tsranges3<-matrix(c(2, 3, 3, 4, 4, Inf, 2,4, 4,Inf, 2,Inf), ncol=2, byrow=TRUE)
a<-1
ans<-vector("list", 4)
names(ans)<-c("shallow.near","shallow.off","deep.near","deep.off")
for(i in 1:2){  # shallow or deep
  for(j in 1:2){ # near or off
    x<-ds.two[[i]][[1]][i.dist[[j]],]  #T
    y<-ds.two[[i]][[3]][i.dist[[j]],]  #Chla
    sig.test<-cohtestfast(y,x, nsurrogs = 10000, min.scale = 2, max.scale = NULL, tsranges = tsranges3)
    sig.test<-sig.test[-(1:2)]
    sig.test$p.timescale<-1-sig.test$emp.rank/10000
    sig.test$Arg<-Arg(sig.test$emp.coh)
    
    x1<-ds.two[[i]][[2]][i.dist[[j]],]  #N
    sig.test1<-cohtestfast(y,x1, nsurrogs = 10000, min.scale = 2, max.scale = NULL, tsranges = tsranges3)
    sig.test1<-sig.test1[-(1:2)]
    sig.test1$p.timescale<-1-sig.test1$emp.rank/10000
    sig.test1$Arg<-Arg(sig.test1$emp.coh)
    
    ans[[a]]<-list(temp=sig.test, no3=sig.test1)
    a<-a+1
  }
}
saveRDS(ans,"res_coherence_factors_clean.RDS")


############ plot test: phase angle vs. freq ###################
ans<-readRDS("res_coherence_factors_clean.RDS")
### plot
tiff("Fig_test_phase_vs_freq_factors.tif", width=6, height=10, units="in",res=300,compression = "lzw")
op<-par(mfrow=c(4,2), oma=c(3,5,2,3), mar=c(2,2,1,1),mgp=c(2,0.5,0))
a<-1
for(i in 1:4){
  X<-ans[[i]]
  for(j in 1:2){
    x<-X[[j]]
    plot(1/x$timescales, x$Arg, ylim=c(-pi, pi), pch=20, xlab=NA, ylab=NA,cex=1.5)
    
    par(usr=c(par("usr")[1:2], c(0,1)))
    lines(1/x$timescales, x$p.timescale,  col="red")
    lines(range(1/x$timescales), c(0.05,0.05), lty="dashed", col="red")
    axis(4,col.axis="red",col.lab="red")
    
    mtext(paste0("(",letters[a],")"), side=3, line=-2.2, adj=0.05, cex=0.9)
    mtext(paste0("p(low-freq): ",round(x$pvals[5],3)), side=3, line=-3.2, adj=0.95,cex=0.8)
    mtext(paste0("p(high-freq): ",round(x$pvals[4],3)), side=3, line=-4.5, adj=0.95,cex=0.8)
    a<-a+1
  }
}
par(fig = c(0, 1, 0, 1), oma=c(1,3,0,0), mar = c(3, 3, 0, 0), new = TRUE)
plot(NA, xlim=c(0,1),ylim=c(0,1), xlab="frequency",ylab="phase angle",
     type = "n", bty = "n", xaxt = "n", yaxt = "n", cex.lab=1.2, font.lab=1)
mtext("Temperature", side=3, line=-2, adj=0.2, cex=0.9)
mtext("Nitrate", side=3, line=-2, adj=0.7, cex=0.9)
mtext("shallow.near", side=2, line=4, adj=0.9, cex=0.9)
mtext("shallow.off", side=2, line=4, adj=0.6, cex=0.9)
mtext("deep.near", side=2, line=4, adj=0.36, cex=0.9)
mtext("deep.off", side=2, line=4, adj=0.09, cex=0.9)
mtext("p-values", side=4, line=-2, cex=0.9, col="red")
par(op)
dev.off()




########### Fig4: rose diagram (chla vs. factors) ##############
ans<-readRDS("res_coherence_factors_clean.RDS")
library(circular)
tiff("Fig4_rose_diagram_factors.tif", width=5, height=5, units="in",res=300,compression = "lzw")
op<-par(mfrow=c(2,2), oma=c(1,1,2,1), mar=c(2,2,1,1),mgp=c(2,0.5,0))
a<-1
for(f in 1:3){
  for(j in 1:2){
    x<-ans[[a]]
    if(a==1){
      rose.phase<-Arg(x$emp.coh[x$timescales<4])  #discard low frequency for chla vs chla nearshore
    }else{rose.phase<-Arg(x$emp.coh)}
    
    rose.diag(circular(rose.phase), bins=16, col="grey")
    mtext(paste0("(",letters[a],")"), side=3, line=-1.2, adj=0.05, cex=0.9)
    a<-a+1
  }
}
par(fig = c(0, 1, 0, 1), oma=c(1,1,0,0), mar = c(3, 3, 0, 0), new = TRUE)
plot(NA, xlim=c(0,1),ylim=c(0,1), xlab=NA,ylab=NA,
     type = "n", bty = "n", xaxt = "n", yaxt = "n", cex.lab=1.2, font.lab=2)
mtext("near-shore", side=3, line=-2, adj=0.1, cex=0.9)
mtext("off-shore", side=3, line=-2, adj=0.78, cex=0.9)
mtext("Chl-a", side=2, line=1, adj=0.82, cex=0.9)
mtext("temperature", side=2, line=1, adj=0.45, cex=0.9)
mtext("nitrate", side=2, line=1, adj=0.08, cex=0.9)
par(op)
dev.off()




############### Fig4: rose diagram (Chla vs. T or NO3) #############
ans<-readRDS("res_coherence_factors_clean.RDS")
library(circular)
tiff("Fig4_rose_diagram_factors.tif", width=5, height=5, units="in",res=300,compression = "lzw")
op<-par(mfrow=c(2,2), oma=c(1,1,2,1), mar=c(2,1,1,0),mgp=c(2,1,0))
a<-1
for(i in 1:4){
  X<-ans[[a]]
  x1<-X[[1]] #T
  x2<-X[[2]] #N
  if(a==1){
    rose.phase1<-Arg(x1$emp.coh[x1$timescales<4])  #discard low frequency for chla vs chla nearshore
    rose.phase2<-Arg(x2$emp.coh[x2$timescales<4])
  }else{rose.phase1<-Arg(x1$emp.coh); rose.phase2<-Arg(x2$emp.coh)}
  
  rose.diag(circular(rose.phase1), bins=16, col="red", cex=0.9)
  rose.diag(circular(rose.phase2), bins=16, col="blue",cex=0.9, add=T)
  mtext(paste0("(",letters[a],")"), side=3, line=-1.2, adj=0.05, cex=0.9)
  a<-a+1
}
par(fig = c(0, 1, 0, 1), oma=c(1,1,0,0), mar = c(0, 0, 0, 0), new = TRUE)
plot(NA, xlim=c(0,1),ylim=c(0,1), xlab=NA,ylab=NA,
     type = "n", bty = "n", xaxt = "n", yaxt = "n", cex.lab=1.2, font.lab=2)
mtext("near-shore", side=3, line=-2, adj=0.22, cex=0.9)
mtext("off-shore", side=3, line=-2, adj=0.79, cex=0.9)
mtext("shallow", side=2, line=-1, adj=0.73, cex=0.9)
mtext("deep", side=2, line=-1, adj=0.25, cex=0.9)
legend(0.38,0.08, c("Chl-a vs. T", "Chl-a vs. NO3"), pch=15, col=c("red","blue"), cex=0.9)
par(op)
dev.off()
