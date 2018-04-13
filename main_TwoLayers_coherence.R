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
ans<-vector("list", 6)
names(ans)<-c("C.near","C.off","T.near","T.off","N.near","N.off")
for(i in c(3,1,2)){
  for(j in 1:2){
    x.top<-ds.two[[1]][[i]][i.dist[[j]],]
    x.bottom<-ds.two[[2]][[i]][i.dist[[j]],]
    sig.test<-cohtestfast(x.bottom, x.top, nsurrogs = 10000, min.scale = 2, max.scale = NULL, tsranges = tsranges3)
    sig.test<-sig.test[-(1:2)]
    sig.test$p.timescale<-1-sig.test$emp.rank/10000
    sig.test$Arg<-Arg(sig.test$emp.coh)
    ans[[a]]<-sig.test
    a<-a+1
  }
}
saveRDS(ans,"res_coherence_clean.RDS")


############ plot test: phase angle vs. freq ###################
ans<-readRDS("res_coherence_clean.RDS")
### plot
tiff("Fig_test_phase_vs_freq.tif", width=6, height=6, units="in",res=300,compression = "lzw")
op<-par(mfrow=c(3,2), oma=c(3,5,2,3), mar=c(2,2,1,1),mgp=c(2,0.5,0))
a<-1
for(f in 1:3){
  for(j in 1:2){
    x<-ans[[a]]
    plot(x$timescales, x$Arg, ylim=c(-pi, pi), pch=20, xlab=NA, ylab=NA,cex=1.5)

    par(usr=c(par("usr")[1:2], c(0,1)))
    lines(x$timescales, x$p.timescale,  col="red")
    lines(range(x$timescales), c(0.05,0.05), lty="dashed", col="red")
    axis(4,col.axis="red",col.lab="red")
    
    lines(c(4,4),c(-4,4), lty="dashed", col="darkgray")
    
    p<-as.character(round(x$pvals[4:5],3))
    p[p!=0]<-paste0("=",p[p!=0])
    p[p==0]<-"<0.001"

    mtext(paste0("(",letters[a],")"), side=3, line=-2.5, adj=0.05, cex=0.9)
    mtext(paste0("p",p[1]), side=3, line=-4, adj=0.05,cex=0.8)
    mtext(paste0("p",p[2]), side=3, line=-4, adj=0.95,cex=0.8)
    a<-a+1
  }
}
par(fig = c(0, 1, 0, 1), oma=c(1,3,0,0), mar = c(3, 3, 0, 0), new = TRUE)
plot(NA, xlim=c(0,1),ylim=c(0,1), xlab="frequency",ylab="phase angle",
     type = "n", bty = "n", xaxt = "n", yaxt = "n", cex.lab=1, font.lab=1)
mtext("near-shore", side=3, line=-2, adj=0.2, cex=0.9, col="blue")
mtext("off-shore", side=3, line=-2, adj=0.7, cex=0.9, col="blue")
mtext("Chl-a", side=2, line=4, adj=0.82, cex=0.9, col="blue")
mtext("temperature", side=2, line=4, adj=0.5, cex=0.9, col="blue")
mtext("nitrate", side=2, line=4, adj=0.12, cex=0.9, col="blue")
mtext("p-values", side=4, line=-2, cex=0.9, col="red")
par(op)
dev.off()




########### Fig3: rose diagram (shallow vs. deep) ##############
ans<-readRDS("res_coherence_clean.RDS")
library(circular)
tiff("Fig3_rose_diagram.tif", width=4, height=6, units="in",res=300,compression = "lzw")
op<-par(mfrow=c(3,2), oma=c(1,1,2,1), mar=c(2,2,1,1),mgp=c(2,0.5,0))
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


