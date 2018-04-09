
####### get shallow and deep data #########
ds<-readRDS("Data_all_seasons_interpolate.RDS")

# get two layers (75m-layer contributed to both shallow and deep)  
dif.depth<-c(5, 10, 10, 15, 22.5, 25, 25, 25, 25)
DS.two<-vector("list", 4)
names(DS.two)<-names(ds)
for(i.season in 1:4){
  X.top<-vector("list",3)
  names(X.top)<-names(ds[[i.season]])
  X.bottom<-X.top
  for(f in 1:3){
    x<-ds[[i.season]][[f]]
    x.top<-0
    for(i in 1:5){
      x.top<-x.top+x[[i]]*dif.depth[i]
    }
    x.top<-x.top+x[[6]]*dif.depth[6]/2
    x.bottom<-(x[[7]]+x[[8]]+x[[9]])*25+x[[6]]*dif.depth[6]/2
    x.top<-x.top/(sum(dif.depth[1:5])+dif.depth[6]/2)
    x.bottom<-x.bottom/(sum(dif.depth[7:9])+dif.depth[6]/2)
    
    index<-c(1, 2, 8, 14, 15, 22:25, 32:35, 43:47, 3:7, 9:13, 16:21, 26:31, 36:42, 48:55)
    X.top[[f]]<-x.top[index,]
    X.bottom[[f]]<-x.bottom[index,]
  }
  ds.two<-list(shallow=X.top, deep=X.bottom)
  DS.two[[i.season]]<-ds.two
}
saveRDS(DS.two,"Data_TwoLayers.RDS")


###### the detrend version ##########
for(i.season in 1:4){
  for(f in 1:3){
    tmp<-CleanData(DS.two[[i.season]][[1]][[f]], normalize=F, detrend=T, rescale=F, do.plot=F)
    DS.two[[i.season]][[1]][[f]]<-tmp$cleandat
    tmp<-CleanData(DS.two[[i.season]][[2]][[f]], normalize=F, detrend=T, rescale=F, do.plot=F)
    DS.two[[i.season]][[2]][[f]]<-tmp$cleandat
  }
}
saveRDS(DS.two,"Data_TwoLayers_detrend.RDS")


###### the clean version ###########
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