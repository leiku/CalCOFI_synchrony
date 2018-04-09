
######### get detrend data ############
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

