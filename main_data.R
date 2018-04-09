## ERDDAP: CalCOFI NOAAHydros ###

rm(list=ls())

########### basic station ###############
name.site<-c(767000+c(100,90,80,70,60,55,51),
             800000+c(100,90,80,70,60,55,51),820047,
             833000+c(110,100,90,80,70,60,55,51,42),
             867000+c(110,100,90,80,70,60,55,50,45,40,35),
             900000+c(120,110,100,90,80,70,60,53,45,35,28),
             933000+c(120,110,100,90,80,70,60,55,50,45,40,35,28))
costal.site<-c(800050.5,817043.5,833039.4,868032.5,
               885030.1,900027.7,934026.4,917026.4,854035.8)

###########   bottle data #############
D0<-read.csv(file="org_bottle.csv")
D0$year<-as.numeric(substr(D0$time,1,4))
D0$month<-as.numeric(substr(D0$time,6,7))
D0$line<-as.numeric(substr(D0$sta_id,1,5))
D0$station<-as.numeric(substr(D0$sta_id,7,11))

D0<-subset(D0, line>70 & line<100 & station>20 & station<125)
D0$site <- D0$line*10000+D0$station   #use site to indicate line and station

D0<-subset(D0, year>=1984 & year<=2011)

D0$monthid <- (D0$year-1984)*12+D0$month

D0<-D0[,c(29,30,5,6,8,12,14:19)]
colnames(D0)<-c("site","monthid","depth","temperature",
                "salinity","o2","sio3","po4","no3","no2",
                "chla","phaeo")
std.dep<-c(0,10,20,30,50,75,100,125,150,200,250,300,400,500)
x<-cut(D0$r_depth,c(-1,std.dep+c(diff(std.dep)/2,50)),labels=std.dep)
D0$depth.bin<-x


D0$line<-D0$site %/%1000 /10
D0$station<-D0$site %% 1000
D0$seasonid<- ceiling((D0$monthid-2)/3)
D0$season<-D0$seasonid %% 4   # 0- winter; 1-spring; 2-summer; 3-autumn
D0<-D0[-which(D0$site %in% costal.site),]  #delete costal sites


########### get monthly data #############
name.monthid<-sort(unique(D0$monthid))
n.monthid<-length(name.monthid)
name.line<-(name.site %/% 1000)/10
name.station<-name.site %% 1000
n.site<-length(name.site)
name.depth<-sort(unique(D0$depth.bin))
n.depth<-length(name.depth)

Data<-data.frame(site=rep(NA,n.monthid*n.site*n.depth),monthid=rep(NA,n.monthid*n.site*n.depth),
                 depth.bin=rep(NA,n.monthid*n.site*n.depth),
                 temperature=rep(NA,n.monthid*n.site*n.depth), salinity=rep(NA,n.monthid*n.site*n.depth),
                 o2=rep(NA,n.monthid*n.site*n.depth), sio3=rep(NA,n.monthid*n.site*n.depth),
                 po4=rep(NA,n.monthid*n.site*n.depth), no3=rep(NA,n.monthid*n.site*n.depth),
                 no2=rep(NA,n.monthid*n.site*n.depth), chla=rep(NA,n.monthid*n.site*n.depth),
                 phaeo=rep(NA,n.monthid*n.site*n.depth))
a=1;
for(i in 1:n.site){
  D1<-subset(D0, line>=name.line[i]-0.666 & line<=name.line[i]+0.666 
             & station>=name.station[i]-2 & station<=name.station[i]+2)
  for(j in 1:n.monthid){
    D2<-subset(D1,monthid==name.monthid[j])
    for(k in 1:n.depth){
      D3<-subset(D2,depth.bin==name.depth[k])
      Data$site[a]=name.site[i]
      Data$monthid[a]=name.monthid[j]
      Data$depth.bin[a]=name.depth[k]
      
      if(nrow(D3)==0){
        Data[a,4:12]=NA
      }else{
        for(s in 4:12){
          Data[a,s]=mean(D3[,s],na.rm=T)
        }
      }
      a=a+1
    }
  }
}
rm(D0)
#write.csv(Data,'Data_month_bottle_depth.csv',row.names=F)

########### get seasonal data (csv) #############
Data$seasonid<- ceiling((Data$monthid-2)/3)
# 0 -> 195101-195102; 1->195103-195105; 2->195106-195108 ...
Data$season<-Data$seasonid %% 4   # 0- winter; 1-spring; 2-summer; 3-autumn
name.seasonid<-unique(Data$seasonid)
n.seasonid<-length(name.seasonid)

D.season<-data.frame(site=rep(NA,n.seasonid*n.site*n.depth),seasonid=rep(NA,n.seasonid*n.site*n.depth),
                     depth.bin=rep(NA,n.seasonid*n.site*n.depth),
                     temperature=rep(NA,n.seasonid*n.site*n.depth), salinity=rep(NA,n.seasonid*n.site*n.depth),
                     o2=rep(NA,n.seasonid*n.site*n.depth), sio3=rep(NA,n.seasonid*n.site*n.depth),
                     po4=rep(NA,n.seasonid*n.site*n.depth), no3=rep(NA,n.seasonid*n.site*n.depth),
                     no2=rep(NA,n.seasonid*n.site*n.depth), chla=rep(NA,n.seasonid*n.site*n.depth),
                     phaeo=rep(NA,n.seasonid*n.site*n.depth))
a=1
for(i in 1:n.site){
  D2<-subset(Data, site==name.site[i])
  for(k in 1:n.depth){
    D3<-subset(D2,depth.bin==name.depth[k])
    for(j in 1:n.seasonid){
      D4<-subset(D3, seasonid==name.seasonid[j])
      
      D.season$site[a]<-name.site[i]
      D.season$seasonid[a]<-name.seasonid[j]
      D.season$depth.bin[a]<-name.depth[k]
      for(s in 1:9){
        D.season[a,s+3]<-mean(D4[,s+3], na.rm=T)
      }
      a=a+1
    }
  }
}
rm(Data)
#write.csv(D.season,'Data_season_bottle_depth.csv',row.names=F)

############ get seasonal data (RDS) ##########

D.season<-D.season[-which(D.season$depth>200),]   #<<<<<<<<<<<<<<---------------
D.season$season <- D.season$seasonid %% 4
name.season<-c("winter","spring","summer","autumn")

site.shallow<-c(800051,867050,833051,833042) # delete these sites because: the first four don't have enough depth; 

D.season<-D.season[!D.season$site %in% site.shallow,]
D.season$yearid<- D.season$seasonid %/% 4 +1

spec<-c("temp","no3","chla")
n.spec<-length(spec)
Site<-sort(unique(D.season$site))
n.site<-length(Site)
Depth<-sort(unique(D.season$depth))
n.depth<-length(Depth)
Year<-range(D.season$year)[1]:range(D.season$year)[2]
n.year<-length(Year)
Station<-Site %% 1000
group<-cut(Station,c(0,55,70,90,120),c(1,2,3,4))

i.spec<-c(4,9,11)

D.org<-vector("list",4)
names(D.org)<-name.season
for(i.season in 1:4){
  D1<-subset(D.season, season==(i.season-1))    # 
  d1<-vector("list",n.spec)
  names(d1)<-spec
  for(s in 1:n.spec){
    d2<-vector("list",n.depth)
    names(d2)<-Depth
    for(d in 1:n.depth){
      D2<-subset(D1,depth.bin==Depth[d])
      d3<-matrix(NA,nrow=n.site, ncol=n.year)
      row.names(d3)<-Site
      colnames(d3)<-Year
      for(i in 1:n.site){
        D3<-subset(D2,site==Site[i])
        D3<-D3[order(D3$yearid),]
        
        tmp<-rep(NA,28)
        tmp[Year%in%D3$yearid]<-D3[,i.spec[s]]
        d3[i,]<-tmp
      }
      d2[[d]]<-d3
    }
    d1[[s]]<-d2
  }
  D.org[[i.season]]<-d1
}

#saveRDS(D.org,"Data_all_seasons_original.RDS")


############### interpolate ##############

D.int<-D.org

for(i.season in 1:4){
  for(s in 1:n.spec){
    for(d in 1:n.depth){
      for(i in 1:n.site){
        
        tmp<-approx(Year, D.org[[i.season]][[s]][[d]][i,], xout=Year,rule=2)
        
        D.int[[i.season]][[s]][[d]][i,]<-tmp$y
      }
    }
  }
}

saveRDS(D.int,"Data_all_seasons_interpolate.RDS")


