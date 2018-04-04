#Just returns the synchrony matrix from whatever method used.
#
wsynmatSimple.mod.pearson<-function(indata, times, tsrange=c(0,Inf),  
                            scale.min=2, scale.max.input=NULL, sigma=1.05, f0 = 0.5, alpha=0.05)
{
  #tiny bit of error checking
  if (any(!is.finite(indata)))
  {
    stop("Error in wsynmatSimple: non-finite entries in indata")
  }
  
  nlocs<-nrow(indata)
  ntimes<-ncol(indata)


  synmat<-cor(t(indata), method="pearson")
  diag(synmat)<-0

  return(synmat)
} 

#***Should be done (enough) with the above, need to do below

SynNetModulesSimple.mod.pearson<-function(indata, times, coords, tsrange=c(0,Inf), 
                                  scale.min=2, scale.max.input=NULL, sigma=1.05, f0 = 1,
                                  filename,alpha=0.05, do.plot=FALSE)
{
  library(igraph)
  
  #get the synchrony matrix and corresponding network and modules
  synmat<-wsynmatSimple.mod.pearson(indata,times,
                            tsrange=tsrange,scale.min,scale.max.input,
                            sigma,f0,alpha)
  synnet<-graph.adjacency(synmat, weighted=T, mode="undirected", diag=F)
  #modules<-cluster_leading_eigen(synnet)
  modules<-newman_eigenvector_mod(synmat)
  
  
  #get mean time series and wmfs and wpmfs for modules
  module.means<-list()
  module.wmfs<-list()
  module.wpmfs<-list()
  for(mm in 1:max(modules$membership)){
    if(sum(modules$membership==mm)>=3){
      name<-paste0("module",mm)
      module.means[[name]]<-colMeans(indata[modules$membership==mm,])
     
    }
  }
  
  if(do.plot){
    #make the module map plot
    #pdf(file=paste(filename,"_ModuleMap.pdf",sep=''))
    #par(mfrow=c(1,1))
    #wt<-E(synnet)$weight
    #if (length(wt)>1000)
    #{
    #  swt<-sort(wt)
    #  whe<-which(wt<swt[length(wt)-1000])
    #  dsynnet<-delete.edges(synnet,whe)
    #}
    #plot(dsynnet,
    #     vertex.color=modules$membership,layout=coords,vertex.label=NA,
    #     vertex.size=2,vertex.frame.color=modules$membership)
    #legend(x="bottomleft",legend=names(module.means))
    #dev.off()
    
    #make the plot of module mean time series
    pdf(file=paste(filename,"_ModuleMeans.pdf",sep=''))
    par(mfrow=c(length(module.means),1))
    for(mm in 1:length(module.means)){
      plot(times, module.means[[mm]], type="l", main=names(module.means)[mm])
    }
    dev.off()
    
  }
  
  
  return(list(synnet=synnet, modules=modules, module.means=module.means, synmat=synmat))
}
