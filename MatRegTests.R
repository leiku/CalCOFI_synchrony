#A function for doing matrix regression tests of a model against a nested model.
#
#Args
#resp         A response variable matrix, square
#preds        A list of predictor matrices of the same dimensions as resp
#drop         A vector of indices of the predictors in preds to drop in the 
#               simpler model
#numperm      Number of permutations on which to base the test
#
#Output - a named list with these elements
#ssr_dat      sum of squared residuals for the linear regression using all predictors
#ssr_perm     vector of sum of squares of residuals for linear regressions using
#               randomized matrices for the indices "drop"
#p            p value of the test of the simpler model as a null hypothesis against
#               the more complex model as an alternative 
#
#Note: no error checking is done so this is a function for our own use
#
matregtest<-function(resp,preds,drop,numperm)
{
  n<-dim(resp)[1]
  
  #A) do a regression of triang(resp) against the vectors triang(preds[[n]]) for 
  #n in 1:length(preds), and save the sum of squared residuals
  resp.vec<-triang(resp)
  preds.vec<-lapply(preds, triang)
  form<-"resp.vec~preds.vec[[1]]"
  if(length(preds.vec)>1){
    for(i in 2:length(preds.vec)){
      form<-paste(form, paste0("preds.vec[[",i,"]]"), sep="+")
    }
  }
  form<-as.formula(form)
  res.lm<-lm(formula=form, na.action=na.omit)
  ssres<- sum(res.lm$residuals^2)
  
  #B) apply the randomization to the matrices that are dropped in the simpler model 
  #and then repeat the regression, numperm times
  ssres_perm<-NA*numeric(numperm)
  for (permcount in 1:numperm)
  {
    #randomize the dropped variables
    perm<-sample.int(n,n) #this is the randomization
    for (predcount in 1:length(drop))
    {
      preds[[drop[predcount]]]<-preds[[drop[predcount]]][perm,perm]
    }
    preds.vec[drop]<-lapply(preds[drop], triang)
    res.lm<-lm(formula=form, na.action=na.omit)
    ssres_perm[permcount]<-sum(res.lm$residuals^2)
  }
  
  #C) the results as follows
  p<-sum(ssres_perm<ssres)/numperm
  return(list(ssr_dat=ssres,ssr_perm=ssres_perm,p=p))
}
  
#Extracts the lower-triangular portion of the square matrix m (no diagonal)
#as a vector, UNLESS the lower-triangle is all 0s or NAs, in which case it returns 
#the upper triangle.
#
#Note: no error checking is done so this is a function for our own use
#
triang<-function(m)
{
  if(all(m[lower.tri(m)]==0) | all(is.na(m[lower.tri(m)]))){
    return(as.vector(m[upper.tri(m)])) #if m is a upper triangular matrix, return the upper triangular matrix
  }else{return(as.vector(m[lower.tri(m)]))} #else return the lower triangular matrix
}
