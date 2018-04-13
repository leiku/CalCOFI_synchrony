#This does unit tests for functions used by master.R or the scripts it calls.
#
#Lei Zhao
#Started 2018 04 09

#*** function to print message
print.message<-function(name="test", condit=T){
  if(condit){
    print(paste0("passed ",name))
  }else{print(paste0("failed ",name))}
}


#***Test functins in MatRegTests.R
rm(list=ls())
source("MatRegTests.R")

#first test triang
m<-matrix(c(1,2,3,4,5,6,7,8,9),3,3)
h<-triang(m)
condit1<-(length(h)==3) && (sum(h==c(2,3,6))==3)

#a test case for the case where the lower triangle is all 0s
m<-matrix(1:9,3,3)
m[lower.tri(m)]<-0
h<-triang(m)
condit2<-sum(h==c(4,7,8))==3

#a test case for the case where the lower triangle is all NAs
m<-matrix(1:9,3,3)
m[lower.tri(m)]<-NA
h<-triang(m)
condit3<-sum(h==c(4,7,8))==3

print.message(condit=condit1 && condit2 && condit3, name="triang in MatRegTests.R")

#a unit test for matregtest
x1<-matrix(runif(25,0,1),5,5)
x2<-matrix(runif(25,0,1),5,5)
y<-2*x1
z1<-matregtest(y,list(x1,x2), 1, 1000)
z2<-matregtest(y,list(x1,x2), 2, 1000)
condit1<-(z1$p<0.001) && (z2$p>0.05)

y<-2*x1+3*x2
z1<-matregtest(y,list(x1,x2), 1, 1000)
z2<-matregtest(y,list(x1,x2), 2, 1000)
condit2<-(z1$p<0.001) && (z2$p<0.001)

print.message(condit=condit1 && condit2, name="matregtest in MagRegTests.R")

#***Lei to add a unit tests for the new version of summed.weights (even though that is in the package)
modelnames<-list(2, c(3,4), c(2,3,4))
x1<-matrix(rnorm(100), 10, 10)
x2<-matrix(rnorm(100), 10, 10)
x3<-matrix(rnorm(100), 10, 10)
y<-2*x2+x3
X<-list(y=y,x1=x1,x2=x2,x3=x3)
z<-lno.weights(X,model.names=modelnames,n=3,maxruns=5,nrand=10) 
z1<-summed.weights(varnames=names(X),weights=z)

print.message(condit=, name="summed.weights in Reumannplatz.R")

#***Lei to add unit tests for any other code for which he has doubts

