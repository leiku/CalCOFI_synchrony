#This does unit tests for functions used by master.R or the scripts it calls.
#
#Lei Zhao
#Started 2018 04 09

#***Test functins in MatRegTests.R
rm(list=ls())
source("MatRegTests.R")

#first test triang
m<-matrix(c(1,2,3,4,5,6,7,8,9),3,3)
h<-triang(m)
if ((length(h)==3) && (sum(h==c(2,3,6))==3))
{
  print("passed")
} else
{
  print("failed triang in MatRegTests.R")
}

#Lei to add a test case for the case where the lower triangle is all 0s


#Lei to also add a unit test for matregtest


#***Lei to add a unit tests for the new version of summed.weights (even though that is in the package)


#***Lei to add unit tests for any other code for which he has doubts

