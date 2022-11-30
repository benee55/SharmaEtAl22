# Generate Pre-calibration Samples

rm(list=ls())
setwd("/glade/u/home/sanjib/FamosHydroModel/lowDim")
# setwd("~/Dropbox/FamosHydroModel/lowDim/")
source("run/mcmc_source_Tr.R")
source("run/rWrapper_Continuous.R")

J=5
ensembleN<-J^4
seqInterval<-t(apply(boundMat,1,function(x,j){seq(x[1],x[2],length.out = j)}, j=J))
testN<-71
parMat<-expand.grid(seqInterval[1,],seqInterval[2,], seqInterval[3,],seqInterval[4,])
testMat<-apply(boundMat,1,function(x,ens){runif(ens,min=x[1],max=x[2])},ens=testN)
save(parMat, testMat,file="input/design.RData")

