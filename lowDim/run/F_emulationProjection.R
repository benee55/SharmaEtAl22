rm(list=ls())
# setwd("~/Dropbox/FamosHydroModel/lowDim/")
setwd("/glade/u/home/sanjib/FamosHydroModel/lowDim/")
# Load data
load("output/modelRuns_projection.RData")
# process output
modelRuns<-matrix(NA,nrow=length(outputMat[1,]), ncol=length(unlist(outputMat[1,1])))
for(k in 1:length(outputMat[1,])){
  modelRuns[k,]<-  unlist(outputMat[1,k])
}
# summary(modelRuns)

# load Design Matrix
load("input/design.RData")
parMat<-round(parMat,3)
# Build Emulator via paralellization
if(TRUE){
# Install on cheyenne by loading gnu module
library(snow);library(snowfall);library(mlegp)
sfInit(parallel=TRUE, cpus=7, type='PSOCK')
# How to do zero-mean GP?
gpEmulator_CM<-mlegp(X=parMat,
                  Z=modelRuns,
                  constantMean = 1,
                  nugget = 0,
                  parallel = TRUE)

save(gpEmulator_CM, modelRuns, parMat, 
     file="output/GPEmulator_projection.RData")  
# How to do zero-mean GP?
gpEmulator<-mlegp(X=parMat,
                  Z=modelRuns,
                  constantMean = 0,
                  nugget = 0,
                  parallel = TRUE)
sfStop()

save(gpEmulator,gpEmulator_CM, modelRuns, parMat, 
     file="output/GPEmulator_projection.RData")  
}

