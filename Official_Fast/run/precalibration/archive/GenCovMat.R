rm(list=ls())
setwd("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/famos/Official_Fast/")
source("run/mcmc_source_Tr.R")
source("run/rWrapper.R")
load("precalibration//preCalibrationResults.RData")
load("input/obsData.RData")
load("precalibration//mhParameters_0.RData")

# LIkeihood
llhd<-apply(outputMat, 1, function(x){sum(dnorm(x = obs, mean = x , sd = 1000, log = TRUE))})
covIndex<-which(llhd>quantile(llhd,probs = 0.9))

keepParMat<-parMat[covIndex,]
# 
# par(mfrow=c(4,3), mar=c(2,2,2,2))
# for(i in 2:13){
#   plot(density(keepParMat[,i]) , xlim=range(boundMat[i-1,1:2]),  main=parNames[i])
#   abline(v=boundMat[i-1,1:2], col="red")
# }



CovMat<-cov(keepParMat[,-1])*((2.38^2)/ncol(keepParMat[,-1]))

save(CovMat, file="~/scratch/famos/output/BeginCovMat_Tr.RData")
