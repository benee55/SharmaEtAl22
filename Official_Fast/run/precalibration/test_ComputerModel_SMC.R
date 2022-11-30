rm(list=ls())
####################################################################################################
# Toy Example for Famos with resampling
####################################################################################################
# Load Libraries and Data
####################################################################################################
setwd("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/famos/Official_Fast/run/precalibration/")
library(snow);library(doParallel);library(foreach)
load("PrelimOutput.RData") # Data and Precalibration
load("MCMCOutput.RData") # MCMC Data
####################################################################################################
# Baby FAMOS
ens<-5000
famParMat<-cbind(runif(ens,rangeMat[1,1],rangeMat[1,2]) ,
              runif(ens,rangeMat[2,1],rangeMat[2,2]) ,
              rinvgamma(ens , shape = 0.2 , rate = 0.2))
famOutputMat<-apply(famParMat,1,compModel, x=inputX)

# Log Likelihood function
llhd_function<-function(x, obs){
  sum(dnorm(x = obs, mean = x[-1] , sd = sqrt(x[1]), log = TRUE))
}

# Calculate Weights for 0.1 tempering
foo<-rbind(famParMat[,3],famOutputMat)
llhd<-apply(foo, 2, llhd_function,obs=obs)
fullWeight<- llhd# first weights are log likelihood
weights<-fullWeight*0.5 # Same as L(theta)^0.1
weights<-exp(weights-max(weights))
weights<-weights/sum(weights)

# Resample particles 
reSampID<-sample(x=1:ens,size = ens, prob= weights, replace = TRUE)
length(unique(reSampID)) # Number of unique particles 
keepParMat<-famParMat[reSampID,] # Resampled particles
keepOutputMat<-apply(keepParMat[,1:2],1,compModel, x=inputX)

par(mfrow=c(1,2))
plot(x=inputX,y=keepOutputMat[,1],typ="n" , ylim=range(keepOutputMat),
     xlim=range(inputX),main="Posterior Model Runs")
thinID<-seq(1,nrow(keepParMat), length.out = 1000)
for(i in 1:length(thinID)){
  lines(x=inputX,y=keepOutputMat[,thinID[i]], lwd=0.01, col="blue")
}
lines(x=inputX, y= modelOutput , col="red" , lwd=2)


lines(x=inputX, y= meanPost , col="black" , lwd=2)
lines(x=inputX, y= hpdPost[1,] , col="black" , lwd=2 , lty=2)
lines(x=inputX, y= hpdPost[2,] , col="black" , lwd=2, lty=2)
points(x=inputX,y=obs,pch=16, col="red")


