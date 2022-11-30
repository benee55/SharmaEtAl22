# Precalibration Summary
rm(list=ls())
library(vioplot)
# Set working directory and load Full observations
setwd("~/Dropbox/FamosHydroModel/Official_Fast/output_validation/")
load("../input/fullObservations.RData") # Load Full observation
load("../precalibration/output/mhParameters_0.RData")
load("precalibration_goodRuns.RData")
# Compile precalibration data
for(i in 1:14){
  
  if(i==1){
    load(paste("preCalibrationResults",i,".RData",sep=""))
    bar<-unlist(outputMat[1,])
    modelOutput<-matrix(bar, nrow=length(outputMat[1,]), ncol=length(outputMat[1,1][[1]]), byrow = TRUE)
    
  }else{
    load(paste("preCalibrationResults",i,".RData",sep=""))
    bar<-unlist(outputMat[1,])
    foo<-matrix(bar, nrow=length(outputMat[1,]), ncol=length(outputMat[1,1][[1]]), byrow = TRUE)
    modelOutput<-rbind(modelOutput,foo)
  }
}


# Format Date
dateVect<-paste(sprintf("%04d",as.numeric(obs[,1])),
                sprintf("%02d",as.numeric(obs[,2])),
                sprintf("%02d",as.numeric(obs[,3])),sep="-")
dateVect<-as.Date(dateVect, format = "%Y-%m-%d")
  
# observation Index
extremeDate<-dateVect[obsInd] # Extreme Dates
extremeObs<-subsetFinalObs # Extreme Values
validationDate<-dateVect[validationInd] # Extreme Dates
validationObs<-subsetFinalValidation # Extreme Values

# modelStart<-which(dateVect==as.Date("2010-01-01")) # Remove Spin-up
# modelEnd<-which(dateVect==as.Date("2011-11-01")) # Trim Last Days
# keepInd<-1:ncol(modelOutput) # Can use sample() to select subset
# newModelOutput<-modelOutput[modelStart:modelEnd,keepInd] # Trimmed Data
# newDateVect<-dateVect[modelStart:modelEnd] # Trimmed Dates



# ADD HandTune
load("handTuneResults.RData")
handTuneOutput<-outputMat[[1]]
rm(modelRun, outputMat)
# ADD FAMOS
load("famosResults_rep_2k.RData")
famosOutput<-matrix(unlist(outputMat[1,]), nrow=length(outputMat[1,]), ncol=length(handTuneOutput), byrow = TRUE)


# Plot Extreme Dates
par(mfrow=c(4,5), mar=c(2,2,2,2))
for(k in 1:ncol(famosOutput)){
  d1<-density(famosOutput[,k])
  d2<-density(modelOutput[goodRuns,k])
  plot(d1, main=validationDate[k] , 
       xlim=range(d1$x, d2$x , validationObs[k] , handTuneOutput[k]), 
       ylim=range(d1$y, d2$y))
  lines(d2, col="gray")
  abline(v=validationObs[k] , col="blue")
  abline(v=handTuneOutput[k] , col="red")
}



