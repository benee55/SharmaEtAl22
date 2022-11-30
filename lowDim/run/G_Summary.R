
rm(list=ls())
library(mlegp)
setwd("~/Dropbox/FamosHydroModel/lowDim/")

# Try out Emulator

load("output/mcmcResults.RData")
sampID<-round(seq(50000,100000, length.out = 10000))
testMat<-amcmc.out$samples[sampID,1:4]
# PRedict on test parameters
CMpredTest<-predTest<-matrix(NA, nrow=nrow(testMat) , ncol=21)
CMpredTest_proj<-predTest_proj<-matrix(NA, nrow=nrow(testMat) , ncol=18)
load("output/GPEmulator_Full.RData")
for(k in 1:21){
  print(k)
  predTest[,k]<-apply(testMat, 1, predict, object=gpEmulator[[k]])
  CMpredTest[,k]<-apply(testMat, 1, predict, object=gpEmulator_CM[[k]])
}
rm(gpEmulator_CM)
load("output/GPEmulator_projection.RData")
for(k in 1:18){
  print(k)
  predTest_proj[,k]<-apply(testMat, 1, predict, object=gpEmulator[[k]])
  CMpredTest_proj[,k]<-apply(testMat, 1, predict, object=gpEmulator_CM[[k]])
}


save(testMat,
     predTest,CMpredTest , predTest_proj , CMpredTest_proj, file="output/finalEmulationResults.RData")


#load observations
load("input/fullObservations.RData")
subsetFinalObs

# Format Date
dateVect<-paste(sprintf("%04d",as.numeric(obs[,1])),
                sprintf("%02d",as.numeric(obs[,2])),
                sprintf("%02d",as.numeric(obs[,3])),sep="-")
dateVect<-as.Date(dateVect, format = "%Y-%m-%d")

# observation Index
extremeDate<-dateVect[obsInd] # Extreme Dates
validationDate<-dateVect[validationInd] # Extreme Dates

load("../Official_Fast/output_validation/handTuneResults.RData")
handTuneOutput<-outputMat[[1]]
# PLot Results - Projection
par(mfrow=c(4,5), mar=c(2,2,2,2))
for(i in 1:18){
  d1<-density(predTest_proj[,i])
  plot(d1, xlim=range(d1$x,subsetFinalValidation[i],handTuneOutput[i], na.rm = TRUE), 
       main = validationDate[i])
  abline(v=subsetFinalValidation[i], col="blue" ,pch=16)
  abline(v=handTuneOutput[i], col="red" ,lwd=1 , lty=2) # ACtion Stage
}

load("../Official_Fast/output_validation/handTuneResults.RData")
handTuneOutput<-modelRun[[1]]
# PLot Results - Model Fitting
par(mfrow=c(5,5), mar=c(2,2,2,2))
for(i in 1:21){
  d1<-density(predTest[,i])
  plot(d1, xlim=range(d1$x,subsetFinalObs[i],handTuneOutput[i], na.rm = TRUE), 
       main = extremeDate[i])
  abline(v=subsetFinalObs[i], col="blue" ,pch=16)
  abline(v=handTuneOutput[i], col="red" ,lwd=1 , lty=2) # ACtion Stage
}
