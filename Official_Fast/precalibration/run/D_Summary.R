# Precalibration Summary
rm(list=ls())
library(vioplot)
# Set working directory and load Full observations
setwd("~/Dropbox/FamosHydroModel/Official_Fast/output/")
load("../input/fullObservations.RData") # Load Full observation
load("../precalibration/output/mhParameters_0.RData")
# source("../run/mcmc_source_Tr.R")
# Compile precalibration data
for(i in 1:14){
  load(paste("../precalibration/output/preCalibrationResults",i,".RData",sep=""))
  if(i==1){
    modelOutput<-outputMat[-nrow(outputMat),]
  }else{
    modelOutput<-cbind(modelOutput,outputMat[-nrow(outputMat),])
  }
}
modelOutput<-modelOutput[,-ncol(modelOutput)] # Remove last model run (overlap with below)


# Format Date
dateVect<-paste(sprintf("%04d",as.numeric(obs[,1])),
                sprintf("%02d",as.numeric(obs[,2])),
                sprintf("%02d",as.numeric(obs[,3])),sep="-")
dateVect<-as.Date(dateVect, format = "%Y-%m-%d")

# observation Index
extremeDate<-dateVect[obsInd] # Extreme Dates
extremeObs<-subsetFinalObs # Extreme Values
modelStart<-which(dateVect==as.Date("2004-08-01")) # Remove Spin-up
modelEnd<-which(dateVect==as.Date("2008-03-31")) # Trim Last Days
keepInd<-1:ncol(modelOutput) # Can use sample() to select subset
newModelOutput<-modelOutput[modelStart:modelEnd,keepInd] # Trimmed Data
newDateVect<-dateVect[modelStart:modelEnd] # Trimmed Dates


########################
par(mfrow=c(5,5), mar=c(2,2,2,2))
for(i in 1:21){
  k<-obsInd[i]
  vioplot(modelOutput[k,], ylim=range(modelOutput[k,],extremeObs[i],4950.55, na.rm=TRUE), 
          main = extremeDate[i])
  points(x=1, y=extremeObs[i], col="red" ,pch=16)
  abline(h=4950.55, col="red" ,lwd=1 , lty=2) # ACtion Stage
}


########################
parMat<-parMat[1:ncol(modelOutput),-1]
par(mfrow=c(4,3), mar=c(2,2,2,2))
parNames<-c("PCTIM" , "ADIMP" , "UZTWM" ,"LZTWM" , 
            "LZFSM" , "LZFPM" , "LZSK" , "snow_SCF" ,
            "REXP" , "UZK" , "Q0CHN" , "QMCHN")
boundMat<-rbind(c(0, 5) , # PCTIM 0.3=original maximum
                c(0 , 2), # ADIMP 0.5=original maximum
                c(-50 , -0.1), # UZTWM
                c(-70 , -0.1), # LZTWM
                c(-100 , -0.1), # LZFSM
                c(-100 , -0.1), # LZFPM
                c(-3.8 , -0.1), # LZSK
                c(0.5 , 1.5), # snow_SCF
                c(-3.5 , -0.1), # REXP
                c(-3.5 , -0.1), # UZK
                c(0.5,4.5), # rutpix_Q0CHN
                c(0.3,2.25)) # rutpix_QMCHN Use 2.25 instead of 3.4

########################
########################
extremeModelOuput<-modelOutput[obsInd,]
MSE<-apply(extremeModelOuput,2,function(x){mean((x-extremeObs)^2)})
goodRuns<-which(MSE<quantile(MSE, probs=0.05, na.rm = TRUE))
goodModelOutput<-modelOutput[modelStart:modelEnd,goodRuns] # Trimmed Data

########################
########################
par(mfrow=c(5,5), mar=c(2,2,2,2))
for(i in 1:21){
  k<-obsInd[i]
  vioplot(modelOutput[k,goodRuns], ylim=range(modelOutput[k,goodRuns],extremeObs[i],4950.55), 
          main = extremeDate[i])
  points(x=1, y=extremeObs[i], col="red" ,pch=16)
  abline(h=4950.55, col="red" ,lwd=1 , lty=2) # ACtion Stage
}



########################
par(mfrow=c(4,3), mar=c(2,2,2,2))
for(i in 1:11){
  plot(density(parMat[goodRuns,i]) , xlim=range(boundMat[i,1:2]),  main=parNames[i])
  abline(v=boundMat[i,1:2], col="red")
}
