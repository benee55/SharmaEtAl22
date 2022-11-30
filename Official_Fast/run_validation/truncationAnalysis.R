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
# Compile precalibration data - Validation

for(i in 1:14){
  
  if(i==1){
    load(paste("~/Dropbox/FamosHydroModel/Official_Fast/output_validation/preCalibrationResults",i,".RData",sep=""))
    bar<-unlist(outputMat[1,])
    modelOutputValidation<-matrix(bar, nrow=length(outputMat[1,]), ncol=length(outputMat[1,1][[1]]), byrow = TRUE)
    
  }else{
    load(paste("~/Dropbox/FamosHydroModel/Official_Fast/output_validation/preCalibrationResults",i,".RData",sep=""))
    bar<-unlist(outputMat[1,])
    foo<-matrix(bar, nrow=length(outputMat[1,]), ncol=length(outputMat[1,1][[1]]), byrow = TRUE)
    modelOutputValidation<-rbind(modelOutputValidation,foo)
  }
}


# Format Date
dateVect<-paste(sprintf("%04d",as.numeric(obs[,1])),
                sprintf("%02d",as.numeric(obs[,2])),
                sprintf("%02d",as.numeric(obs[,3])),sep="-")
dateVect<-as.Date(dateVect, format = "%Y-%m-%d")
2003/06/01-2008/03/31

modelStart<-which(dateVect==as.Date("2003-06-01")) # Remove Spin-up
modelEnd<-which(dateVect==as.Date("2008-03-31")) # Trim Last Days
newModelOutput<-modelOutput[modelStart:modelEnd,] # Trimmed Data
newDateVect<-dateVect[modelStart:modelEnd] # Trimmed Dates
flowObs<-obs[modelStart:modelEnd,4]

# observation Index
extremeDate<-dateVect[obsInd] # Extreme Dates
extremeObs<-subsetFinalObs # Extreme Values
validationDate<-dateVect[validationInd] # Extreme Dates
validationObs<-subsetFinalValidation # Extreme Values
extremeIndex<-which(newDateVect%in%extremeDate)

####################################################################################
# Which model runs are above activation?
greater5k<-apply(modelOutput[extremeIndex,],2,function(x){sum(x>4950.55)})
largeModelIndex<-which(greater5k>=15)
summary(largeModelIndex)
largeParMat<-parMat[largeModelIndex,-1]


extremeLargeModelOutput<-modelOutput[extremeIndex,largeModelIndex]
extremeLargeValidationOutput<-t(modelOutputValidation[largeModelIndex,])
scoreVal<-apply(extremeLargeModelOutput,2,function(x){
  sum(dnorm(x=extremeObs , mean=x , sd=2000, log=TRUE))
})
goodRuns<-which(scoreVal>quantile(scoreVal, probs=0.95, na.rm = TRUE))



# Parameters
boundMat<-rbind(c(0, 5) , # PCTIM 0.3=original maximum
                c(0 , 2), # ADIMP 0.5=original maximum
                c(-50 , -0.1), # UZTWM
                c(-70 , -0.1), # LZTWM
                c(-100 , -0.1), # LZFSM
                c(-100 , -0.1), # LZFPM Old -120
                c(-3.8 , -0.1), # LZSK
                c(0.5 , 1.5), # snow_SCF
                c(-3.5 , -0.1), # REXP
                c(-3.5 , -0.1), # UZK
                c(0.5,4.5), # rutpix_Q0CHN
                c(0.3,1.9)) # rutpix_QMCHN Use Original: 3.4 ; BPrior: 2.25 
parNames<-c("PCTIM" , "ADIMP" , "UZTWM" ,"LZTWM" , 
            "LZFSM" , "LZFPM" , "LZSK" , "snow_SCF" , 
            "REXP" , "UZK" , "Q0CHN" , "QMCHN")
par(mfrow=c(3,4), mar=c(2,2,2,2))
for(k in 1:ncol(largeParMat)){
  d1<-density(largeParMat[,k])
  plot(d1, xlim=range(d1$x,boundMat[k,]), main=parNames[k])
  abline(v=boundMat[k,], col="red")
}

par(mfrow=c(3,4), mar=c(2,2,2,2))
for(k in 1:ncol(largeParMat)){
  d1<-density(largeParMat[goodRuns,k])
  plot(d1, xlim=range(d1$x,boundMat[k,]), main=parNames[k])
  abline(v=boundMat[k,], col="red")
}

par(mfrow=c(5,5), mar=c(2,2,2,2))
for(h in 1:length(extremeDate)){
  useDat<-extremeLargeModelOutput[h,goodRuns]
  plot(density(useDat), xlim=range(useDat, extremeObs[h], 4950.55), main=extremeDate[h])
  abline(v=4950.55, col="red", lty=2)
  abline(v=extremeObs[h], col="blue", lty=2)
  
}

par(mfrow=c(5,4), mar=c(2,2,2,2))
for(h in 1:length(validationDate)){
  useDat<-extremeLargeValidationOutput[h,goodRuns]
  plot(density(useDat), xlim=range(useDat, extremeObs[h], 4950.55), main=validationDate[h])
  abline(v=4950.55, col="red", lty=2)
  abline(v=validationObs[h], col="blue", lty=2)
  
}
