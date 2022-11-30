# Methods
# 1. Famous Full
# 2. Famous 4-D
# 3. Emulation-Calibration 4-D
# 4. Hand Tune
# 5. Pre-calibration
########################################################################################
########################################################################################
#Preliminaries
########################################################################################
########################################################################################
rm(list=ls())
setwd("~/Dropbox/FamosHydroModel/")
load("Official_Fast/input/fullObservations.RData") # Load Full observation
# Format Date
dateVect<-paste(sprintf("%04d",as.numeric(obs[,1])),
                sprintf("%02d",as.numeric(obs[,2])),
                sprintf("%02d",as.numeric(obs[,3])),sep="-")
dateVect<-as.Date(dateVect, format = "%Y-%m-%d")

dateVect<-dateVect[which(dateVect=="2003-06-01"):which(dateVect=="2008-03-31")]



##########################################################################################
##########################################################################################
# Compile Data
##########################################################################################
##########################################################################################

### 1. Famous Full
load("Official_Fast/output_rep/mhParameters_4.RData")
famosOutput<-matrix(NA,
                    nrow=length(initResultsList[[2]]), 
                    ncol=length(initResultsList[[2]][[1]][[2]]), byrow = TRUE)
for(k in 1:length(initResultsList[[2]])){
  famosOutput[k,]<-initResultsList[[2]][[k]][[2]]
}

### 3. Famous Four Parameter
load("lowDim/output_f/mhParameters_4.RData")
famosOutput4Par<-matrix(NA,
                    nrow=length(initResultsList[[2]]), 
                    ncol=length(initResultsList[[2]][[1]][[2]]), byrow = TRUE)
for(k in 1:length(initResultsList[[2]])){
  famosOutput4Par[k,]<-initResultsList[[2]][[k]][[2]]
}

###  4. Hand Tune
load("Official_Fast/output_validation/handTuneResults.RData")
handTuneOutput<-modelRun[[2]]

###  5. Pre-calibration
# Compile precalibration data
for(i in 1:14){
  load(paste("Official_Fast/precalibration/output/preCalibrationResults",i,".RData",sep=""))
  if(i==1){
    modelOutput<-outputMat[-nrow(outputMat),]
  }else{
    modelOutput<-cbind(modelOutput,outputMat[-nrow(outputMat),])
  }
}

modelOutput<-t(modelOutput)
load("Official_Fast/output_validation/precalibration_goodRuns.RData")
precalibrationOutput_MSE<-modelOutput[goodRuns_MSE,]
precalibrationOutput_Window<-modelOutput[goodRuns_Window,]

##########################################################################################
##########################################################################################
# Visualization - Carlibration Dates Densities
##########################################################################################
##########################################################################################
load("Official_Fast/input/fullObservations.RData") # Load Full observation
##########################################################################################
##########################################################################################
# Save files (for Iman and Sanjib)
names(handTuneOutput)<-colnames(precalibrationOutput_MSE)<-colnames(precalibrationOutput_Window)<-colnames(famosOutput)<-dateVect
save(dateVect , subsetFinalObs , famosOutput4Par , 
     famosOutput , precalibrationOutput_MSE , precalibrationOutput_Window,
     handTuneOutput,
     file="Analysis/resultsStremflow_calibration_full.RData") 
