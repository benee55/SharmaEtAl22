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
dateVect<-dateVect[which(dateVect=="2009-01-01"):which(dateVect=="2011-10-01")]

##########################################################################################
##########################################################################################
# Compile Data
##########################################################################################
##########################################################################################

### 1. Famous Full
load("Official_Fast/output_validation/famosResults_rep_2k.RData")
famosOutput<-matrix(unlist(outputMat[2,]),
                    nrow=length(outputMat[2,]), 
                    ncol=length(outputMat[2,][[1]]), byrow = TRUE)

### 1. Famous 4 parameter
for(i in 1:58){
  print(i)
  if(i==1){
    load(paste("lowDim/output_f/validationResults",i,".RData",sep=""))
    bar<-unlist(outputMat[2,])
    famosOutput4Par<-matrix(bar, nrow=length(outputMat[2,]), ncol=length(outputMat[2,1][[1]]), byrow = TRUE)
    
  }else{
    load(paste("lowDim/output_f/validationResults",i,".RData",sep=""))
    bar<-unlist(outputMat[2,])
    foo<-matrix(bar, nrow=length(outputMat[2,]), ncol=length(outputMat[2,1][[1]]), byrow = TRUE)
    famosOutput4Par<-rbind(famosOutput4Par,foo)
  }
}

dim(famosOutput4Par)

###  4. Hand Tune
load("Official_Fast/output_validation/handTuneResults.RData")
handTuneOutput<-outputMat[[2]]

###  5. Pre-calibration
# Compile precalibration data
for(i in 1:14){
  
  if(i==1){
    load(paste("Official_Fast/output_validation/preCalibrationResults",i,".RData",sep=""))
    bar<-unlist(outputMat[2,])
    modelOutput<-matrix(bar, nrow=length(outputMat[2,]), ncol=length(outputMat[2,1][[1]]), byrow = TRUE)
    
  }else{
    load(paste("Official_Fast/output_validation/preCalibrationResults",i,".RData",sep=""))
    bar<-unlist(outputMat[2,])
    foo<-matrix(bar, nrow=length(outputMat[2,]), ncol=length(outputMat[2,1][[1]]), byrow = TRUE)
    modelOutput<-rbind(modelOutput,foo)
  }
}
load("Official_Fast/output_validation/precalibration_goodRuns.RData")
precalibrationOutput_MSE<-modelOutput[goodRuns_MSE,]
precalibrationOutput_Window<-modelOutput[goodRuns_Window,]


##########################################################################################
##########################################################################################
# Save files (for Iman)
names(handTuneOutput)<-colnames(precalibrationOutput_MSE)<-colnames(precalibrationOutput_Window)<-colnames(famosOutput)<-dateVect
save(dateVect , subsetFinalValidation , famosOutput4Par,
     famosOutput , famosOutput , precalibrationOutput_MSE , precalibrationOutput_Window,
     handTuneOutput,
     file="Analysis/resultsStremflow_validation_full.RData") 

# save(dateVect , subsetFinalValidation , famosOutput4Par,
#      file="Analysis/famos4parameter_validation.RData") 