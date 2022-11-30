rm(list=ls())
setwd("~/Dropbox/FamosHydroModel/")
load("Official_Fast/input/fullObservations.RData") # Load Full observation

load("Analysis/resultsStremflow_calibration_full.RData")
dateVect<-dateVect[obsInd]
famosOutput4Par<-famosOutput4Par[,obsInd]
famosOutput<-famosOutput[,obsInd]
precalibrationOutput_MSE<-precalibrationOutput_MSE[,obsInd]
precalibrationOutput_Window<-precalibrationOutput_Window[,obsInd]
handTuneOutput<-handTuneOutput[obsInd]
  
save(dateVect , subsetFinalObs , famosOutput4Par ,
     famosOutput , precalibrationOutput_MSE , precalibrationOutput_Window,
     handTuneOutput,
     file="Analysis/resultsStremflow_calibration_extreme.RData")


rm(list=ls())
setwd("~/Dropbox/FamosHydroModel/")
load("Official_Fast/input/fullObservations.RData") # Load Full observation

load("Analysis/resultsStremflow_validation_full.RData")
# Format Date
dateVect<-paste(sprintf("%04d",as.numeric(obs[,1])),
                sprintf("%02d",as.numeric(obs[,2])),
                sprintf("%02d",as.numeric(obs[,3])),sep="-")
dateVect<-as.Date(dateVect, format = "%Y-%m-%d")
validationDate<-dateVect[validationInd] # Extreme Dates

dateVect<-dateVect[which(dateVect=="2009-01-01"):which(dateVect=="2011-10-01")]
# observation Index
keepValidation<-which(dateVect%in%validationDate)
dateVect<-dateVect[keepValidation]
famosOutput4Par<-famosOutput4Par[,keepValidation]
famosOutput<-famosOutput[,keepValidation]
precalibrationOutput_MSE<-precalibrationOutput_MSE[,keepValidation]
precalibrationOutput_Window<-precalibrationOutput_Window[,keepValidation]
handTuneOutput<-handTuneOutput[keepValidation]

save(dateVect , subsetFinalValidation , famosOutput4Par,
     famosOutput , precalibrationOutput_MSE , precalibrationOutput_Window,
     handTuneOutput,
     file="Analysis/resultsStremflow_validation_extreme.RData")

rm(list=ls())
