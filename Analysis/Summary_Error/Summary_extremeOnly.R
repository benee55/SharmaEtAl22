rm(list=ls())
setwd("~/Dropbox/FamosHydroModel/")
load("Official_Fast/input/fullObservations.RData") # Load Full observation

load("manuscript/revisionCode/resultsStremflow_calibration_full.RData")
dateVect<-dateVect[obsInd]
famosOutput_rep_var<-famosOutput_rep_var[,obsInd]
famosOutput<-famosOutput[,obsInd]
precalibrationOutput_MSE<-precalibrationOutput_MSE[,obsInd]
precalibrationOutput_Window<-precalibrationOutput_Window[,obsInd]
handTuneOutput<-handTuneOutput[obsInd]
  
save(dateVect , subsetFinalObs , famosOutput_rep_var ,
     famosOutput , precalibrationOutput_MSE , precalibrationOutput_Window,
     handTuneOutput,
     file="manuscript/revisionCode/resultsStremflow_calibration_extreme.RData")


rm(list=ls())
setwd("~/Dropbox/FamosHydroModel/")
load("Official_Fast/input/fullObservations.RData") # Load Full observation

load("manuscript/revisionCode/resultsStremflow_validation_full.RData")
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
famosOutput_rep_var<-famosOutput_rep_var[,keepValidation]
famosOutput<-famosOutput[,keepValidation]
precalibrationOutput_MSE<-precalibrationOutput_MSE[,keepValidation]
precalibrationOutput_Window<-precalibrationOutput_Window[,keepValidation]
handTuneOutput<-handTuneOutput[keepValidation]

save(dateVect , subsetFinalValidation , famosOutput_rep_var,
     famosOutput , precalibrationOutput_MSE , precalibrationOutput_Window,
     handTuneOutput,
     file="manuscript/revisionCode/resultsStremflow_validation_extreme.RData")

rm(list=ls())
