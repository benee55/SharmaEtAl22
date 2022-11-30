# Methods
# 1. Famous Full
# 2. Famous Variable Variance
# 3. Hand Tune
# 4. Pre-calibration
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
# observation Index
extremeDate<-dateVect[obsInd] # Extreme Dates
extremeObs<-subsetFinalObs # Extreme Values

##########################################################################################
##########################################################################################
# Compile Data
##########################################################################################
##########################################################################################

### 1. Famous Full
load("Official_Fast/output_rep/mhParameters_4.RData")
famosOutput<-matrix(NA,
                    nrow=length(initResultsList[[2]]), 
                    ncol=length(initResultsList[[2]][[1]][[1]]), byrow = TRUE)
for(k in 1:length(initResultsList[[2]])){
  famosOutput[k,]<-initResultsList[[2]][[k]][[1]]
}

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
rep2orig<-function(par,boundMat){ # Convert from reparameterized to original parameters
  c(par[1],par[-1]*(boundMat[,2]-boundMat[,1])/10+boundMat[,1])
}
famosParMat<-t(apply(parMat,1,rep2orig,boundMat=boundMat))


### 2. Famous - Estimated Variance
load("Official_Fast/output_rep_var/mhParameters_6.RData")
famosParMat_rep_var<-t(apply(parMat,1,rep2orig,boundMat=boundMat))[,-1]

###  3. Hand Tune
load("Official_Fast/output_validation/handTuneResults.RData")
handTuneOutput<-modelRun[[1]]
handTuneParMat<-c(0 , # PCTIM 
               0.1 , # ADIMP 
               -0.7 , # UZTWM
               -0.7 , # LZTWM
               -1.3 , # LZFSM
               -1.2 , # LZFPM
               -2.3 , # LZSK
               1.0 , # snow_SCF
               -1.5 , # REXP
               -1.3 , # UZK
               3.0 , # rutpix_Q0CHN
               1.0) # rutpix_QMCHN



###  4. Pre-calibration
load("Official_Fast/output_validation/precalibration_goodRuns.RData")
# precalibrationOutput<-modelOutput[goodRuns,]
load("Official_Fast/precalibration/output/mhParameters_0.RData")
precalibrationParMat_MSE<-parMat[goodRuns_MSE,-1]
precalibrationParMat_Window<-parMat[goodRuns_Window,-1]


##########################################################################################
##########################################################################################
# Save files (for Iman and Sanjib)
# names(handTuneOutput)<-colnames(precalibrationOutput)<-colnames(emulationOutput)<-colnames(famosOutput)<-names(extremeObs)<-extremeDate
# save(extremeDate , extremeObs , 
#      famosOutput ,  emulationOutput , precalibrationOutput , handTuneOutput,
#      file="Analysis/resultsStremflow_calibration.RData") 

save(famosParMat , 
     famosParMat_rep_var , 
     handTuneParMat , 
     precalibrationParMat_MSE,
     precalibrationParMat_Window,
     file="manuscript/revisionCode/calibrationParameters.RData") 

