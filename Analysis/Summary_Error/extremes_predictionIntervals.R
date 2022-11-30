rm(list=ls())
setwd("~/Dropbox/FamosHydroModel/")
load("Official_Fast/input/fullObservations.RData") # Load Full observation

load("manuscript/revisionCode/resultsStremflow_calibration_full.RData")
load("manuscript/revisionCode/calibrationParameters.RData")
obsErrVar<-mean(famosParMat[,1])
obsErr<-rnorm(n=nrow(famosParMat), mean=0, sd=sqrt(obsErrVar))
dateVect<-dateVect[obsInd]
famosOutput_pred<-famosOutput[,obsInd]+matrix(rnorm(n=length(famosOutput[,obsInd]), mean=0, sd=sqrt(obsErrVar)), 
                                              nrow=nrow(famosOutput[,obsInd]), ncol=ncol(famosOutput[,obsInd]))
precalibrationOutput_Window_pred<-precalibrationOutput_Window[,obsInd]+matrix(rnorm(n=length(precalibrationOutput_Window[,obsInd]), mean=0, sd=sqrt(obsErrVar)), 
                                                                              nrow=nrow(precalibrationOutput_Window[,obsInd]), ncol=ncol(precalibrationOutput_Window[,obsInd]))
handTuneOutput<-handTuneOutput[obsInd]


## Using Ming-Hui Chen's paper in Journal of Computational and Graphical Stats.
hpd <- function(samp,p=0.05){
  ## to find an approximate (1-p)*100% HPD interval from a
  ## given posterior sample vector samp
  
  r <- length(samp)
  samp <- sort(samp)
  rang <- matrix(0,nrow=trunc(p*r),ncol=3)
  dimnames(rang) <- list(NULL,c("low","high","range"))
  for (i in 1:trunc(p*r)) {
    rang[i,1] <- samp[i]
    rang[i,2] <- samp[i+(1-p)*r]
    rang[i,3] <- rang[i,2]-rang[i,1]
  }
  hpd <- rang[order(rang[,3])[1],1:2]
  return(hpd)
}

par(mfrow=c(5,5), mar=c(4,4,2,1))
plot.new()
legend("center", legend=c("95% Prediction \n Interval","Observation"), 
       lty=c(2,1), col=c("black","blue"), cex=1, bty="n", y.intersp = 0.2)
for(k in 1:21){
  plot(density(famosOutput_pred[,k]), col="black", main=dateVect[k], 
       xlab = "Streamflow")
  abline(v=subsetFinalObs[k], col="blue",lwd=2)
  # abline(v=handTuneOutput[k], col="red",lwd=2, lty=1)
  abline(v=hpd(famosOutput_pred[,k]), col="black",lwd=2, lty=2)
  # lines(density(precalibrationOutput_Window_pred[,k]), col="red")
}



save(dateVect , obsErrVar , obsErr , 
     subsetFinalObs , famosOutput_pred ,famosOutput , 
     precalibrationOutput_Window,precalibrationOutput_Window_pred,
     handTuneOutput,
     file="manuscript/revisionCode/predictionInterval_calibration.RData")

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
load("manuscript/revisionCode/calibrationParameters.RData")
obsErrVar<-mean(famosParMat[,1])

famosOutput_pred<-famosOutput[,keepValidation]+matrix(rnorm(n=length(famosOutput[,keepValidation]), mean=0, sd=sqrt(obsErrVar)), 
                                        nrow=nrow(famosOutput[,keepValidation]), ncol=ncol(famosOutput[,keepValidation]))
precalibrationOutput_Window_pred<-precalibrationOutput_Window[,keepValidation]+matrix(rnorm(n=length(precalibrationOutput_Window[,keepValidation]), mean=0, sd=sqrt(obsErrVar)), 
                                                                              nrow=nrow(precalibrationOutput_Window[,keepValidation]), ncol=ncol(precalibrationOutput_Window[,keepValidation]))

## Using Ming-Hui Chen's paper in Journal of Computational and Graphical Stats.
hpd <- function(samp,p=0.05){
  ## to find an approximate (1-p)*100% HPD interval from a
  ## given posterior sample vector samp
  
  r <- length(samp)
  samp <- sort(samp)
  rang <- matrix(0,nrow=trunc(p*r),ncol=3)
  dimnames(rang) <- list(NULL,c("low","high","range"))
  for (i in 1:trunc(p*r)) {
    rang[i,1] <- samp[i]
    rang[i,2] <- samp[i+(1-p)*r]
    rang[i,3] <- rang[i,2]-rang[i,1]
  }
  hpd <- rang[order(rang[,3])[1],1:2]
  return(hpd)
}

par(mfrow=c(4,5), mar=c(4,4,2,1))
plot.new()
legend("center", legend=c("95% Prediction \n Interval","Observation"), 
       lty=c(2,1), col=c("black","blue"), cex=1, bty="n", y.intersp = 0.2)
for(k in 1:18){
  plot(density(famosOutput_pred[,k]), col="black", main=dateVect[k], 
       xlab = "Streamflow")
  abline(v=subsetFinalValidation[k], col="blue",lwd=2)
  # abline(v=handTuneOutput[k], col="red",lwd=2, lty=1)
  abline(v=hpd(famosOutput_pred[,k]), col="black",lwd=2, lty=2)
  # lines(density(precalibrationOutput_Window_pred[,k]), col="red")
}



save(dateVect , obsErrVar , obsErr , 
     subsetFinalObs , famosOutput_pred ,famosOutput , 
     precalibrationOutput_Window,precalibrationOutput_Window_pred,
     handTuneOutput,
     file="manuscript/revisionCode/predictionInterval_validation.RData")


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
