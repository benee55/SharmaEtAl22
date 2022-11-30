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

# observation Index
extremeDate<-dateVect[obsInd] # Extreme Dates
extremeObs<-subsetFinalObs # Extreme Values
validationDate<-dateVect[validationInd] # Extreme Dates
validationObs<-subsetFinalValidation # Extreme Values

##########################################################################################
##########################################################################################
# Compile Data
##########################################################################################
##########################################################################################

### 1. Famous Full
load("Official_Fast/output_validation/famosResults_rep_2k.RData")
famosOutput<-matrix(unlist(outputMat[1,]),
                    nrow=length(outputMat[1,]), 
                    ncol=length(outputMat[1,][[1]]), byrow = TRUE)
load("Official_Fast/output_rep/mhParameters_4.RData")
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
famosParMat<-t(apply(parMat,1,rep2orig,boundMat=boundMat))[,-1]
### 2. Famous 4 Parameter
load("~/Dropbox/FamosHydroModel/lowDim/output_f/mhParameters_4.RData")
famos4ParMat<-t(apply(parMat,1,rep2orig,boundMat=boundMat[c(1,2,11,12),]))[,-1]

### 3. Emulation-Calibration
load("lowDim/output/finalEmulationResults.RData")
emulationOutput<-predTest_proj
load("lowDim/output/mcmcResults.RData")
emulationParMat<-testMat

###  4. Hand Tune
load("Official_Fast/output_validation/handTuneResults.RData")
handTuneOutput<-outputMat[[1]]
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

###  5. Pre-calibration
# Compile precalibration data
for(i in 1:14){
  
  if(i==1){
    load(paste("Official_Fast/output_validation/preCalibrationResults",i,".RData",sep=""))
    bar<-unlist(outputMat[1,])
    modelOutput<-matrix(bar, nrow=length(outputMat[1,]), ncol=length(outputMat[1,1][[1]]), byrow = TRUE)
    
  }else{
    load(paste("Official_Fast/output_validation/preCalibrationResults",i,".RData",sep=""))
    bar<-unlist(outputMat[1,])
    foo<-matrix(bar, nrow=length(outputMat[1,]), ncol=length(outputMat[1,1][[1]]), byrow = TRUE)
    modelOutput<-rbind(modelOutput,foo)
  }
}
load("Official_Fast/output_validation/precalibration_goodRuns.RData")
precalibrationOutput<-modelOutput[goodRuns,]
load("Official_Fast/precalibration/output/mhParameters_0.RData")
precalibrationParMat<-parMat[goodRuns,-1]




##########################################################################################
##########################################################################################
# Point Estimation: RMSPE for Validation 
##########################################################################################
##########################################################################################
rmspeTable<-vector("numeric")
### 1. Famous Full
pred<-apply(famosOutput,2,mean)
rmspeTable[1]<-sqrt(mean(pred-validationObs)^2)
### 2. Famous 4 Parameter

### 3. Emulation-Calibration
pred<-apply(emulationOutput,2,mean)
rmspeTable[3]<-sqrt(mean(pred-validationObs)^2)
###  4. Hand Tune
pred<-handTuneOutput
rmspeTable[4]<-sqrt(mean(pred-validationObs)^2)
###  5. Pre-calibration
pred<-apply(precalibrationOutput,2,mean)
rmspeTable[5]<-sqrt(mean(pred-validationObs)^2)

##########################################################################################
##########################################################################################
# Interval Estimation: 95% CI + Coverage 
##########################################################################################
##########################################################################################
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

intervalTable<-vector("numeric")
### 1. Famous Full
rangeMat<-rbind(validationObs, apply(famosOutput, 2, hpd))
intervalTable[1]<-mean(apply(rangeMat, 2, function(x){x[1]<x[3]&x[1]>x[2]}))
### 2. Famous 4 Parameter

### 3. Emulation-Calibration
rangeMat<-rbind(validationObs, apply(emulationOutput, 2, hpd))
intervalTable[3]<-mean(apply(rangeMat, 2, function(x){x[1]<x[3]&x[1]>x[2]}))

###  4. Hand Tune
intervalTable[4]<-NA
###  5. Pre-calibration
rangeMat<-rbind(validationObs, apply(precalibrationOutput, 2, hpd))
intervalTable[5]<-mean(apply(rangeMat, 2, function(x){x[1]<x[3]&x[1]>x[2]}))
##########################################################################################
##########################################################################################
# Combined Table 
##########################################################################################
##########################################################################################
comboTable<-t(rbind(rmspeTable, intervalTable))
colnames(comboTable)<-c("rmspe" , "coverage")
rownames(comboTable)<-c("Famous-12" , "Famous-4" , "Emulation-4" , "HandTune-12" , "Precalibration-12")
library(xtable)
print(xtable(comboTable))

##########################################################################################
##########################################################################################
# Visualization - Validation Dates Densities
##########################################################################################
##########################################################################################
load("Official_Fast/input/fullObservations.RData") # Load Full observation
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

png(file = "Projections.png", height = 750, width=1000)
# Plot Extreme Dates
par(mfrow=c(4,5), mar=c(2,2,2,2))
plot(x=1,y=1,typ="n")
# legend("center" , legend=c("Famous - 12 par","Famous - 4 par" , "Emulation - 4 par" , 
#                            "Precalibration - 12 par" , "Handtune - 12 par" , "Truth"), 
#        lty=c(1,1,2,1,1,1), col=c("black","red","red","gray","green","blue"), 
#        lwd=c(2,2,2,2,2,2) , cex=1.5, bty = "n")
legend("center" , legend=c("Famous - 12 par","Emulation - 4 par" ,
                           "Precalibration - 12 par" , "Handtune - 12 par" , "Truth"),
       lty=c(1,2,1,1,1), col=c("black","red","gray","green","blue"),
       lwd=c(2,2,2,2,2) , cex=1.5, bty = "n")
for(k in 1:ncol(famosOutput)){
  d1<-density(famosOutput[,k])
  # d2<-density(emulationOutput[,k])
  d3<-density(emulationOutput[,k])
  d4<-density(precalibrationOutput[,k])
  plot(d1, main=validationDate[k] , 
       xlim=range(d1$x, d2$x ,d3$x, d4$x , validationObs[k] , handTuneOutput[k]), 
       ylim=range(d1$y, d2$y,d3$y, d4$y))
  # lines(d2, col="red")
  lines(d3, col="red", lty=2)
  lines(d4, col="gray")
  abline(v=validationObs[k] , col="blue")
  abline(v=handTuneOutput[k] , col="green")
}
dev.off()
# Plot Parameters

# Parameter Names
parNames<-c("PCTIM" , "ADIMP" , "UZTWM" ,"LZTWM" , 
            "LZFSM" , "LZFPM" , "LZSK" , "snow_SCF" ,
            "REXP" , "UZK" , "Q0CHN" , "QMCHN")

png(file = "Parameters.png", height = 500, width=750)
par(mfrow=c(4,4), mar=c(2,2,2,2))
plot(x=1,y=1,typ="n")
legend("center" , legend=c("Famous - 12 par","Famous - 4 par" , "Emulation - 4 par" ,
                           "Precalibration - 12 par" , "Handtune - 12 par" , "Prior"),
       lty=c(1,1,2,1,1,1), col=c("black","red","red","gray","green","blue"),
       lwd=c(2,2,2,2,2,2) , cex=1 , bty = "n")
# legend("center" , legend=c("Famous - 12 par" , "Emulation - 4 par" , 
#                            "Precalibration - 12 par" , "Handtune - 12 par" , "Prior"), 
#        lty=c(1,2,1,1,1), col=c("black","red","gray","green","blue"), 
#        lwd=c(2,2,2,2,2) , cex=1 , bty = "n")
for(k in 1:ncol(famosParMat)){
  d1<-density(famosParMat[,k])
  d4<-density(precalibrationParMat[,k])
  if(k %in% c(1,2,11,12)){
    hk<-ifelse(k%in%c(1,2), k, k-8)
    d2<-density(famos4ParMat[,hk])
    d3<-density(emulationParMat[,hk])
    plot(d1, main=parNames[k] , 
         xlim=range(d1$x, d3$x, d4$x , handTuneParMat[k],boundMat[k,]), 
         ylim=range(d1$y, d3$y, d4$y))  
    lines(d2, col="red")
    lines(d3, col="red", lty=2)
  }else{
plot(d1, main=parNames[k] , 
         xlim=range(d1$x, d4$x , handTuneParMat[k],boundMat[k,]), 
         ylim=range(d1$y, d4$y))

  }
  lines(d4, col="gray")
  abline(v=handTuneParMat[k] , col="green")
  abline(v=boundMat[k,] , col="blue")
}

dev.off()


##########################################################################################
##########################################################################################
# Save files (for Iman)
names(handTuneOutput)<-colnames(precalibrationOutput)<-colnames(emulationOutput)<-colnames(famosOutput)<-names(validationObs)<-validationDate
save(validationDate , validationObs , 
     famosOutput ,  emulationOutput , precalibrationOutput , handTuneOutput,
     file="resultsStremflow.RData") 
