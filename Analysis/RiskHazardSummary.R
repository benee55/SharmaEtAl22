rm(list=ls())
setwd("~/Dropbox/FamosHydroModel/Analysis/")
load("hazard_risk_calib_results.RData")
# Create Index
emulationInd<-which(risk_output$type=="emulation")
famosInd<-which(risk_output$type=="famos")
precalibrationInd<-which(risk_output$type=="precalibration")

#Risk
emulationRisk<-risk_output[emulationInd,]
famosRisk<-risk_output[famosInd,]
precalibrationRisk<-risk_output[precalibrationInd,]
#Hazard
emulationHazard<-hazard_output[emulationInd,]
famosHazard<-hazard_output[famosInd,]
precalibrationHazard<-hazard_output[precalibrationInd,]

#Dates
dates<-names(validation_hazard_output)

# Hazard
png(file = "HazardComparison.png", 
    height = 600, width=800)
par(mfrow=c(4,5), mar=c(2,2,2,2))
plot.new()
legend("center" , legend=c("Famous - 12 par" , "Emulation - 4 par" ,
                           "Precalibration - 12 par" , "Handtune - 12 par" , "Observation"),
       lty=c(1,1,1,1,1,1), col=c("blue","red","gray","green" , "black"),
       lwd=c(2,2,2,2,2,2) , cex=1 , bty = "n")
for(indVal in 1:18){
  d1<-density(precalibrationHazard[,indVal])
  d2<-density(famosHazard[,indVal])
  d3<-density(emulationHazard[,indVal])
  plot(d1,ylim=range(d1$y, d2$y, d3$y) ,
       xlim=range(d1$x, d2$x, d3$x ,
                  validation_hazard_output[indVal] ,
                  handtuned_hazard_output[indVal]), 
       main=paste("Hazard:",dates[indVal]) , col="gray")
  lines(d2 , col="blue") # Famos
  lines(d3 , col="red") # Emulation
  abline(v=validation_hazard_output[indVal], col="black")
  abline(v=handtuned_hazard_output[indVal], col="green")
}
dev.off()

png(file = "RiskComparison.png", 
    height = 600, width=800)
par(mfrow=c(4,5), mar=c(2,2,2,2))
plot.new()
legend("center" , legend=c("Famous - 12 par" , "Emulation - 4 par" ,
                           "Precalibration - 12 par" , "Handtune - 12 par" , "Observation"),
       lty=c(1,1,1,1,1,1), col=c("blue","red","gray","green" , "black"),
       lwd=c(2,2,2,2,2,2) , cex=1 , bty = "n")

for(indVal in 1:18){
d1<-density(precalibrationRisk[,indVal])
d2<-density(famosRisk[,indVal])
d3<-density(emulationRisk[,indVal])
plot(d1,ylim=range(d1$y, d2$y, d3$y) ,
     xlim=range(d1$x, d2$x, d3$x , validation_risk_output[indVal] , handtuned_risk_output[indVal]), 
     main=paste("Risk:",dates[indVal]) , col="gray")
lines(d2 , col="blue") # Famos
lines(d3 , col="red") # Emulation
abline(v=validation_risk_output[indVal], col="black")
abline(v=handtuned_risk_output[indVal], col="green")
}
dev.off()




rm(list=ls())
load("~/Dropbox/FamosHydroModel/Official_Fast/input/fullObservations.RData") # Load Full observation
dateVect<-paste(sprintf("%04d",as.numeric(obs[,1])),
               sprintf("%02d",as.numeric(obs[,2])),
               sprintf("%02d",as.numeric(obs[,3])),sep="-")
dateVect<-as.Date(dateVect, format = "%Y-%m-%d")
# observation Index
extremeDate<-dateVect[obsInd] # Extreme Dates

load("~/Dropbox/FamosHydroModel/Analysis/resultsStremflow_calibration_full.RData")
dim(famosOutput4Par)


extremeFourPar<-famosOutput4Par[,obsInd]
extremeFamos<-famosOutput[,obsInd]
extremePrecalib<-precalibrationOutput[,obsInd]
extremeHand<-handTuneOutput[obsInd]


par(mfrow=c(5,5), mar=c(2,2,2,2))
plot.new()
legend("center", legend=c("Famos4", "Famos12"  , "precalibration" , "HandTune" , "Truth"), 
       lty=c(1,1,1,1,1), col=c("red","blue","gray","green","black"))
for(k in 1:length(subsetFinalObs)){
  d1<-density(extremeFourPar[,k])
  d2<-density(extremeFamos[,k])
  d3<-density(extremePrecalib[,k])
  
  plot(d1, col="red" , ylim=range(d1$y,d2$y,d3$y), 
       xlim=range(d1$x,d2$x,d3$x,subsetFinalObs[k] ,extremeHand[k]), 
       main=extremeDate[k])
  lines(density(extremeFamos[,k]), col="blue")
  lines(density(extremePrecalib[,k]), col="gray")
  abline(v=subsetFinalObs[k], col="black")
  abline(v=extremeHand[k], col="green")
  
}














