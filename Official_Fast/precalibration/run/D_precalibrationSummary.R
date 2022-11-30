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

############################################################################################################
############################################################################################################
# Failing runs
############################################################################################################
############################################################################################################
foo<-apply(modelOutput,2,function(x)sum(is.na(x)))
badRuns<-which(foo!=0)
badParMat<-parMat[badRuns,]
goodRuns<-which(foo==0)
goodParMat<-parMat[goodRuns,]
failCol<-ifelse(foo==0,"blue","red")
# Plot failing runs with respect to ADIMP and PCTIM
parNames<-c("PCTIM" , "ADIMP" , "UZTWM" ,"LZTWM" , 
            "LZFSM" , "LZFPM" , "LZSK" , "snow_SCF" ,
            "REXP" , "UZK" , "Q0CHN" , "QMCHN")
parA<-"QMCHN"
parB<-"PCTIM"
plot(x=parMat[-badRuns,parA] , y= parMat[-badRuns,parB] , pch=16, cex=0.5,
     xlab=parA , ylab=parB , col="black")
points(x=badParMat[,parA] , y= badParMat[,parB] , pch=16, col="red")
# lines(x=c(0.5, 4.5) , y=c(0.3,0.3), col="blue", lwd=3)
# lines(x=c(0.5, 4.5) , y=c(2.25,2.25), col="blue", lwd=3)
# lines(x=c(0, 0) , y=c(0,0.5), col="blue", lwd=3)
# lines(x=c(0.3, 0.3) , y=c(0,0.5), col="blue", lwd=3)
legend("topright" , 
       legend=c("Success" , "Failed" , "Original Prior"),
       lty=c(NA,NA,1), 
       pch=c(16,16,NA), 
       col=c("black","red","blue"),
       lwd=c(NA,NA,2),cex=1)

summary(badParMat[,-1])
############################################################################################################
# Violin Plots of runs where "QMCHN" is less than 1.9
############################################################################################################
keepNew<-which(goodParMat[,"QMCHN"]<1.9)
par(mfrow=c(5,5), mar=c(2,2,2,2))
for(i in 1:21){
  k<-obsInd[i]
  vioplot(modelOutput[k,keepNew], ylim=range(modelOutput[k,keepNew],extremeObs[i],4950.55, na.rm = TRUE), 
          main = extremeDate[i])
  points(x=1, y=extremeObs[i], col="red" ,pch=16)
  abline(h=4950.55, col="red" ,lwd=1 , lty=2) # ACtion Stage
}
############################################################################################################
# Model Runs where :
# 1. QMCHN<1.9
# 2. Model Runs completed
# 3. Model outputs are all higher than the observations
############################################################################################################

goodParMat<-goodParMat[,-1]
bar<-modelOutput[obsInd,goodRuns]
sdVect<-apply(modelOutput[obsInd,goodRuns],1,sd)
sdVect<-0
aboveObs<-apply(bar,2, function(x){sum(x>=(extremeObs-sdVect))})
extremeIndex<-which(aboveObs==21)
extremeParMat<-(goodParMat[extremeIndex,])
extremeOutput<-(bar[,extremeIndex])

par(mfrow=c(5,5), mar=c(2,2,2,2))
for(k in 1:21){
  vioplot(extremeOutput[k,], ylim=range(extremeOutput[k,],extremeObs[k],4950.55, na.rm = TRUE), 
          main = extremeDate[k])
  points(x=1, y=extremeObs[k], col="red" ,pch=16)
  abline(h=4950.55, col="red" ,lwd=1 , lty=2) # ACtion Stage
}
names(extremeObs)<-rownames(extremeOutput)<-extremeDate
colnames(extremeOutput)<-paste("ModelRun",1:ncol(extremeOutput),sep="")
save(extremeParMat,extremeOutput, extremeObs, extremeDate,  file="captureExtremes.RData")

which(parMat[,2]==extremeParMat[1,1])
############################################################################################################
############################################################################################################
# All runs Violin Plots and Paramter densities
############################################################################################################
par(mfrow=c(5,5), mar=c(2,2,2,2))
for(i in 1:21){
  k<-obsInd[i]
  vioplot(modelOutput[k,], ylim=range(modelOutput[k,],extremeObs[i],4950.55, na.rm=TRUE), 
          main = extremeDate[i])
  points(x=1, y=extremeObs[i], col="red" ,pch=16)
  abline(h=4950.55, col="red" ,lwd=1 , lty=2) # ACtion Stage
}

############################################################################################################
############################################################################################################
# Scoring
############################################################################################################
############################################################################################################
extremeModelOuput<-modelOutput[obsInd,]
MSE<-apply(extremeModelOuput,2,function(x){mean((x-extremeObs)^2)})
goodRuns<-which(MSE<quantile(MSE, probs=0.05, na.rm = TRUE))
goodModelOutput<-modelOutput[modelStart:modelEnd,goodRuns] # Trimmed Data
save(goodRuns, file="~/Dropbox/FamosHydroModel/Official_Fast/output_validation/precalibration_goodRuns.RData")
############################################################################################################
############################################################################################################
# Figures
############################################################################################################
############################################################################################################
# All runs with good runs in green
############################################################################################################
par(mfrow=c(1,1), mar=c(5,4,2,2))
plot(x=newDateVect, y= obs[modelStart:modelEnd,4], typ="n", 
     ylim=range(newModelOutput,na.rm = TRUE), xlim=c(as.Date("2004-08-01") , as.Date("2008-03-31")),
     ylab="Streamflow" , xlab="Date " , 
     main="Streamflow")
for(k in 1:ncol(newModelOutput)){
  lines(x=newDateVect, y= newModelOutput[,k] , col="gray" , lwd=0.5)
}
for(k in goodRuns){
  lines(x=newDateVect, y= newModelOutput[,k] , col="green" , lwd=0.5)
}
lines(x=newDateVect, y= obs[modelStart:modelEnd,4] , col="blue", lwd=1)
points(x=extremeDate , y = extremeObs , col="red" , pch=16, cex=1.5)
abline(h=4950.55, col="red", lty=2)
legend("topright" , legend=c("Observations" , "Model Output" , "Good Model Runs","Extreme Points" , "Action Stage"),
       lty=c(1,1,1,NA,2) , pch=c(NA,NA,NA,16,NA), col=c("blue","gray","green","red","red"),
       lwd=rep(2,2,2,NA,1),cex=0.75)
############################################################################################################
############################################################################################################
# Good runs Violin Plots
############################################################################################################
par(mfrow=c(5,5), mar=c(2,2,2,2))
for(i in 1:21){
  k<-obsInd[i]
  vioplot(modelOutput[k,goodRuns], ylim=range(modelOutput[k,goodRuns],extremeObs[i],4950.55), 
          main = extremeDate[i])
  points(x=1, y=extremeObs[i], col="red" ,pch=16)
  abline(h=4950.55, col="red" ,lwd=1 , lty=2) # ACtion Stage
}
############################################################################################################
############################################################################################################
# Good runs Parameters
############################################################################################################
par(mfrow=c(4,3), mar=c(2,2,2,2))
parNames<-c("PCTIM" , "ADIMP" , "UZTWM" ,"LZTWM" , 
            "LZFSM" , "LZFPM" , "LZSK" , "snow_SCF" ,
            "REXP" , "UZK" , "Q0CHN" , "QMCHN")
boundMat<-rbind(c(0, 5) , # PCTIM 0.3=original maximum
                c(0 , 2), # ADIMP 0.5=original maximum
                c(-50 , -0.1), # UZTWM
                c(-70 , -0.1), # LZTWM
                c(-100 , -0.1), # LZFSM
                c(-120 , -0.1), # LZFPM
                c(-3.8 , -0.1), # LZSK
                c(0.5 , 1.5), # snow_SCF
                c(-3.5 , -0.1), # REXP
                c(-3.5 , -0.1), # UZK
                c(0.5,4.5), # rutpix_Q0CHN
                c(0.3,2.25)) # rutpix_QMCHN Use 2.25 instead of 3.4
for(i in 1:11){
  plot(density(goodParMat[,i]) , xlim=range(boundMat[i,1:2]),  main=parNames[i])
  abline(v=boundMat[i,1:2], col="red")
}


############################################################################################################
############################################################################################################

# Figure - Streamflow 
par(mfrow=c(1,1), mar=c(5,4,2,2))
plot(x=newDateVect, y= obs[modelStart:modelEnd,4], typ="n", 
     ylim=range(newModelOutput, na.rm = TRUE), xlim=c(as.Date("2004-08-01") , as.Date("2008-03-31")),
     ylab="Streamflow" , xlab="Date " , 
     main="Pre-calibration Streamflow")
for(k in 1:ncol(newModelOutput)){
  lines(x=newDateVect, y= newModelOutput[,k] , col="gray" , lwd=0.5)
}
lines(x=newDateVect, y= obs[modelStart:modelEnd,4] , col="blue", lwd=1)
points(x=extremeDate , y = extremeObs , col="red" , pch=16, cex=1.5)
abline(h=4950.55, col="red", lty=2)
legend("topright" , legend=c("Observations" , "Model Output" , "Extreme Points" , "Action Stage"),
       lty=c(1,1,NA,2) , pch=c(NA,NA,16,NA), col=c("blue","gray","red","red"),
       lwd=rep(2,2,NA,1),cex=0.75)

# Figure - Violin Plots of Extreme Dates 
par(mfrow=c(5,5), mar=c(2,2,2,2))
for(i in 1:21){
  k<-obsInd[i]
  vioplot(modelOutput[k,], ylim=range(modelOutput[k,],extremeObs[i],4950.55, na.rm = TRUE), 
          main = extremeDate[i])
  points(x=1, y=extremeObs[i], col="red" ,pch=16)
  abline(h=4950.55, col="red" ,lwd=1 , lty=2) # ACtion Stage
}


# Figure - 2004-2005 of observations vs. max value of model runs
maxVals<-apply(modelOutput,1,max, na.rm = TRUE)
plot(x=dateVect, y=maxVals, lwd=0.5,
     typ="l",xlim=c(as.Date("2004-09-01") , as.Date("2005-04-30")),
     ylim=c(0,20000))
points(x=extremeDate , y = extremeObs , col="red" , pch=16, cex=1.5)
points(x=dateVect , y = maxVals , col="black" , pch=16, cex=1)



# Figure - Streamflow + Includes Spinup Time
useModelOutput<-modelOutput[,keepInd]
par(mfrow=c(1,1), mar=c(5,4,2,2))
plot(x=dateVect, y= obs[,4], typ="n", 
     ylim=range(newModelOutput, na.rm = TRUE), 
     ylab="Streamflow" , xlab="Date " , 
     main="Streamflow with Spinup")
for(k in 1:ncol(useModelOutput)){
  lines(x=dateVect, y= useModelOutput[,k] , col="gray" , lwd=0.5)
}
lines(x=dateVect, y= obs[,4] , col="blue", lwd=1)
points(x=extremeDate , y = extremeObs , col="red" , pch=16, cex=1.5)
abline(h=4950.55, col="red", lty=2)
legend("topright" , legend=c("Observations" , "Model Output" , "Extreme Points" , "Action Stage"),
       lty=c(1,1,NA,2) , pch=c(NA,NA,16,NA), col=c("blue","gray","red","red"),
       lwd=rep(2,2,NA,1),cex=0.75)



# Figure - Streamflow of 2005 observations and model runs
par(mfrow=c(1,1), mar=c(5,4,2,2))
plot(x=newDateVect, y= obs[modelStart:modelEnd,4], typ="n", 
     ylim=c(0,15000), xlim=c(as.Date("2005-01-01") , as.Date("2005-04-15")),
     ylab="Streamflow" , xlab="Date " , 
     main="2005 Streamflow")
for(k in 1:ncol(newModelOutput)){
  lines(x=newDateVect, y= newModelOutput[,k] , col="gray" , lwd=0.5)
}
lines(x=newDateVect, y= obs[modelStart:modelEnd,4] , col="blue", lwd=1)
points(x=extremeDate , y = extremeObs , col="red" , pch=16, cex=1.5)
abline(h=4950.55, col="red", lty=2)
legend("topright" , legend=c("Observations" , "Model Output" , "Extreme Points" , "Action Stage"),
       lty=c(1,1,NA,2) , pch=c(NA,NA,16,NA), col=c("blue","gray","red","red"),
       lwd=rep(2,2,NA,1),cex=0.75)

