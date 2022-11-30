
########################################################################################
rm(list=ls())
setwd("~/Dropbox/FamosHydroModel/")
load("Official_Fast/input/fullObservations.RData") # Load Full observation
load("Official_Fast/input/fullObservations.RData") # Load Full observation
# Format Date
dateVect<-paste(sprintf("%04d",as.numeric(obs[,1])),
                sprintf("%02d",as.numeric(obs[,2])),
                sprintf("%02d",as.numeric(obs[,3])),sep="-")
dateVect<-as.Date(dateVect, format = "%Y-%m-%d")

dateVect<-dateVect[which(dateVect=="2003-06-01"):which(dateVect=="2008-03-31")]

# Format Date
# Precalibration Data
# Compile precalibration data
timeVect<-vector("numeric")
for(i in 1:14){
  load(paste("Official_Fast/precalibration/output/preCalibrationResults",i,".RData",sep=""))
  if(i==1){
    modelOutput<-outputMat[-nrow(outputMat),]
    timeVect<-outputMat[nrow(outputMat),]
  }else{
    modelOutput<-cbind(modelOutput,outputMat[-nrow(outputMat),])
    timeVect<-c(timeVect,outputMat[nrow(outputMat),])
  }
}

summary(timeVect/60)


startInd<-which(obs[,1]==2003 & obs[,2]==6 & obs[,3]==1)
endInd<-which(obs[,1]==2008 & obs[,2]==3 & obs[,3]==31)
obs<-obs[startInd:endInd,]
summary(obs[,4])


extremeModelOuput<-modelOutput[obsInd,]
# lowRange<-4950.55 # ACtivation
# lowRange<-0.75*min(obs[,4]) # ACtivation
lowRange<-0.50*4950.55 # ACtivation
highRange<-1.5*max(obs[,4])

lowRange<-0.25*(subsetFinalObs)
highRange<-1.75*(subsetFinalObs)

windowResult<-matrix(NA,nrow(extremeModelOuput), ncol(extremeModelOuput))
for(k in 1:nrow(extremeModelOuput)){
  windowResult[k,]<-extremeModelOuput[k,]>lowRange[k] & extremeModelOuput[k,]<highRange[k]
}
preCalibConditions<-apply(windowResult,2, sum)
table(preCalibConditions)
goodRuns_Window<-which(preCalibConditions==21)


MSE<-apply(extremeModelOuput,2,function(x){mean((x-subsetFinalObs)^2)})
goodRuns_MSE<-which(MSE<quantile(MSE, probs=0.05, na.rm = TRUE))

save(goodRuns_Window, goodRuns_MSE, 
     file="~/Dropbox/FamosHydroModel/Official_Fast/output_validation/precalibration_goodRuns.RData")

extremeDate<-dateVect[obsInd]
par(mfrow=c(5,5), mar=c(2,2,2,2))
for(k in 1:21){
  plot(density(extremeModelOuput[k,]), col="gray" , main=extremeDate[k])
  abline(v=subsetFinalObs[k], col="blue")
  abline(v=lowRange[k], col="red")
  abline(v=highRange[k], col="red")
}


par(mfrow=c(5,5), mar=c(2,2,2,2))
for(k in 1:21){
  d1<-density(extremeModelOuput[k,goodRuns_Window])
  plot(d1, xlim=range(d1$x, lowRange, highRange), col="gray" , main=extremeDate[k])
  abline(v=subsetFinalObs[k], col="blue")
  abline(v=lowRange[k], col="red")
  abline(v=highRange[k], col="red")
}

bestGuess<-apply(extremeModelOuput[,goodRuns_Window],1,mean)


library(hydroGOF)
KGE(sim=bestGuess, obs = subsetFinalObs)
NSE(sim=bestGuess, obs = subsetFinalObs)
cor(bestGuess, subsetFinalObs)
pbias(sim=bestGuess, obs = subsetFinalObs)

