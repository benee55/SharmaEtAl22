
rm(list=ls())
setwd("~/Dropbox/FamosHydroModel/lowDim/")
# Train Data
load("output/modelRuns.RData")
# process output
modelRuns<-matrix(NA,nrow=length(outputMat[1,]), ncol=length(unlist(outputMat[1,1])))
for(k in 1:length(outputMat[1,])){
  modelRuns[k,]<-  unlist(outputMat[1,k])
}
# Test Data
load("output/testModelRuns.RData")
testRuns<-matrix(NA,nrow=length(outputMat[1,]), ncol=length(unlist(outputMat[1,1])))
for(k in 1:length(outputMat[1,])){
  testRuns[k,]<-  unlist(outputMat[1,k])
}

# load Design Matrix
load("input/design.RData")
parMat<-round(parMat,3)


# Try out Emulator
load("output/GPEmulator_Full.RData")
library(mlegp)
# PRedict on test parameters
predTest_CM<-predTest<-matrix(NA, nrow=nrow(testMat) , ncol=21)
for(k in 1:21){
  print(k)
  predTest[,k]<-apply(testMat, 1, predict, object=gpEmulator[[k]])
  predTest_CM[,k]<-apply(testMat, 1, predict, object=gpEmulator_CM[[k]])
}
#Validation
rmse<-sqrt(apply((testRuns-predTest)^2,1,mean))
rmse_CM<-sqrt(apply((testRuns-predTest_CM)^2,1,mean))
mean(rmse)
mean(rmse_CM)


#load observations
load("input/fullObservations.RData")
subsetFinalObs

# Format Date
dateVect<-paste(sprintf("%04d",as.numeric(obs[,1])),
                sprintf("%02d",as.numeric(obs[,2])),
                sprintf("%02d",as.numeric(obs[,3])),sep="-")
dateVect<-as.Date(dateVect, format = "%Y-%m-%d")

# observation Index
extremeDate<-dateVect[obsInd] # Extreme Dates

dev.off()
# sampInd<-sample(1:71, 25)
sampInd<-order(rmse,decreasing = TRUE)[1:11]
par(mfrow=c(3,4),mar=c(5,4,2,2))
plot.new()
legend("center" , legend=c("Truth", "Emulator") , 
       pch=c(16,16) , col=c("black","red"))
for(k in 1:length(sampInd)){
  hj<-sampInd[k]
  plot(x=1:21, y=predTest[hj,],typ="p" , 
       ylim=range(predTest[hj,],testRuns[hj,]), 
       main=paste("Test Run", hj), 
       ylab="Streamflow" , xlab="", xaxt="none")
  axis(1, at=1:21, labels=extremeDate,las=2)
  points(x=1:21, y=predTest[hj,], pch=16, col="black")
  points(x=1:21, y=testRuns[hj,], pch=16, col="red")
}

par(mfrow=c(2,3))
plot.new()
legend("center", legend=c("Truth","Emulator" , "Design Point"), 
       col=c("black","red","black"), 
       lty=c(1,1,NA), pch=c(NA,NA,16))

parNames<-c("PCTIM" , "ADIMP" , "Q0CHN" , "QMCHN")
# First parameter
apply(parMat,2,unique)
testSeq<-seq(0,5, length.out = 100)
resVect_CM<-resVect<-vector("numeric")
for(j in 1:length(testSeq)){
  # print(j)
  jobPar<-c(testSeq[j],1.0,2.5,1.1)
  resVect[j]<-predict(object = gpEmulator[[1]], newData = jobPar)
  resVect_CM[j]<-predict(object = gpEmulator_CM[[1]], newData = jobPar)
}

plot(x=testSeq, y=resVect , typ="l" , ylim=range(resVect), 
     ylab="Streamflow" , xlab=parNames[1] , main=parNames[1]);
lines(x=testSeq , y=resVect_CM, col="red")
foo<-which(round(parMat[,2],1)==1.0 & round(parMat[,3],1)==2.5 & round(parMat[,4],1)==1.1)
points(x=parMat[foo,1] , y=modelRuns[foo,1] , pch=16)


# Second parameter
testSeq<-seq(0,2, length.out = 100)
resVect_CM<-resVect<-vector("numeric")
for(j in 1:length(testSeq)){
  # print(j)
  jobPar<-c(2.5,testSeq[j],2.5,1.1)
  resVect[j]<-predict(object = gpEmulator[[1]], newData = jobPar)
  resVect_CM[j]<-predict(object = gpEmulator_CM[[1]], newData = jobPar)
}
plot(x=testSeq, y=resVect , typ="l" , ylim=range(resVect), 
     ylab="Streamflow" , xlab=parNames[2] , main=parNames[2]);
lines(x=testSeq , y=resVect_CM, col="red")
apply(parMat,2,unique)
foo<-which(round(parMat[,1],1)==2.5 & round(parMat[,3],1)==2.5 & round(parMat[,4],1)==1.1)
points(x=parMat[foo,2] , y=modelRuns[foo,1], pch=16)


# Third parameter
apply(parMat,2,unique)
testSeq<-seq(0.5,4.5, length.out = 100)
resVect_CM<-resVect<-vector("numeric")
for(j in 1:length(testSeq)){
  # print(j)
  jobPar<-c(2.5,1,testSeq[j],1.1)
  resVect[j]<-predict(object = gpEmulator[[1]], newData = jobPar)
  resVect_CM[j]<-predict(object = gpEmulator_CM[[1]], newData = jobPar)
}

plot(x=testSeq, y=resVect , typ="l" , ylim=range(resVect), 
     ylab="Streamflow" , xlab=parNames[3] , main=parNames[3]);
lines(x=testSeq , y=resVect_CM, col="red")
foo<-which(round(parMat[,1],1)==2.5 & round(parMat[,2],1)==1.0 & round(parMat[,4],1)==1.1)
points(x=parMat[foo,3] , y=modelRuns[foo,1], pch=16)

# Fourth parameter
apply(parMat,2,unique)
testSeq<-seq(0.3,1.9, length.out = 100)
resVect_CM<-resVect<-vector("numeric")
for(j in 1:length(testSeq)){
  # print(j)
  jobPar<-c(2.5,1,2.5,testSeq[j])
  resVect[j]<-predict(object = gpEmulator[[1]], newData = jobPar)
  resVect_CM[j]<-predict(object = gpEmulator_CM[[1]], newData = jobPar)
}
plot(x=testSeq, y=resVect , typ="l" , ylim=range(resVect), 
     ylab="Streamflow" , xlab=parNames[4] , main=parNames[4]);
lines(x=testSeq , y=resVect_CM, col="red")
apply(parMat,2,unique)
foo<-which(round(parMat[,1],1)==2.5 & round(parMat[,2],1)==1.0 & round(parMat[,3],1)==2.5)
points(x=parMat[foo,4] , y=modelRuns[foo,1], pch=16)






# Plot
library(vioplot)
par(mfrow=c(5,5), mar=c(2,2,2,2))
for(i in 1:21){
  vioplot(modelRuns[,i], ylim=range(modelRuns[,i],subsetFinalObs[i],4950.55, na.rm = TRUE), 
          main = extremeDate[i])
  points(x=1, y=subsetFinalObs[i], col="red" ,pch=16)
  abline(h=4950.55, col="red" ,lwd=1 , lty=2) # ACtion Stage
}

# Good runs
foo<-apply(modelRuns, 1, function(x){sum(x>4950.55)})
goodRuns<-which(foo==21)

par(mfrow=c(5,5), mar=c(2,2,2,2))
for(i in 1:21){
  d1<-density(modelRuns[goodRuns,i])
  plot(d1, xlim=range(d1$x,subsetFinalObs[i],4950.55, na.rm = TRUE), 
       main = extremeDate[i])
  abline(v=subsetFinalObs[i], col="blue" ,pch=16)
  abline(v=4950.55, col="red" ,lwd=1 , lty=2) # ACtion Stage
}
