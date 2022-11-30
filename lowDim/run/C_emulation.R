rm(list=ls())
# setwd("~/Dropbox/FamosHydroModel/lowDim/")
setwd("/glade/u/home/sanjib/FamosHydroModel/lowDim/")
# Load data
load("output/modelRuns.RData")
# process output
modelRuns<-matrix(NA,nrow=length(outputMat[1,]), ncol=length(unlist(outputMat[1,1])))
for(k in 1:length(outputMat[1,])){
  modelRuns[k,]<-  unlist(outputMat[1,k])
}
# summary(modelRuns)

load("output/testModelRuns.RData")
testRuns<-matrix(NA,nrow=length(outputMat[1,]), ncol=length(unlist(outputMat[1,1])))
for(k in 1:length(outputMat[1,])){
  testRuns[k,]<-  unlist(outputMat[1,k])
}
# summary(testRuns)


# load Design Matrix
load("input/design.RData")
parMat<-round(parMat,3)
# Build Emulator via paralellization
if(FALSE){
# Install on cheyenne by loading gnu module
library(snow);library(snowfall);library(mlegp)
sfInit(parallel=TRUE, cpus=7, type='PSOCK')
# How to do zero-mean GP?
gpEmulator_CM<-mlegp(X=parMat,
                  Z=modelRuns,
                  constantMean = 1,
                  nugget = 0,
                  parallel = TRUE)

save(gpEmulator_CM, modelRuns, parMat, 
     file="output/GPEmulator_Full.RData")  
# How to do zero-mean GP?
gpEmulator<-mlegp(X=parMat,
                  Z=modelRuns,
                  constantMean = 0,
                  nugget = 0,
                  parallel = TRUE)
sfStop()

save(gpEmulator,gpEmulator_CM, modelRuns, parMat, 
     file="output/GPEmulator_Full.RData")  
}


# Try out Emulator
load("output/GPEmulator_Full.RData")

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


dev.off()
# sampInd<-sample(1:71, 25)
sampInd<-order(rmse,decreasing = TRUE)[1:25]
par(mfrow=c(5,5),mar=c(2,2,2,2))
for(k in 1:length(sampInd)){
  hj<-sampInd[k]
plot(x=1:21, y=predTest[hj,],typ="l" , 
     ylim=range(predTest[hj,],testRuns[hj,]), 
     main=paste("Test Run", hj))
# lines(x=1:21, y=predTest_CM[1,], col="red")
points(x=1:21, y=predTest[hj,], pch=16, col="black")
# points(x=1:21, y=predTest_CM[1,], pch=16, col="red")
lines(x=1:21, y=testRuns[hj,], col="blue")
points(x=1:21, y=testRuns[hj,], pch=16, col="blue")
# if(k==1){
#   legend("topright", legend=c("Truth","Emulation"), lty=c(1,1) , 
#          col=c("blue","black"))
# }

}

apply(parMat,2,unique)
keepInd<-which(parMat[,3]==2.5 &parMat[,4]==1.1)
bar<-parMat[keepInd,1:2]
baz<-modelRuns[keepInd,1]
source("~/Dropbox/BenJaewoo/SMB_SGLMM/RevisionMaterials/Code/FullModis/source/sharedFunctions.R")
library(fields); library(classInt)
plotRF(dat=baz, location = bar)



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
par(mfrow=c(1,1))
plot(x=testSeq, y=resVect , typ="l" , ylim=range(resVect));
lines(x=testSeq , y=resVect_CM, col="red")
apply(parMat,2,unique)
foo<-which(round(parMat[,2],1)==1.0 & round(parMat[,3],1)==2.5 & round(parMat[,4],1)==1.1)
points(x=parMat[foo,1] , y=modelRuns[foo,1])

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

# Plot
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
