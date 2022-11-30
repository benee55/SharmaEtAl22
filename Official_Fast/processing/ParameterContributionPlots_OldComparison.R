rm(list=ls())
source("../run/mcmc_source.R")
setwd("../output/")
pdf(file="../processing/SummaryFiles.pdf",height=8.5,width=11)
load("mhParameters_0.RData")
priorNewparMat<-parMat
priorNewparMat[,c(1,2,4,9,11)]<-log(parMat[,c(1,2,4,9,11)],10)
priorNewparMat[,10]<-log(parMat[,c(10)]/0.3,10)

load("mhParameters_1.RData")
# load("~/Box Sync/LeeEtal-PSUICE-3D/calibration/MPITest/discr_output/InitialRun/mhParameters_5.RData")
bigNewparMat<-parMat
bigNewparMat[,c(1,2,4,9,11)]<-log(parMat[,c(1,2,4,9,11)],10)
bigNewparMat[,10]<-log(parMat[,c(10)]/0.3,10)
load("../FirstRunOutput/mhParameters_1.RData")
newparMat<-parMat
newparMat[,c(1,2,4,9,11)]<-log(parMat[,c(1,2,4,9,11)],10)
newparMat[,10]<-log(parMat[,c(10)]/0.3,10)

par(mfrow=c(3,4),mar=c(2,2,2,2))
for(i in 1:11){
  d1<-density(bigNewparMat[,i])
  d2<-density(newparMat[,i])
  d3<-density(priorNewparMat[,i])
  plot(d1,xlim=range(d1$x,d2$x,d3$x),ylim=range(d1$y,d2$y,d3$y),col="blue",main=parNames[i])
  lines(d2,col="red")
  lines(d3,col="black",lty=2)
  abline(v=quantile(bigNewparMat[,i],probs = c(0.025,0.975)),col="blue")
  abline(v=quantile(priorNewparMat[,i],probs = c(0.025,0.975),na.rm=TRUE),col="black",lty=2)
  
  
    # rug(newparMat[,i],ticksize = 0.2,lwd = 0.1,col="red")
}

contributionNames<-c("Plio","LIG","LGM","Mod-Vol","Mod-Area","2100","2200","2300","2400","2500")
load("rsParameters_1.RData")
priorFinalcontributionMat<-FinalcontributionMat
load("contribution_1.RData")
bigFinalcontributionMat<-FinalcontributionMat
load("../FirstRunOutput/contribution_1.RData")
par(mfrow=c(3,4),mar=c(2,2,2,2))
for(i in 1:10){
  d1<-density(bigFinalcontributionMat[,i])
  d2<-density(FinalcontributionMat[,i])
  d3<-density(priorFinalcontributionMat[,i],na.rm=TRUE)
  
  plot(d1,xlim=range(d1$x,d2$x,d3$x),ylim=range(d1$y,d2$y,d3$y),col="blue",main=contributionNames[i])
  lines(d2,col="red")
  lines(d3,col="black",lty=2)
  abline(v=quantile(bigFinalcontributionMat[,i],probs = c(0.025,0.975)),col="blue")
  abline(v=quantile(priorFinalcontributionMat[,i],probs = c(0.025,0.975),na.rm=TRUE),col="black",lty=2)
  # hist(newparMat[samp,i])
  # rug(FinalcontributionMat[,i],ticksize = 0.2,lwd = 0.1,col="red")
}
dev.off()