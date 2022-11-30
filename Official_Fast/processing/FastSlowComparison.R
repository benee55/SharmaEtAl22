rm(list=ls())
setwd("~/work/LeeEtal-PSUICE-3D/calibration/Official_Fast/output/")
source("../run/mcmc_source_Tr.R")

dev.off()
pdf(file="../processing/FastRegularComparison.pdf",height=8.5,width=11)

load("~/work/LeeEtal-PSUICE-3D/calibration/Official/output/mhParameters_6.RData")
colnames(parMat)<-parNames
parMat<-t(apply(parMat,1,transformLogBase,log=TRUE))
fullparMat<-parMat
load("mhParameters_2.RData")
colnames(parMat)<-parNames
parMat<-t(apply(parMat,1,transformLogBase,log=TRUE))

par(mfrow=c(3,4),mar=c(2,2,2,2))
plot(x=1,y=1,typ="n",main="Posterior Densities - Parameters",cex.main=1.5)
legend("topright",legend=c("Adaptive - 6.5 hours","Standard - 127 hours \n(5.3 days)"),
       col=c(1:2),lty=c(1,1),cex=1.5)
for(i in 1:ncol(parMat)){
  d1<-density(parMat[,i])
  d2<-density(fullparMat[,i])
  plot(d1,col="black",lwd=2,ylim=range(d1$y,d2$y),main=parNames[i],cex.main=1.5)
  lines(d2,col="red",lwd=2)
}



load("~/work/LeeEtal-PSUICE-3D/calibration/Official/output/contribution_6.RData")
fullContribution<-FinalcontributionMat
load("contribution_2.RData")



projNames<-c("Pliocene","LIG","LGM","Modern Volume","Modern Area","2100","2200","2300","2400","2500")
par(mfrow=c(3,4),mar=c(2,2,2,2))
plot(x=1,y=1,typ="n",main="Posterior Densities - Projections",cex.main=1.5)
legend("topright",legend=c("Adaptive - 6.5 hours","Standard - 127 hours \n(5.3 days)"),
       col=c(1:2),lty=c(1,1),cex=1.5)
for(i in 1:ncol(FinalcontributionMat)){
  d1<-density(FinalcontributionMat[,i])
  d2<-density(fullContribution[,i])
  plot(d1,col="black",lwd=2,ylim=range(d1$y,d2$y),main=projNames[i],cex.main=1.5)
  lines(d2,col="red",lwd=2)
}
dev.off()
