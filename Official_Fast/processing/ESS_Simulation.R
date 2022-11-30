rm(list=ls())
setwd("~/work/LeeEtal-PSUICE-3D/calibration/Official_Fast/output/")
source("../run/mcmc_source_Tr.R")
load("isPost.RData")
load("rsParameters_1.RData")

pdf(file="../processing/ESS_PSU3DICE.pdf",height=8.5,width=11)
par(mfrow=c(1,2))
for(cycle in 1:2){


n.seq<-100
ESS<-vector("numeric")
for(h in 1:n.seq){
  if(cycle==1){
    seqTemp<-seq(0,1,length.out = n.seq)
    newTemper<-seqTemp[h]  
    llhd<-newTemper*(isPost[,cycle])
    llhd<-exp(llhd-max(llhd))
    llhd<-llhd/sum(llhd)
    
  }else{
    seqTemp<-seq(0.1,1,length.out = n.seq)
    newTemper<-seqTemp[h]  
    llhd<-newTemper*(isPost[,cycle])
    llhd<-exp(llhd-max(llhd))
    llhd<-llhd/sum(llhd)
    }

  
  newWeight<-(llhd)/sum(llhd)
  ESS[h]<-1/sum((newWeight^2))
  
}

plot(x=seqTemp,y=ESS,typ="l",main=paste("Likelihood Power vs. ESS \n cycle",cycle),ylab="ESS",xlab="Likelihood Power",lwd=2,
     ylim=range(ESS[-1],1007/2))
abline(v=c(seqTemp[1]+0.1,seqTemp[length(seqTemp)]),col="blue",lwd=2)
abline(h=1007/2,col="red")
optimPoint<-c(0.1,1)[cycle]
llhd<-optimPoint*(isPost[,cycle])
llhd<-exp(llhd-max(llhd))
llhd<-llhd/sum(llhd)
newWeight<-(llhd)/sum(llhd)
newESS<-1/sum((newWeight^2))
points(x=optimPoint,y=newESS,pch=18,col="orange",cex=2)

}

dev.off()