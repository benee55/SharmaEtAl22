rm(list=ls())
setwd("../output/")
source("../run/mcmc_source.R")
SummaryMat<-list()
sampID<-1:1007
for(kk in 1:11){

SummaryMat[[kk]]<-matrix(NA,nrow=length(sampID),ncol=49)
# Prior
load("mhParameters_0.RData")
SummaryMat[[kk]][,1]<-parMat[sampID,kk]
# PF  
for(i in 1:length(sampID)){
  j=sampID[i]
load(paste("PF_1_",j,".RData",sep=""))  
  SummaryMat[[kk]][i,2]<-jobPar[kk]
}

# RS Summary
load("rsParameters_1.RData")
SummaryMat[[kk]][,3]<-parMat[sampID,kk]
# MCMC
for(i in 1:length(sampID)){
  j=sampID[i]
  load(paste("MCMC_1_1_",j,".RData",sep=""))  
  SummaryMat[[kk]][i,4:48]<-amcmc.out[[1]][,kk]
}
# MCMC Summary
load("mhParameters_1.RData")
SummaryMat[[kk]][,49]<-parMat[sampID,kk]
}

checkList<-function(x){
  x<-round(x,digits = 10)
  c(all(x[,1]==x[,2]),
    all(x[,3]==x[,4]),
    all(x[,48]==x[,49])
  )
}

lapply(SummaryMat,checkList)
pdf(file="TestParticles.pdf")

cbind(SummaryMat[[1]][,3],SummaryMat[[1]][,4],
      SummaryMat[[1]][,3]==SummaryMat[[1]][,4])
checkList(SummaryMat[[1]])

for(kk in 1:11){
  par(mfrow=c(3,4),mar=c(2,2,2,2))  
  apply(SummaryMat[[kk]],1,function(x){plot.ts(x[4:48],main=kk)})
}
dev.off()

