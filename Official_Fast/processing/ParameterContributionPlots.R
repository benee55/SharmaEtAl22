rm(list=ls())
source("../run/mcmc_source.R")
setwd("../output/")
pdf(file="../processing/SummaryFiles.pdf",height=8.5,width=11)
# Priors
load("pollardResults.RData")
TransformParPollard<-transformLogBase(parPollard,log=TRUE)
# First Run
# load("mhParameters_3.RData")
load("mhParameters_6.RData")
NewparMat<-parMat
# load("rsParameters_5.RData")
# ISNewparMat<-parMat
# load("masterTotalParticles.RData")
# NewparMat<-masterTotalParticles[-(1:100000),]
NewparMat[,c(1,2,4,9,11)]<-log(NewparMat[,c(1,2,4,9,11)],10)
NewparMat[,10]<-log(NewparMat[,c(10)]/0.3,10)
ISNewparMat[,c(1,2,4,9,11)]<-log(ISNewparMat[,c(1,2,4,9,11)],10)
ISNewparMat[,10]<-log(ISNewparMat[,c(10)]/0.3,10)

par(mfrow=c(3,4),mar=c(2,2,2,2))
for(i in 1:11){
  x<-seq(boundMat[1,i]-0.1*(boundMat[3,i]-boundMat[1,i]),
         boundMat[3,i]+0.1*(boundMat[3,i]-boundMat[1,i]),length.out = 1000)
  y<-dunif(x=x,min = boundMat[1,i],max = boundMat[3,i])
  d1<-density(NewparMat[,i])
  # d2<-density(ISNewparMat[,i])
  
  plot(x=x,y=y,xlim=range(d1$x,x),ylim=range(d1$y,y),col="red",main=parNames[i],typ="l",lwd=2)
  lines(d1,col="blue")
  # lines(d2,col="black")
  abline(v=TransformParPollard[i],lwd=2,lty=1,col="darkgreen")
  # 
  # hist(NewparMat[,i],freq = FALSE,col = "blue",main=parNames[i],lwd=2,breaks=20,plot = TRUE,
       # xlim=range(d1$x,x),ylim=range(d1$y,y))
  # lines(x=x,y=y,col="red",lwd=3)
  
  abline(v=quantile(NewparMat[,i],probs = c(0.025,0.975)),col="blue",lwd=2,lty=2)
    # rug(NewparMat[,i],ticksize = 0.1,lwd = 0.1,col="blue")
}

# Pollard Parameters
parOutput<-outputToResults(resultsPollard)

contributionNames<-c("Plio","LIG","LGM","Mod-Vol","Mod-Area","2100","2200","2300","2400","2500")
load("rsParameters_1.RData")
priorFinalcontributionMat<-FinalcontributionMat
# load("rsParameters_4.RData")
load("contribution_6.RData")
bigFinalcontributionMat<-FinalcontributionMat
# load("rsParameters_5.RData")
# ISFinalcontributionMat<-FinalcontributionMat[reSampleInd,]
par(mfrow=c(3,4),mar=c(2,2,2,2))
for(i in 1:10){
  d1<-density(bigFinalcontributionMat[,i])
  d2<-density(priorFinalcontributionMat[,i],na.rm=TRUE)
  # d3<-density(ISFinalcontributionMat[,i])
  
  plot(d1,xlim=range(d1$x,d2$x),ylim=range(d1$y,d2$y),
       col="blue",main=contributionNames[i],lwd=2)
  lines(d2,col="red",lwd=2)
  # lines(d3,col="black",lwd=2)
  abline(v=parOutput[i],col="darkgreen",lwd=2)
  abline(v=quantile(bigFinalcontributionMat[,i],probs = c(0.025,0.975)),col="blue",lty=2)
  abline(v=quantile(priorFinalcontributionMat[,i],probs = c(0.025,0.975),na.rm=TRUE),col="red",lty=2)
  # rug(FinalcontributionMat[,i],ticksize = 0.2,lwd = 0.1,col="red")
}
dev.off()



par(mfrow=c(4,3),mar=c(2,2,2,2))
for(i in 1){
  load(paste("mhParameters_",i,".RData",sep=""))
  apply(parMat,2,function(x,col){plot(density(x),col)},col=i)

}

