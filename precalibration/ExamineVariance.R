rm(list=ls())
setwd('/gpfs/group/kzk10/default/private/hydrocalib/SGrove/')
# Observation Data
load("LHSoutput/input/LHS_Inputs.RData")
dat<-read.table("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/SBYP1_obs.txt",stringsAsFactors = FALSE)
obsFull<-dat[(30*365+10):(30*365+10+365-1),]
obs<-as.numeric(obsFull[,5])*0.0283 # Observations

# LHS Input
flowMat<-matrix(NA,nrow=365,ncol=1000)
for(i in 1:1000){
  print(i)
  load(paste("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/LHSoutput/output/output",i,".RData",sep=""))
  flowMat[,i]<-lpostVal[[2]][,2] ;
}

squaredDiff<-apply(flowMat,2,function(x){sum((x-obs)^2)})
alpha<-0.002+(365*0.5)
beta<-(squaredDiff+2*0.002)/2

plot(x=1:365, y=flowMat[,1], typ="l")
lines(x=1:365, y=obs, col="red")


rinvgamma(n = 100, shape = alpha, rate=beta)
