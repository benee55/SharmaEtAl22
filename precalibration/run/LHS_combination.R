rm(list=ls())
setwd("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/LHSoutput/output/")


for(i in 1:500){
  load(paste("output",i,".RData",sep=""))
  if(i==1){
    outputMat<-matrix(NA,nrow=nrow(lpostVal[[2]]),ncol=500)
    lpostVect<-vector("numeric")
  }
  outputMat[,i]<-lpostVal[[2]][,2]
  lpostVect[i]<-lpostVal[[1]]
}



save(outputMat,lpostVect, file="LHS_output_Full.RData")

dat<-read.table("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/SBYP1_obs.txt")
obsFull<-dat[(30*365+10):(30*365+10+365-1),]
obs<-as.numeric(obsFull[,5])*0.0283 # Need to convert 
rmseVect<-apply(outputMat,2,function(x){sqrt(mean((x-obs)^2))})

top50Ind<-which(rmseVect<quantile(rmseVect,probs = c(0.1)))
top50Output<-outputMat[,top50Ind] ; dim(top50Output)
top20Ind<-which(rmseVect<quantile(rmseVect,probs = c(0.04)))
top20Output<-outputMat[,top20Ind] ; dim(top20Output)
save(top50Output,top20Output, file="LHS_output_Top.RData")
