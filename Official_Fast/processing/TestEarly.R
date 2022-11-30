rm(list=ls())
setwd("~/work/LeeEtal-PSUICE-3D/calibration/Official_Fast/output/")
load("mhParameters_0.RData")
masterParMat<-list()
masterParMat[[1]]<-parMat
masterParMat[[2]]<-matrix(NA,nrow=2015,ncol=11)
for(i in 1:2015){
  load(paste("PF_1_",i,".RData",sep=""))
  masterParMat[[2]][i,]<-jobPar 
}
load("rsParameters_1.RData")
masterParMat[[3]]<-parMat
masterParMat[[4]]<-matrix(NA,nrow=2015,ncol=11)
masterParMat[[5]]<-matrix(NA,nrow=2015,ncol=11)
for(i in 1:2015){
  load(paste("MCMC_1_1_",i,".RData",sep=""))
  masterParMat[[4]][i,]<-amcmc.out[[1]][1,]
  masterParMat[[5]][i,]<-amcmc.out[[1]][15,]
}
load("mhParameters_1.RData")
masterParMat[[6]]<-parMat
#Cycle 2
masterParMat[[7]]<-matrix(NA,nrow=2015,ncol=11)
for(i in 1:2015){
  load(paste("PF_2_",i,".RData",sep=""))
  masterParMat[[7]][i,]<-jobPar 
}
load("rsParameters_2.RData")
masterParMat[[8]]<-parMat
masterParMat[[9]]<-matrix(NA,nrow=2015,ncol=11)
masterParMat[[10]]<-matrix(NA,nrow=2015,ncol=11)
for(i in 1:2015){
  load(paste("MCMC_2_1_",i,".RData",sep=""))
  masterParMat[[9]][i,]<-amcmc.out[[1]][1,]
  masterParMat[[10]][i,]<-amcmc.out[[1]][15,]
}
load("mhParameters_2.RData")
masterParMat[[11]]<-parMat


foo<-cbind(masterParMat[[1]][,1],
           masterParMat[[2]][,1],
           masterParMat[[3]][,1],
           masterParMat[[4]][,1],
           masterParMat[[5]][,1],
           masterParMat[[6]][,1],
           masterParMat[[7]][,1],
           masterParMat[[8]][,1],
           masterParMat[[9]][,1],
           masterParMat[[10]][,1],
           masterParMat[[11]][,1])
all(round(foo[,1],digits = 4)==round(foo[,2],digits = 4)) #Prior vs. Resampled Init
all(round(foo[,2],digits = 4)==round(foo[,3],digits = 4)) # resmpled vs. MH start
all(round(foo[,3],digits = 4)==round(foo[,4],digits = 4)) # MH Start vs. MH finish
all(round(foo[,4],digits = 4)==round(foo[,5],digits = 4)) 
all(round(foo[,5],digits = 4)==round(foo[,6],digits = 4))
all(round(foo[,6],digits = 4)==round(foo[,7],digits = 4))
all(round(foo[,7],digits = 4)==round(foo[,8],digits = 4))
all(round(foo[,8],digits = 4)==round(foo[,9],digits = 4))
all(round(foo[,9],digits = 4)==round(foo[,10],digits = 4))
all(round(foo[,10],digits = 4)==round(foo[,11],digits = 4))
save(masterParMat,foo,file="CheckTest.RData")

