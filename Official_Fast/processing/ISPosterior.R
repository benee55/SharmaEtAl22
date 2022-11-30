rm(list=ls())
setwd("~/work/LeeEtal-PSUICE-3D/calibration/Official_Fast/output/")
source("../run/mcmc_source_Tr.R")

ISPosterior<-function(cycle){
  fileDirLoad<-list.files(pattern = paste("PF_",cycle,sep=""))
  
  orderNum<-sapply(fileDirLoad,function(x) as.numeric(strsplit(x,split="_|[/.]")[[1]][3])) # Need to order it by numerical JobNum
  fileDirLoad<-fileDirLoad[order(orderNum)]
  weightVect<-vector(("numeric"))
  # Load Files
  for(i in 1:length(fileDirLoad)){
    load(fileDirLoad[i])
   weightVect[i]<-llhd_t[[1]]
  }
  return(weightVect)
}

isPost<-matrix(NA,nrow=1007,ncol=2)
for(k in 1:2){
  isPost[,k]<-ISPosterior(k)  
}
save(isPost,file="isPost.RData")
