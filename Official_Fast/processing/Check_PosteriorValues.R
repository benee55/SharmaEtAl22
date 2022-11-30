rm(list=ls())
##########################################################################################
setwd("~/Box Sync/LeeEtal-PSUICE-3D/calibration/Official/output/")
source("../run/mcmc_source_Tr.R")
logLlhd_tempered<-function(par,temper,transform=FALSE,job,ResultsList){
  sle<-57.554
  names(par)<-parNames
  parEval<-par
  if(transform==TRUE){ # Tranform if necessary
    parEval[logPars[which(logPars!="ENHANCESHELF")]]<-(1*10)^(parEval[logPars[which(logPars!="ENHANCESHELF")]])
    parEval[logPars[which(logPars=="ENHANCESHELF")]]<- 0.3*10^(parEval[logPars[which(logPars=="ENHANCESHELF")]])
  }else{
    par[logPars[which(logPars!="ENHANCESHELF")]]<-log(par[logPars[which(logPars!="ENHANCESHELF")]],10)
    par[logPars[which(logPars=="ENHANCESHELF")]]<- log((par[logPars[which(logPars=="ENHANCESHELF")]]/0.3),10)
  }
  
  if(any(par[1:11]<boundMat[1,]) | any(par[1:11]>boundMat[3,])) {
    llhd<- -Inf
    results<-NA
  }else{
    results<-ResultsList[[job]]
    # results<-psu3dRun(parEval, runDir)
    
    if(any(results[[1]]==-1)){
      llhd<--Inf
    }else{
      output<-c(sle-results[[2]]["plio", "sle"],
                sle-results[[2]]["lig",  "sle"],
                
                results[[2]]["mod", "sle"]
                - results[[2]]["lgm", "sle"],
                
                results[[2]]["mod", "total"]/1e+15,
                results[[2]]["mod", "grounded"]/1e+12)
      
      # Model Output
      # jacobian<-sum(par[logPars]*log(10)) # jacobian from change of variables. Note how many of the terms cancel out
      llhd<-temper*sum(dtnorm(output[1],mean=windows[1,1]+abs(windows[2,1]-windows[1,1])/2,
                              sd=30,lower=windows[1,1],upper=windows[2,1],log=TRUE),
                       dtnorm(output[2],mean=windows[1,2]+abs(windows[2,2]-windows[1,2])/2,
                              sd=10,lower=windows[1,2],upper=windows[2,2],log=TRUE),
                       dtnorm(output[3],mean=windows[1,3]+abs(windows[2,3]-windows[1,3])/2,
                              sd=20,lower=windows[1,3],upper=windows[2,3],log=TRUE),
                       dtnorm(output[4],mean=windows[1,4]+abs(windows[2,4]-windows[1,4])/2,
                              sd=abs((windows[2,4]-windows[1,4])/4),
                              lower=windows[1,4],upper=windows[2,4],log=TRUE), # assume window contains 95%% so STD is window/4
                       dtnorm(output[5],mean=windows[1,5]+abs(windows[2,5]-windows[1,5])/2,
                              sd=abs((windows[2,5]-windows[1,5])/4),
                              lower=windows[1,5],upper=windows[2,5],log=TRUE))
      
      
      
    }}
  return(list(llhd,results))
}
################################################################################################
load("rsParameters_1.RData")
RsParMat<-parMat;colnames(RsParMat)<-parNames
llhd<-initResultsList[[1]]
Results<-initResultsList[[2]]
llhdTemper<-tempSched[1]
temper<-tempMCMC[1]
RsParMat[,logPars[which(logPars!="ENHANCESHELF")]]<-log(RsParMat[,logPars[which(logPars!="ENHANCESHELF")]],10)
RsParMat[,logPars[which(logPars=="ENHANCESHELF")]]<- log((RsParMat[,logPars[which(logPars=="ENHANCESHELF")]]/0.3),10)

lpostA<-temper*(llhd/llhdTemper)+apply(RsParMat,1,logPrior)# Calculate via simpled calculation
lpostB<-vector("numeric")
for(i in 1:nrow(RsParMat)){
  lpostB[i]<-logLlhd_tempered(par=RsParMat[i,],temper=temper,transform=TRUE,job=i,ResultsList=Results)[[1]]+logPrior(RsParMat[i,])
}

lpostC<-vector("numeric")
# Check if this matches up with MCMC Initial Condition
  for(k in 1:nrow(parMat)){
    print(k)
    load(paste("MCMC_1_1_",k,".RData",sep=""))
    lpostC[k]<-amcmc.out[[4]] [1]   
  }
posteriorMat<-round(cbind(lpostA,lpostB,lpostC),digits=10)
sum(posteriorMat[,1]-posteriorMat[,2])
sum(posteriorMat[,1]-posteriorMat[,3])
sum(posteriorMat[,2]-posteriorMat[,3])

#Results show that particles are being initialized correctly


################################################################
#Check if all parameters are continuing correctly on MCMC restarts
lastFirst<-list()
lpostMat<-matrix(NA,nrow=1007,ncol=2)
for(k in 1:1007){
  lastFirst[[k]]<-matrix(NA,nrow=11,ncol=2)
  print(k)
  load(paste("MCMC_1_1_",k,".RData",sep=""))
  lastFirst[[k]][,1]<-amcmc.out[[1]][45,]   
  lpostMat[k,1]<-amcmc.out[[4]][45]
  load(paste("MCMC_1_2_",k,".RData",sep=""))
  lastFirst[[k]][,2]<-amcmc.out[[1]][1,]   
  lpostMat[k,2]<-amcmc.out[[4]][1]
}

sum(unlist(lapply(lastFirst,function(x){all(round(x[,1])==round(x[,2]))})))
sum(apply(lpostMat,1,function(x) all(round(x[1])==round(x[2]))))
