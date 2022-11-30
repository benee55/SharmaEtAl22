
# Battacharyya Distance 
countTable<-function(x,dat1,dat2){
  d1<-sum(dat1>=x[1]&dat1<x[2],na.rm = TRUE)/length(dat1)
  d2<-sum(dat2>=x[1]&dat2<x[2],na.rm = TRUE)/length(dat2)
  BC<-(sqrt(d1*d2))
}

# SourceCode for Importance Sampling with Discrepancy
forwardmodel<-function(mu,X){
  return(5*exp(-mu*(X)))
}

forwardmodelDisc<-function(mu,X){
  return(5*exp(-mu*(X))-(1.5*X))
}

sqeCov<-function(distMat,phi,sigma2){
  sigma2*exp(-0.5*(distMat/phi)^2)+diag(nrow(distMat))*0.000001
}
expCov<-function(distMat,phi,sigma2){
  sigma2*exp(-distMat/phi)+diag(nrow(distMat))*0.000001
}
hpd <- function(samp,p=0.05){
  ## to find an approximate (1-p)*100% HPD interval from a
  ## given posterior sample vector samp
  
  r <- length(samp)
  samp <- sort(samp)
  rang <- matrix(0,nrow=trunc(p*r),ncol=3)
  dimnames(rang) <- list(NULL,c("low","high","range"))
  for (i in 1:trunc(p*r)) {
    rang[i,1] <- samp[i]
    rang[i,2] <- samp[i+(1-p)*r]
    rang[i,3] <- rang[i,2]-rang[i,1]
  }
  hpd <- rang[order(rang[,3])[1],]
  return(hpd)
}

# Likelihood
llhd.temper<-function(par,obs,Xvar,distMat,covFun, tempExp){
  compMod<-forwardmodel(mu=par[1],X=Xvar)
  covMat<-covFun(distMat = distMat,phi = par[2],sigma2 = par[3])
  llhd<-dmvnorm(x=as.numeric(obs),mean = compMod,
                sigma = covMat+par[4]*diag(length(obs)),log = TRUE)
  return((tempExp)*llhd)
}

# Prior
lprior.mu<-function(parameter,prior1,prior2){
  lprior<-dnorm(x = parameter,mean = prior1, sd = sqrt(prior2), log = TRUE)
  return(lprior)
}

lprior.phi<-function(parameter,prior1,prior2){
  lprior<-dunif(x = parameter,min = prior1, max = prior2, log = TRUE)
  return(lprior)
}

lprior.sigma2<-function(parameter,prior1,prior2){
  lprior<-dinvgamma(x = parameter,shape = prior1,scale=prior2,log = TRUE)
  return(lprior)
}

lprior.tau2<-function(parameter,prior1,prior2){
  lprior<-dinvgamma(x = parameter,shape = prior1,scale=prior2,log = TRUE)
  return(lprior)
}



logPost<-function(par,obs,Xvar,distMat,covFun, tempExp,prior1,prior2){
  if(par[2]<prior1[2]|par[2]>prior2[2]|par[3]<0|par[4]<0){
    lpost=-Inf
  }else{
    llhd<-llhd.temper(par=par,obs=obs,Xvar=Xvar,distMat=distMat,
                            covFun=covFun,tempExp=tempExp)
    lprior<-sum(lprior.mu(parameter=par[1],prior1=prior1[1],prior2=prior2[1]),
                lprior.phi(parameter=par[2],prior1=prior1[2],prior2=prior2[2]),
                lprior.sigma2(parameter=par[3],prior1=prior1[3],prior2=prior2[3]),
                lprior.tau2(parameter=par[4],prior1=prior1[4],prior2=prior2[4])
    )
    lpost<-llhd+lprior
  }
  return(lpost)
}

################
# Turn Full WEights to tempered weights and ESS
fullToESS<-function(fullWeight,tempSeq,showESS=FALSE){
  weights<-fullWeight*tempSeq
  weights<-exp(weights-max(weights))
  weights<-weights/sum(weights)
  ESS<-1/sum((weights^2))
  
  if(showESS==TRUE){
    return(ESS)
  }else{
    return(list(weights=weights,ESS=ESS))
  }
  }

#Maximize ESS
maxESS<-function(fullWeight,tempSeq,ens){
  ESS<-fullToESS(fullWeight,tempSeq,showESS=TRUE)
  (ESS-0.5*ens)^2
}

#Bhattacharrya Distance
# Count Table Function
countTable<-function(x,dat1,dat2){
  d1<-sum(dat1>=x[1]&dat1<x[2],na.rm = TRUE)/length(dat1)
  d2<-sum(dat2>=x[1]&dat2<x[2],na.rm = TRUE)/length(dat2)
  BC<-(sqrt(d1*d2))
}
# Function to evaluate Battacharrya Distance 
battaDistance<-function(dat.current,dat.previous,breaks){
  
  penEnd<-quantile(dat.previous,probs = c(0.005,0.995),na.rm = TRUE)
  partition<-c(-Inf,seq(penEnd[1],penEnd[2],length.out = breaks-2),Inf)
  partitionTable<-rbind(partition[-length(partition)],partition[-1])
  
  sampBD<--log(sum(apply(partitionTable,2,countTable,
                         dat1=dat.previous,
                         dat2=dat.current)))
  return(sampBD)
}

################
pSpider<-function(ens,initPar,cores,type,obs,Xvar,distMat,covFun,prior1,prior2,
                  tempSched,pfIter,parTruth,mcmcComp,TuningMat){
  
  #Add Bhattacharyya MC sample
  sampBD<-vector("numeric");sampBD[1]<-NA
  #Initialization
  parNames<-c("mu","phi","sigma2","tau2")
  parList<-list()
  parList[[1]]<-initPar
  uniqueParticles<-vector("numeric");uniqueParticles[1]<-ens
  totIter<-length(tempSched)+1
  tempMCMC<-cumsum(tempSched)
  # Initialize Parallelization
  cl <- parallel::makeCluster(cores, type=type)
  doParallel::registerDoParallel(cl)
  
  # Calculate Bhattacharrya Distance via parallelization
  breaks=1000
  empirBattaBS<-vector("numeric")
  print("start MC for Bhattacharrya Distance")
  for(j in 1:1000){
    dat1<-rnorm(ens,mean=1,sd=1)
    dat2<-rnorm(ens,mean=1,sd=1)
    penEnd<-quantile(dat1,probs = c(0.025,0.975),na.rm = TRUE)
    partition<-c(-Inf,seq(penEnd[1],penEnd[2],length.out = breaks-2),Inf)
    partitionTable<-rbind(partition[-length(partition)],partition[-1])
    empirBattaBS[j]<--log(sum(apply(partitionTable,2,countTable,
                                    dat1=dat1,
                                    dat2=dat2)))
  }
  bdMC<-quantile(empirBattaBS,probs = 0.95,na.rm = TRUE)
  countBD<-0
  print(paste("end MC for Bhattacharrya Distance: Threshold=",bdMC))
  
  # Start Calibration
  for(k in 2:totIter){
    print(paste("cycle",k))
    # PArallelization #1
    currentPar<-parList[[k-1]]
    print('Weighting')
    pt<-proc.time()
    weights<-foreach::foreach(ens.member=1:ens,.combine='c',.packages=c('adaptMCMC','mvtnorm'),
                              .noexport = c("mcmcComp","parList")) %dopar% {
                                source("source_IS_Bayarri.R")
                                llhd.temper(par=currentPar[ens.member,],obs=obs,Xvar=Xvar,
                                            distMat=distMat,covFun=covFun,tempExp=tempSched[k-1])
                              }
    
    print(proc.time()-pt)
    # Reweight and resample
    weights<-exp(weights-max(weights))
    weights<-weights/sum(weights)
    reSampleInd<-sample(x=1:ens,size = ens,replace = TRUE,prob = weights)
    parList[[k]]<-parList[[k-1]][reSampleInd,]
    #########################
    #Mutation
    print('Mutation')
    pt<-proc.time()
    currentPar<-parList[[k]]
    parList[[k]] <- foreach::foreach(ens.member=1:ens,.noexport = c("mcmcComp","parList"),
                                     .packages=c('adaptMCMC','mvtnorm','invgamma'),
                                     .combine='rbind') %dopar% {
                                       source("source_IS_Bayarri.R")
                                       MCMC(p=logPost, n=pfIter, init=currentPar[ens.member,], adapt=FALSE, 
                                            scale = TuningMat,list=FALSE,
                                            obs=obs,Xvar=Xvar,distMat=distMat,covFun=covFun, 
                                            tempExp=tempMCMC[k-1],
                                            prior1=prior1,prior2=prior2,showProgressBar = FALSE)[pfIter,]
                                     }
    
    print(proc.time()-pt)
    #########################
    uniqueParticles[k]<-length(unique(parList[[k]][,1]))
    par(mfrow=c(2,2))
    for(i in 1:4){
      plot(density(parList[[k]][,i]),main=paste(parNames[i],k,sep=": Cycle "),col="black")
      lines(density(mcmcComp[,i]),col="red")
      abline(v=parTruth[i],lwd=2,col="blue")
    }
    #Check Stopping Criterion
    #sampBD<-function(sample1,sample2)
    penEnd<-quantile(parList[[k]][,1],probs = c(0.025,0.975),na.rm = TRUE)
    partition<-c(-Inf,seq(penEnd[1],penEnd[2],length.out = breaks-2),Inf)
    partitionTable<-rbind(partition[-length(partition)],partition[-1])
    
    sampBD[k]<--log(sum(apply(partitionTable,2,countTable,
                                         dat1=parList[[k]][,1],
                                         dat2=parList[[k-1]][,1])))
    print(paste("BD=",sampBD[k]))
    if(sampBD[k]<bdMC & tempMCMC[k-1]>=0.9){
      countBD<-countBD+1
    }else{
        countBD<-0
    }
    print(paste("Cycles under threshold=",countBD))
    if(countBD>=2){
      print("end algorithm")
      break
      }
    
    
  }
  
  # Return important Info
  return(list(parList,uniqueParticles))
}



