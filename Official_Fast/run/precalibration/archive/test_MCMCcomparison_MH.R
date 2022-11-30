rm(list=ls())
setwd("~/Dropbox/hydroFamos/run/precalibration")
source("test_Source.R")
inputDir<-"~/Dropbox/hydroFamos/run/precalibration/input"
outputDir<-"~/Dropbox/hydroFamos/run/precalibration/output"
load("output/mhParameters_0.RData")

######################################################
#############################################################################
# MCMC 
#############################################################################
mcmcManual_tempered_noTrunc<-function(iter,
                              init,
                              propCov,
                              x,
                              obs, 
                              priorPar,
                              temper,
                              llhdTemper,
                              initResults,
                              parNames,
                              additional=FALSE){
  
  # Initialize Values
  n=length(obs)
  names(init)<-parNames
  # Intialize containers
  alphaMat<-parMat<-candMat<-matrix(NA,nrow=iter,ncol=length(init))
  candMat[1,]<-parMat[1,]<-init ; alphaMat[1,]<-1
  colnames(parMat)<-colnames(candMat)<-colnames(alphaMat)<-parNames
  # Add posterior
  lpost<-vector("numeric");
  if(additional==FALSE){
    lpost[1]<-temper*(initResults[[1]]/llhdTemper)+logPrior(par=init,priorPar=priorPar) # Calculate Log Posterior from previous run  - IS step
  }else{
    lpost[1]<-initResults[[1]] # Use log posterior from previous run
  }
  # Note that we are doing this (MCMC Weight- Cumulative) x (Prev)/likelihoodweight + Log Prior, where Prev=(Log Likelihood *likelihoodweight)
  resultsList<-list();resultsList[[1]]<-initResults[[2]]
  curResultsList<-list();curResultsList[[1]]<-initResults[[2]]
  #Begin MCMC
  for(i in 2:iter){
    # Propose on the log scale from the mutlivariate truncated normal
    candMat[i,]<-c(parMat[i-1,1] , mvtnorm::rmvnorm(n = 1,mean = parMat[i-1,-1],
                                            sigma = propCov))
    
    curResults<-logPosterior_temper(par=as.numeric(candMat[i,]) ,
                                    priorPar = priorPar ,
                                    obs = obs ,
                                    x = x ,
                                    temper=temper)
    ############################################################################################################################
    # account for asymmetric proposal
    # Do some simulations to make sure this working properly
    # boundMat<-priorPar[-1,]
    # adj1<-log(mvtnorm::pmvnorm(upper= c(Inf,Inf),lower=c(-Inf,-Inf),
    #                            mean=as.numeric(parMat[i-1,-1]),sigma=propCov))[1]#previous
    # adj2<-log(mvtnorm::pmvnorm(upper= c(Inf,Inf),lower=c(-Inf,-Inf),
    #                            mean=as.numeric(candMat[i,-1]),sigma=propCov))[1]#current 
    # alpha<-min(curResults[[1]]+adj1-lpost[i-1]-adj2,log(1)) # Asymmetric Proposal. Follow Example. Add and subtract constants. 
     # Asymmetric Proposal. Follow Example. Add and subtract constants. 
    #https://journal.r-project.org/archive/2010-1/RJournal_2010-1_Wilhelm+Manjunath.pdf
    ############################################################################################################################
    
    alpha<-min(curResults[[1]]-lpost[i-1],log(1))
    alphaMat[i,]<-c(1,rep(exp(alpha),2)) # Gibbs sampling for variance and exp(alpha) for the thetas
    
    # Accept/Reject
    if(log(runif(n = 1,min = 0,max = 1))<alpha){ # Accept
      parMat[i,]<-candMat[i,];
      lpost[i]<-curResults[[1]]
      resultsList[[i]]<-curResults[[2]]
    }else{ # REject
      parMat[i,]<-parMat[i-1,]
      lpost[i]<-lpost[i-1]
      resultsList[[i]]<-resultsList[[i-1]]
    }
    
    curResultsList[[i]]<-curResults[[2]]
    
    # Update Gamma Parameter - Gibbs Update
    alphaPar<-priorPar[1,1]+0.5*n 
    betaPar<-(sum((resultsList[[i]]-obs)^2)+2*priorPar[1,2])/2
    parMat[i,1]<-rinvgamma(n=1, shape = alphaPar, rate=betaPar)
    
  }
  return(list(parMat,candMat,alphaMat,lpost,resultsList,curResultsList))
}
#######################################################
llhd_t<-logLikelihood_temper(par =parMat[1,], obs = obs ,  x=inputX, temper = 1)
amcmc.out<-mcmcManual_tempered_noTrunc(iter=50000,
                               init=parMat[1,],
                               propCov=diag(2)*0.00001,
                               obs = obs,
                               priorPar = priorPar,
                               temper=1,
                               x=inputX,
                               llhdTemper=1,
                               initResults=llhd_t,
                               parNames=parNames)

amcmc.out2<-mcmcManual_tempered(iter=50000,
                               init=c(1,0.1,0.2),
                               propCov=diag(2)*0.000001,
                               propSigma2=0.1,
                               obs = obs,
                               priorPar = priorPar,
                               temper=1,
                               x=inputX,
                               llhdTemper=1,
                               initResults=llhd_t,
                               parNames=parNames)
plot.ts(amcmc.out2[[1]][-(1:burnin),1])
burnin<-5000
apply(amcmc.out[[1]][-(1:burnin),],2,summary)
apply(amcmc.out2[[1]][-(1:burnin),],2,summary)
apply(amcmc.out[[3]][-(1:burnin),],2, mean)
apply(amcmc.out2[[3]][-(1:burnin),],2, mean)
load("MCMCOutput.RData")
par(mfrow=c(2,2),mar=c(2,2,2,2))
plot.ts(amcmc.out[[1]][-(1:burnin),2])
lines(amcmc.out2[[1]][-(1:burnin),2],col="red")
lines(mcmcOutput[-(1:burnin),1],col="blue")
plot.ts(amcmc.out[[1]][-(1:burnin),3])
lines(amcmc.out2[[1]][-(1:burnin),3],col="red")
plot.ts(amcmc.out[[1]][-(1:burnin),1])
lines(amcmc.out2[[1]][-(1:burnin),1],col="red")
load("MCMCOutput.RData")
par(mfrow=c(2,2))
plot(density(amcmc.out[[1]][-(1:burnin),2]))
lines(density(amcmc.out2[[1]][-(1:burnin),2]),col="red")
lines(density(mcmcOutput[-(1:burnin),1]),col="blue")

plot(density(amcmc.out[[1]][-(1:burnin),3]))
lines(density(amcmc.out2[[1]][-(1:burnin),3]),col="red")
lines(density(mcmcOutput[-(1:burnin),2]),col="blue")

plot(density(amcmc.out[[1]][-(1:burnin),1]))
lines(density(amcmc.out2[[1]][-(1:burnin),1]),col="red")
lines(density(mcmcOutput[-(1:burnin),3]),col="blue")

