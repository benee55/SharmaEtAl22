rm(list=ls())
setwd("~/Dropbox/FamosHydroModel/lowDim/")
# setwd("/glade/u/home/sanjib/FamosHydroModel/lowDim/")
load("input/fullObservations.RData")
load("output/GPEmulator_Full.RData")
source("run/mcmc_source_Tr.R")
# Log prior
logPrior<-function(par, boundMat){
  sum(dunif(x=par[-5],min=boundMat[,1],max=boundMat[,2],log = TRUE),
      dinvgamma(x = par[5] , shape = 2, rate =2, log=TRUE))
}

logLlhd<-function(par,boundMat,gpEmulator,Jy,dat){
  output<-predict(gpEmulator[[1]], newData = matrix(par[-5],nrow=1))
  for(i in 2:Jy){
    output<-c(output,
              predict(gpEmulator[[i]], newData = matrix(par[-5],nrow=1)))
  }
  # Model Output
  sum(dnorm(x=dat,mean=output,sd=(par[5]),log=TRUE))
  
}

# Log Posterior
logPost<-function(par,boundMat,gpEmulator,Jy,dat){
  
  if(any(par[-5]<boundMat[,1]) | any(par[-5]>boundMat[,2]) | par[5]<0) {
    lPost<- -Inf
  }else{
    lPrior<-logPrior(par=par , boundMat=boundMat)
    llhd<-logLlhd(par = par,
                  boundMat = boundMat,
                  gpEmulator = gpEmulator,
                  Jy = Jy,
                  dat = dat)
    
    lPost<-lPrior+llhd
  }
  
  return(lPost)
}

# MCMC
library(adaptMCMC)
accept.mcmc = 0.234										# Optimal acceptance rate as # parameters->infinity
#	(Gelman et al, 1996; Roberts et al, 1997)
niter.mcmc =100000										# number of iterations for MCMC
gamma.mcmc = 0.55										# rate of adaptation (between 0.5 and 1, lower is faster adaptation)
burnin = round(niter.mcmc*0.1)				# how much to remove for burn-in
stopadapt.mcmc = round(niter.mcmc*5)

par.init<-c(apply(boundMat,1,mean),1000)
useDat<-subsetFinalObs
# Run MCMC
pt<-proc.time()
amcmc.out = MCMC(p=logPost, n=niter.mcmc, init=par.init, 
                 boundMat = boundMat,
                 gpEmulator = gpEmulator,
                 Jy = length(useDat),dat = useDat,
                 acc.rate=accept.mcmc,
                 gamma=gamma.mcmc, list=TRUE, 
                 n.start=round(0.01*niter.mcmc),
                 adapt=TRUE)
ptFinal<-proc.time()-pt
ptFinal
amcmc.out$acceptance.rate
apply(amcmc.out$samples,2,mean)
par(mfrow=c(5,2),mar=c(2,2,2,2))
burnin=niter.mcmc*0.25
newParNames<-c(parNames,"s2")
for(i in 1:5){
  plot.ts(amcmc.out$samples[-(1:burnin),i], main=newParNames[i])
  # abline(h=theta.ori[20,i],col="red")
  plot(density(amcmc.out$samples[-(1:burnin),i]), main=newParNames[i])
  # abline(v=theta.ori[20,i],col="red")
}
save(amcmc.out,file="output/mcmcResults.RData")