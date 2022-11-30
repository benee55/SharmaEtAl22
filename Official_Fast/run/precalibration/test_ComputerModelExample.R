rm(list=ls())
####################################################################################################
# Toy Example for Famos with resampling
####################################################################################################
# Load Libraries
####################################################################################################
setwd("~/Dropbox/hydroFamos/run/precalibration/")
library(snow);library(doParallel);library(foreach); library(mvtnorm)
####################################################################################################
# Computer Model from Fangzheng Xie and Yanxun Xu (2019) - Bayesian Projected Calibration of Computer Models (Computer Model 1 Section 5.1)
####################################################################################################
compModel<-function(theta,x){
  7*(sin(2*pi*theta[1]-pi))^2+2*(((2*pi*theta[2]-pi)^2)*sin(2*pi*x-pi))
  # theta[1]+theta[2]*x^2
  # return(5*exp(-theta[1]*(x))-(1.5*theta[2]))
}
####################################################################################################
# Generate Observations
####################################################################################################
set.seed(2021)
n=21
trueTheta<-c(0.1,0.3)
inputX<-seq(0.01,0.99,length.out = n)
obsVar<-0.2
err<-rnorm(n=n,mean=0,sd=sqrt(obsVar))
modelOutput<-compModel(theta=trueTheta , x=inputX) 
obs<-modelOutput+err
plot(x=inputX,y=obs,typ="l")
lines(x=inputX, y=modelOutput)
save(inputX, obs , trueTheta, obsVar, n , compModel,file="observationalData.RData")
####################################################################################################
# Generate ensemble
####################################################################################################
rangeMat<-rbind(c(0,0.3),c(0,0.5))
ens<-5000
parMat<-cbind(runif(ens,rangeMat[1,1],rangeMat[1,2]) ,
              runif(ens,rangeMat[2,1],rangeMat[2,2]) ,
              rinvgamma(ens , shape = 2 , rate = 2))
outputMat<-apply(parMat,1,compModel, x=inputX)
####################################################################################################
# Assess Ensemble and Plot 
####################################################################################################
foo<-rbind(parMat[,3],outputMat)
llhd<-apply(foo,2,function(x,obs){sum(dnorm(x=obs , mean = x[-1] , sd=sqrt(x[1]) , log=TRUE))}, obs=obs)
keepInd<-which(llhd>quantile(llhd, probs = 0.90)) # Keep alpha % of the samples with highest log-likelihood

save(compModel, n , obsVar , trueTheta, inputX, obs, rangeMat, ens,  parMat,
     outputMat , keepInd , llhd, file="PrelimOutput.RData")
# Plot Model Runs
par(mfrow=c(1,1))
plot(x=inputX,y=modelOutput,typ="l" , ylim=range(outputMat))
for(i in 1:ens){
  if(i%in%keepInd){col="blue"}else{col="gray"}
  lines(x=inputX,y=outputMat[,i], lwd=0.2, col=col)
}
lines(x=inputX, y= modelOutput , col="black" , lwd=2)
points(x=inputX,y=obs,pch=16, col="red")
# Plot Parameters
par(mfrow=c(2,2), mar=c(2,2,2,2))
for(i in 1:3){
  d1<-density(parMat[,i]) ; d2<-density(parMat[keepInd,i])
  plot(d1, main = paste("Theta",i), 
       ylim=range(d1$y , d2$y),
       xlim=range(d1$x , d2$x))
  lines(d2, col="blue")
  abline(v=c(trueTheta,obsVar)[i] , lwd=2, col="red")
}
# Plot Parameters - Precalibration
par(mfrow=c(2,2), mar=c(2,2,2,2))
for(i in 1:3){
 d2<-density(parMat[keepInd,i])
  plot(d2, main = paste("Theta",i))

  abline(v=c(trueTheta,obsVar)[i] , lwd=2, col="red")
}
####################################################################################################
# MCMC - Helper Files
####################################################################################################


# Gibbs Sampler for Variance
mcmcCompModel<-function(niter, 
                        initPar, 
                        obs, 
                        x , 
                        priorTheta, 
                        priorVar, 
                        propTheta){
  ###################
  llhdFunc<-function(par , obs , x){
    foo<-compModel(theta=par[1:2] , x=x) 
    llhd<-sum(dnorm(x=obs , mean = foo , sd=sqrt(par[3]) , log = TRUE))
    return(list(llhd,foo))
  }
  
  priorPar<-function(par,prior){
    sum(dunif(x=par[1] , min = prior[1,1] , max = prior[1,2] , log = TRUE),
        dunif(x=par[2] , min = prior[2,1] , max = prior[2,2] , log = TRUE))
  }
  
  lpostTheta<-function(par,prior, obs,x){
    foo<-llhdFunc(par=par , obs=obs , x=x)
    lPost<-priorPar(par=par,prior=prior)+foo[[1]]
    return(list(lPost,foo[[2]]))
    
  }
  ###################
  n<-length(obs)
  k<-length(initPar)
  parMat<-matrix(NA, nrow=niter , ncol=k)
  parMat[1,]<-initPar
  output<-lpostTheta(par=initPar , prior = priorTheta , obs = obs , x=x)[[2]]
  for(i in 2:niter){
    if(i%%1000==0){print(i)}
    parMat[i,]<-parMat[i-1,]
    # Update Theta
    cand<-as.numeric(rmvnorm(n=1 , mean=parMat[i-1,1:2] , sigma = propTheta))
    num<-lpostTheta(par=c(cand, parMat[i-1,3]) , prior = priorTheta , obs = obs , x=x)
    den<-lpostTheta(par=parMat[i-1,1:3] , prior = priorTheta , obs = obs , x=x)
    alpha<-min(num[[1]]-den[[1]], log(1))
    if(alpha>log(runif(1))){parMat[i,1:2]<-cand ; output<-num[[2]]}
    # Update Var
    alphaPar<-priorVar[1]+0.5*n 
    betaPar<-(sum((output-obs)^2)+2*priorVar[2])/2
    parMat[i,3]<-rinvgamma(n=1, shape= alphaPar, rate = betaPar)
    
  }
  return(parMat)
}

## Using Ming-Hui Chen's paper in Journal of Computational and Graphical Stats.
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
  hpd <- rang[order(rang[,3])[1],1:2]
  return(hpd)
}

# Calculate Acceptance Rate for MCMC algorithms
accRateFunc<-function(x){
  accRate<-(length(unique(x))-1)/(length(x)-1)
  return(accRate)
}

####################################################################################################
# MCMC - Run MCMC algorithm
####################################################################################################
# Initial Conditions
niter=5000
initPar<-c(apply(rangeMat,1,mean),1)
priorTheta<-rangeMat
priorVar<-c(0.2,0.2)
propTheta<-diag(2)*0.0001
x = inputX
#Pilot Run
mcmcOutput<-mcmcCompModel(niter=niter, 
              initPar=initPar, 
              obs=obs, 
              x = x, 
              priorTheta=priorTheta, 
              priorVar=priorVar, 
              propTheta=propTheta)

# Generate proposal from rosenthal 
propTheta<-cov(mcmcOutput[,1:2])*((2.38^2)/ncol(mcmcOutput))
niter=50000
mcmcOutput<-mcmcCompModel(niter=niter, 
                          initPar=initPar, 
                          obs=obs, 
                          x = x, 
                          priorTheta=priorTheta, 
                          priorVar=priorVar, 
                          propTheta=propTheta)

# Summary
mcmcOutput<-mcmcOutput[-(1:2000),]
par(mfrow=c(3,2)) # Trace Plots and Posterior Densities
plot.ts(mcmcOutput[,1]) ; plot(density(mcmcOutput[,1]))
plot.ts(mcmcOutput[,2]); plot(density(mcmcOutput[,2]))
plot.ts(mcmcOutput[,3]); plot(density(mcmcOutput[,3]))
# Acceptance Rate, Posterior Mean, and HPD
apply(mcmcOutput,2,accRateFunc)
apply(mcmcOutput,2,mean)
apply(mcmcOutput,2,hpd)

# Posterior Runs - PLOT
posteriorOutputMat<-apply(mcmcOutput[,1:2],1,compModel, x=inputX)
meanPost<-apply(posteriorOutputMat,1,mean)
hpdPost<-apply(posteriorOutputMat,1,hpd)

# Plot Prior Runs
par(mfrow=c(1,2))
plot(x=inputX,y=modelOutput,typ="l" , ylim=range(outputMat),
    xlim=range(inputX),main="All Model Runs")
thinID<-ceiling(seq(1,ncol(outputMat), length.out = 1000))
for(i in 1:length(thinID)){
  if(thinID[i]%in%keepInd){col="blue"}else{col="gray"}
  lines(x=inputX,y=outputMat[,thinID[i]], lwd=0.2, col=col)
}
lines(x=inputX, y= modelOutput , col="blue" , lwd=2)
points(x=inputX,y=obs,pch=16, col="red")
# Plot Posterior Runs
plot(x=inputX,y=posteriorOutputMat[,1],typ="n" , ylim=range(outputMat),
     xlim=range(inputX),main="Posterior Model Runs")
thinID<-seq(1,nrow(mcmcOutput), length.out = 1000)
for(i in 1:length(thinID)){
  lines(x=inputX,y=posteriorOutputMat[,thinID[i]], lwd=0.01, col="blue")
}
lines(x=inputX, y= modelOutput , col="blue" , lwd=2)
lines(x=inputX, y= meanPost , col="black" , lwd=2)
lines(x=inputX, y= hpdPost[1,] , col="black" , lwd=2 , lty=2)
lines(x=inputX, y= hpdPost[2,] , col="black" , lwd=2, lty=2)
points(x=inputX,y=obs,pch=16, col="red")

####################################################################################################
# Save FIles from MCMC Approach
####################################################################################################
save(mcmcOutput, modelOutput , meanPost, hpdPost , file="MCMCOutput.RData")
