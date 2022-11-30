# # Copyright (C) 2018 Ben S. Lee
#email: skl5261@psu.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
####################################################################################################

rm(list=ls())
setwd("~/work/LeeEtal-PSUICE-3D/calibration/AOAS_code/ToyExample/")
source("source_IS.R")

library(fields);library(classInt);library(mvtnorm);library(invgamma)
library(snow);library(doParallel);library(foreach);library(adaptMCMC)

load(file="mcmcGoldStandard.RData")
# Step 1: Draw samples from Prior 
ens=2000
pfIter=10
mp_type = "PSOCK" # PSOCK or MPI
nprocs =2
newParList<-list()
newPar<-list()
rsPar<-list()
currentPar<-list()
tempSeq<-list()
ESS<-list()
cumulTemp<-vector("numeric")
useTemp<-vector("numeric")
useESS<-vector("numeric")


cumulTemp[1]<-0;
mastercumulTemp<-cumulTemp[1]
cycle=1


while(mastercumulTemp<0.99){
  print(paste("Cycle:",cycle))
  print(cumulTemp[cycle])
  if(cycle==1){# Priors
    rsPar[[cycle]]<-cbind(rnorm(ens,mean=prior1[1],sd=sqrt(prior2[1])),
                          runif(ens,min=prior1[2],max=prior2[2]),
                          rinvgamma(ens,shape=prior1[3],scale=prior2[3]),
                          rinvgamma(ens,shape=prior1[4],scale=prior2[4]))
  }else{# Previous cycle
    rsPar[[cycle]]<-newPar[[cycle-1]]
  }
  # Full Results
  print(paste("Weighting: Cycle=",cycle))
  cl <- parallel::makeCluster(nprocs, type=mp_type)
  doParallel::registerDoParallel(cl)
  fullWeight  <- foreach::foreach(ens.member=1:ens,
                                  .packages=c('mvtnorm','invgamma'),
                                  .combine='c') %dopar% {
                                    llhd.temper(par = rsPar[[cycle]][ens.member,],
                                                obs=Z,Xvar=X,distMat=distMat,
                                                covFun=expCov,tempExp=1)
                                  }
  
  
  #Selection Criterion
  optimTemp<-optimize(f = maxESS,fullWeight=fullWeight,
                      ens=ens,interval = c(0.1,1-cumulTemp[cycle]))
  useTemp[cycle]<-optimTemp$minimum;
  useESS[cycle]<-fullToESS(fullWeight = fullWeight,tempSeq = useTemp[cycle],showESS = TRUE)
  useWeights<-fullToESS(fullWeight = fullWeight,tempSeq = useTemp[cycle],showESS = FALSE)$weights
  tempSeq[[cycle]]<-seq(0.1,1-cumulTemp[cycle],length.out = 100)
  ESS[[cycle]]<-sapply(tempSeq[[cycle]],fullToESS,fullWeight = fullWeight,showESS=TRUE)
  
  # Resampling 
  print(paste("Resampling: Cycle=",cycle))
  reSampleInd<-sample(x=1:ens,size = ens,replace = TRUE,prob = useWeights)
  currentPar[[cycle]]<-rsPar[[cycle]][reSampleInd,]
  # Mutation
  
  cl <- parallel::makeCluster(nprocs, type=mp_type)
  doParallel::registerDoParallel(cl)
  cumulTemp[cycle+1]<-useTemp[cycle]+cumulTemp[cycle]
  print(paste("Mutation: Cycle=",cycle))
  print(paste("Cumulative Temp=",round(cumulTemp[cycle+1],digits = 4)))
  newParList[[cycle]] <- foreach::foreach(ens.member=1:ens,
                                          .packages=c('adaptMCMC','mvtnorm','invgamma')) %dopar% {
                                            source("source_IS.R")
                                            MCMC(p=logPost, n=pfIter, init=currentPar[[cycle]][ens.member,], adapt=FALSE, 
                                                 scale = TuningMat,list=FALSE,
                                                 obs=Z,Xvar=X,distMat=distMat,covFun=expCov, 
                                                 tempExp=cumulTemp[cycle+1],
                                                 prior1=prior1,prior2=prior2,showProgressBar = FALSE)
                                          }
  newPar[[cycle]]<-t(sapply(newParList[[cycle]],function(x){x[pfIter,]}))
  
  
  
  cycle<-cycle+1
  mastercumulTemp<-max(cumulTemp)
}
##################################################################################################################
##################################################################################################################
resultsMat<-list()# Process Data
for(k in 1:length(newParList)){
  if(k==1){
    for(i in 1:4){
      resultsMat[[i]]<-t(sapply(newParList[[k]],function(x)x[,i]))
    }
  }else{
    for(i in 1:4){
      resultsMat[[i]]<-cbind(resultsMat[[i]],t(sapply(newParList[[k]],function(x)x[,i])))
    }
  }
}

summaryMat<-list()
for(j in 1:4){
  summaryMat[[j]]<-apply(resultsMat[[j]],2,function(x)c(mean(x),hpd(x)))
}

# Bhattacharrya
#Check Bhattacharrya
# Bhattacharrya Test
MCVal<-rnorm(ens,mean=mean(amcmc.out[[1]][,1]),sd=sd(amcmc.out[[1]][,1]))
bdVect<-vector("numeric")
for(k in 1:1000){
  if(k%%250==0){print(k)}
  bdVect[k]  <-battaDistance(dat.current=rnorm(ens,mean=mean(amcmc.out[[1]][,1]),sd=sd(amcmc.out[[1]][,1])),
                             dat.previous=MCVal,breaks=1000)
}

bdRange<-c(mean(bdVect),hpd(bdVect))
bdRange
battEnd<-vector("numeric")
for(i in 1:ncol(resultsMat[[1]])){
  battEnd[i]<-battaDistance(dat.current=resultsMat[[1]][,i],
                            dat.previous=resultsMat[[1]][,ncol(resultsMat[[1]])],
                            breaks=1000)
}

battaDistance(dat.current=resultsMat[[1]][,5],
              dat.previous=resultsMat[[1]][,10],
              breaks=1000)

battaDistance(dat.current=resultsMat[[1]][,15],
              dat.previous=resultsMat[[1]][,20],
              breaks=1000)

battaDistance(dat.current=resultsMat[[1]][,25],
              dat.previous=resultsMat[[1]][,30],
              breaks=1000)

battaDistance(dat.current=resultsMat[[1]][,35],
              dat.previous=resultsMat[[1]][,40],
              breaks=1000)
##################################################################################################################
##################################################################################################################
# Save Data
save(tempSeq,ESS,cumulTemp,ens,pfIter,useTemp,useESS,newPar,summaryMat,resultsMat,
     MCVal,bdVect,bdRange,battEnd,
     file="UpdateAdaptiveParticle.RData")
##################################################################################################################
##################################################################################################################
# Plot Progression of the Cycles
pdf(file="UpdateadaptiveParticle.pdf",height=8.5,width=11)
par(mfrow=c(4,5),mar=c(2,2,2,2))
for(kk in 1:length(tempSeq)){
  #PLot - Cycles and ESS
  plot(x=tempSeq[[kk]],y=ESS[[kk]],typ="l",ylim=range(ens/2,ESS[[kk]]),
       main="ESS vs. Gamma",ylab="ESS",xlab="Power of Likelihood");
  abline(v=c(0.1,1-cumulTemp[kk]),col="blue",lwd=2);
  abline(h=ens/2,col="red",lwd=2)
  points(x=useTemp[kk],y=useESS[kk],pch=16,cex=2,col="orange")
  #PLot - Densities 
  for(i in 1:4){
    d1<-density(amcmc.out[[1]][,i])
    d2<-density(newPar[[kk]][,i])
    plot(d1,col="blue",ylim=range(d1$y,d2$y),xlim=range(d1$x,d2$x),main=parNames[i])
    lines(d2,col="black")
  }
}
##################################################################################################################


########################################################################
#SUmmary of Densities
par(mfrow=c(2,2),mar=c(2,2,2,2))
colVect<-rep(1:(cycle-1),each=pfIter)
for(j in 1:4){
  plot(x=1:ncol(summaryMat[[j]]),y=summaryMat[[j]][1,],
       typ="n",ylim=range(summaryMat[[j]][1:3,]),
       main=parNames[j])  
  lines(x=1:ncol(summaryMat[[j]]),y=summaryMat[[j]][1,],lty=1,col=colVect)
  lines(x=1:ncol(summaryMat[[j]]),y=summaryMat[[j]][2,],lty=2,col=colVect)
  lines(x=1:ncol(summaryMat[[j]]),y=summaryMat[[j]][3,],lty=2,col=colVect)
  abline(h=hpd(amcmc.out[[1]][,j])[1:2],col="red",lwd=1)
}


par(mfrow=c(1,1))
plot(x=1:length(battEnd),y=battEnd,main="BD-Mu",
     ylim=range(bdRange[1:3],battEnd),lwd=2,typ="n")
for(i in 1:(cycle-1)){
  plotSeq<-((i-1)*pfIter)+(1:pfIter)
  lines(x=plotSeq,y=battEnd[plotSeq],col=i)
  abline(v=plotSeq[pfIter],col=i,lwd=2)
}
abline(h=bdRange[3],lwd=2,lty=2)



dev.off()

