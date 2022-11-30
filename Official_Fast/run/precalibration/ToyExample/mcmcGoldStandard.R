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
set.seed(12345)
n=300
mu=1.7;phi=0.2;sigma2=1;tau2<-0.5
parTruth=c(mu,NA,NA,tau2)
grid.location<-cbind(runif(n,0,1),runif(n,0,1))
distMat <- as.matrix(dist(grid.location))
X<-grid.location[,1]+grid.location[,2]

# Computer Model
ytrue<-forwardmodel(mu=mu,X = X) # Truth
y<-forwardmodelDisc(mu=mu,X = X) # Discrepancy
#Observation error
eps<-rnorm(n,mean=0,sd=sqrt(tau2))
# Discrepancy
W<-y-ytrue
# Observation
Z<-y+eps

# Plot
pdf(file=paste("Discrepancy_n",n,".pdf",sep="'"))
par(mfrow=c(2,2))
breaks <- seq(range(ytrue,y,Z)[1],range(ytrue,y,Z)[2],length.out=20)
pal <- tim.colors(length(breaks)-1)
plotNames<-c("Model Output","Discrepancy", "Model Output + Discrepancy", "Observation")
plotdat<-rbind(ytrue,W,y,Z)
for(i in 1:4){
  fb <- classIntervals(plotdat[i,], n = length(pal),
                       style = "fixed", fixedBreaks = breaks)
  col <- findColours(fb, pal)
  plot(x=grid.location[,1],y=grid.location[,2],col=col,pch=16,cex=1.5,
       main=plotNames[i])
}

dev.off()
####################################################################
#MCMC calibration
# Priors and Initial COnditions
prior1=c(0,0.01,2,2)
prior2=c(100,1.5,2,2)
par.init<-c(1,1,1,1)
covFun<-expCov # Covariance function to model the discrepancy. We choose exponential for a rougher surface

library(adaptMCMC)
accept.mcmc = 0.234										# Optimal acceptance rate as # parameters->infinity
#	(Gelman et al, 1996; Roberts et al, 1997)
niter.mcmc = 200000										# number of iterations for MCMC
gamma.mcmc = 0.7											# rate of adaptation (between 0.5 and 1, lower is faster adaptation)
burnin = round(niter.mcmc*0.5)				# how much to remove for burn-in
stopadapt.mcmc = round(niter.mcmc*1)


# Log LIkelihood
amcmc.out = MCMC(p=logPost, n=niter.mcmc, init=par.init, adapt=TRUE, 
                 acc.rate=accept.mcmc,
                 gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
                 obs=Z,Xvar=X,distMat=distMat,covFun=expCov, tempExp=1,
                 prior1=prior1,prior2=prior2)

amcmc.out$acceptance.rate
chain1<-(amcmc.out$samples)
Fullresults<-apply(chain1,2,function(x)c(mean(x),hpd(x)))

colnames(Fullresults)<-parNames<-c("mu","phi","sigma2","tau2")
rownames(Fullresults)<-c("mean","low95CI","high95CI","rangeCI")
Fullresults
par(mfrow=c(2,2))

for(i in 1:ncol(Fullresults)){
  plot(density(chain1[,i]),main=parNames[i])
  abline(v=Fullresults[2:3,i],col="blue",lty=2)
  abline(v=parTruth[i],col="red",lwd=2)
}

TuningMat<-amcmc.out$cov.jump
save.image(file="mcmcGoldStandard.RData")
##########################################################
