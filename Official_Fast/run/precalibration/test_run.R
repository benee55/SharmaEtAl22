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
# Arguments

####################################################################################################

rm(list=ls())
setwd("~/Dropbox/hydroFamos/run/precalibration")
####################################################################################
# Parallelize
library(snow);library(doParallel);library(foreach)
# library(Rmpi);
nprocs <- 7
mp_type = "PSOCK" # PSOCK or MPI
cl <- parallel::makeCluster(nprocs, type=mp_type)
doParallel::registerDoParallel(cl)

source("test_PreCalibration.R")
for(hj in 1:11){
  # print(j)

# Intialize
cycle=hj # Cycle <- passed in through PBS file
ensembleN=3000 # Total number of particles
niter<-10 # Number of MCMC iterations
####################################################################################################
setwd("~/Dropbox/hydroFamos/run/precalibration")
source("test_Source.R")
inputDir<-"~/Dropbox/hydroFamos/run/precalibration/input"
outputDir<-"~/Dropbox/hydroFamos/run/precalibration/output"
####################################################################################

####################################################################################################
####################################################################################################
# Importance Sampling - Parallelize
if(cycle==1){
  source("test_Initialize.R")  
}else{
  load(paste("output/temperVal_",cycle-1,".RData",sep=""))
  print(temperVal$cumulative)
  if(temperVal$cumulative>0.98){stop("Stopping Criterion Met")}
  load(paste("output/mhParameters_",cycle-1,".RData",sep=""))
}

#Compute Initial Weigths for Cycle
if(cycle==1){
  
  foreach::foreach(jobNum=1:ensembleN) %dopar% {
    source("test_Source.R")
    load(paste("output/mhParameters_",(cycle-1),".RData",sep=""))
    jobPar<-parMat[jobNum,]
    
    llhd_t<-logLikelihood_temper(par =jobPar, obs = obs ,  x=inputX, temper = 1)
    save(jobPar,llhd_t,file=paste("output/PF_",cycle,"_",jobNum,".RData",sep=""))
  }
  
}else{
  load(paste("output/temperVal_",cycle-1,".RData",sep=""))
  MCMCtemperVal<-temperVal$cumulative
  for(jobNum in 1:ensembleN){
    load(paste("output/mhParameters_",(cycle-1),".RData",sep=""))
    jobPar<-parMat[jobNum,]
    llhd_t<-calcPF(cycle=cycle,jobNum=jobNum,llhdTemper=1, # we are tempering this by 1 to get the full likelihood
                   mcmcTemper=MCMCtemperVal, # Cumulative tempered value on cycle-1 MCMC step
                   priorPar=priorPar,
                   initResults=list(initResultsList[[1]][jobNum],initResultsList[[2]][[jobNum]]))
    # Save the file
    save(jobPar,llhd_t,file=paste("output/PF_",cycle,"_",jobNum,".RData",sep=""))
  }
}

print("Stopped")
rm(list=setdiff(ls(), c("ensembleN","cycle","niter","inputDir","outputDir","hj")))
# ####################################################################################################
# ####################################################################################################
# Combine - Central Node
print("Central")
source("test_Source.R")
####################################################################################################
# Combine and Optimize
load(paste("output/temperVal_",cycle-1,".RData",sep=""))
combineIS(cycle=cycle,cumulTemp=temperVal$cumulative,prop=0.5)

print("Central Complete")
rm(list=setdiff(ls(), c("ensembleN","cycle","niter","inputDir","outputDir","hj")))

# ################################################################################################################################################################################################
# # Metropolis Hastings  - Parallelize
# #Parallel Settings

foreach::foreach(jobNum=1:ensembleN,
                 .packages=c('mvtnorm','TruncatedNormal')) %dopar% {
                   setwd("~/Dropbox/hydroFamos/run/precalibration/")
                   source("test_Source.R")
                   ################################################################################
                   load(paste("output/rsParameters_",cycle,".RData",sep=""))
                   # MCMC
                   niter.mcmc = niter
                   par.init<-parMat[jobNum,]
                   ##############################
                   load(paste("output/temperVal_",cycle,".RData",sep=""))
                   MCMCtemperVal<-temperVal$cumulative
                   temperVal<-temperVal$incremental
                   ##############################
                   ##############################
                   # Generate prorposal matrix for first sample
                   ##############################
                   ##############################
                   CovMat<-genPropMat(cycle=cycle,scale=0.001)   # Note that we use a different function. This finds a good proposal based on the sample cov of particles form current cycle.
                   initResults<-list(initResultsList[[1]][jobNum],initResultsList[[2]][[jobNum]])
                   # set.seed(jobNum*1234*cycle) #set seed
                   ##################
                   # TO DO
                   amcmc.out<-mcmcManual_tempered(iter=niter.mcmc,
                                                  init=par.init,
                                                  propCov=CovMat,
                                                  obs = obs,
                                                  priorPar = priorPar,
                                                  temper=MCMCtemperVal,
                                                  x=inputX,
                                                  initResults=initResults,
                                                  parNames=parNames)

                   ##################
                   ##################
                   save(amcmc.out,MCMCtemperVal,temperVal,
                        file=paste("output/MCMC_",cycle,"_1_",jobNum,".RData",sep=""))
                   ################################################################################
                   
                   rm(list=setdiff(ls(), c("ensembleN","cycle","niter","inputDir","outputDir","hj")))
                 }
# 
# 
# ####################################################################################################
# ####################################################################################################
# # Combine MH- Central Node
print("Mutation Complete")
setwd("~/Dropbox/hydroFamos/run/precalibration/output/")
source("../test_Source.R")
combineMH(cycle=cycle,ens=ensembleN,stage=1) # Combine MH
combineTotalParticles(cycle=cycle)# combine Total Particles for Covariance Matrix generation(proposal)
}



rm(list = ls())
load("output/mhParameters_4.RData")

load("MCMCOutput.RData")
mean(acceptVect)
par(mfrow=c(2,2))
plot(density(parMat[,1]))
lines(density(mcmcOutput[,3]))
plot(density(parMat[,2]))
lines(density(mcmcOutput[,1]))
plot(density(parMat[,3]))
lines(density(mcmcOutput[,2]))

rm(list = ls())
load("output/mhParameters_1.RData")
plot(density(parMat[,2]))
for(k in 2:4){
  load(paste("output/mhParameters_",k,".RData",sep=""))
  lines(density(parMat[,2]))  
}

plot(x=1:length(initResultsList[[2]][[1]]), y=initResultsList[[2]][[2]], ylim=range(unlist(initResultsList[[2]])))
for(k in 1:length((initResultsList[[2]]))){
  lines(x=1:length(initResultsList[[2]][[k]]), y=initResultsList[[2]][[k]], col="gray")
}
points(x=1:length(initResultsList[[2]][[1]]), y=initResultsList[[2]][[2]], col="blue", pch=16)
    
    
for(k in 1:5){
  load(paste("output/mhParameters_",k,".RData",sep=""))
  # print(mean(acceptVect))
  print(mean(acceptVect!=0))
}


rm(list = ls())
load("output/rsParameters_4.RData")
plot(density(parMat[,3]), col="red")
load("output/mhParameters_4.RData")
lines(density(parMat[,3]),col="blue")


rm(list=ls())
load("output/masterTotalParticles.RData")
dim(masterTotalParticles)
nrow(masterTotalParticles)/5/10
sum(weights)
hist(weightVect)
hist(weights)
1/sum((weights^2))
length(unique(reSampleInd))

