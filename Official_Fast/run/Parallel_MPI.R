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
#######
# For study: ens = 1007; niter = 10 ; Proposal Matrix scaling 0.001
#######
# Arguments
args = as.numeric(commandArgs(trailingOnly=TRUE))
print(args)
####################################################################################################
# Intialize
cycle=args[1] # Cycle <- passed in through PBS file
ens <-args[2] # Number of particles
niter<-args[3] # Number of MCMC iterations
# cycle=1 ; ens=1007 ; niter=6
####################################################################################################
setwd("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/")
source("run/mcmc_source_Tr.R")
inputDir<-"/glade/scratch/sanjib/runA/input"
outputDir<-"/glade/scratch/sanjib/runA/output"
####################################################################################
# Parallelize
library(snow);library(Rmpi);library(doParallel);library(foreach)
mp_type = "MPI"
cl <- parallel::makeCluster(spec = ens, type=mp_type)
doParallel::registerDoParallel(cl)
####################################################################################

####################################################################################################
####################################################################################################
# Importance Sampling - Parallelize
if(cycle==1){
  source("run/Initialize.R") 
}else{
  load(paste("output/temperVal_",cycle-1,".RData",sep=""))
  if(temperVal$cumulative>0.98){stop("Stopping Criterion Met")}
  MCMCtemperVal<-temperVal$cumulative
  load(paste("output/mhParameters_",cycle-1,".RData",sep=""))
  # Reweight based on cycle
  for (jobNum in 1:ens){ # Faster to do this than parallelize
    source("run/mcmc_source_Tr.R")
    jobPar<-parMat[jobNum,]
    llhd_t<-calcPF(cycle=cycle,jobNum=jobNum,llhdTemper=1, # we are tempering this by 1 to get the full likelihood
                   mcmcTemper=MCMCtemperVal, # This is the tempering value from the previous cycle's mutation stage
                   initResults=list(initResultsList[[1]][jobNum],initResultsList[[2]][[jobNum]]), 
                   priorPar=priorPar)
    # Save the file
    save(jobPar,llhd_t,file=paste(outputDir,"/PF_",cycle,"_",jobNum,".RData",sep=""))
  }
}

print("Stopped")
rm(list=setdiff(ls(), c("ens","cycle","niter","inputDir","outputDir")))
# ####################################################################################################
# ####################################################################################################
# Combine - Central Node
print("Central")
source("run/mcmc_source_Tr.R")
####################################################################################################
# Combine and Optimize
load(paste("output/temperVal_",cycle-1,".RData",sep=""))
combineIS(cycle=cycle,cumulTemp=temperVal$cumulative,ens=ens, prop=0.5)
print("Central Complete")
rm(list=setdiff(ls(), c("ens","cycle","niter","inputDir","outputDir")))

# ################################################################################################################################################################################################
# Metropolis Hastings  - Parallelize

foreach::foreach(jobNum=1:ens) %dopar% {
  setwd("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/")
  source("run/mcmc_source_Tr.R")
                   ################################################################################
                   load(paste("output/rsParameters_",cycle,".RData",sep=""))
                   # MCMC
  niter.mcmc = niter
  par.init<-parMat[jobNum,]
  load(paste("output/temperVal_",cycle,".RData",sep=""))
  MCMCtemperVal<-temperVal$cumulative
  temperVal<-temperVal$incremental
                   ##############################
                   ##############################
                   # Generate prorposal matrix for first sample
                   ##############################
                   ##############################
  if(cycle==1){
    CovMat<-genPropMat(cycle=cycle,scale=0.5)
  }else{
    CovMat<-genPropMat(cycle=cycle,scale=0.01)
  }
  # Note that we use a different function. This finds a good proposal based on the sample cov of particles form current cycle.
  initResults<-list(initResultsList[[1]][jobNum],initResultsList[[2]][[jobNum]])
  amcmc.out<-mcmcManual_tempered(iter=niter.mcmc,
                                 init=par.init,
                                 propCov=CovMat,
                                 inputDir = inputDir,
                                 outputDir = outputDir,
                                 obs = obs,
                                 priorPar = priorPar,
                                 jobNum=jobNum,
                                 temper=MCMCtemperVal,
                                 initResults=initResults,
                                 parNames=parNames)
  
                   ##################
                   save(amcmc.out,MCMCtemperVal,temperVal,
                        file=paste("/glade/scratch/sanjib/runA/output/MCMC_",cycle,"_1_",jobNum,".RData",sep=""))
                   ################################################################################
                   rm(list=setdiff(ls(), c("ens","cycle","niter","inputDir","outputDir")))
                 }


####################################################################################################
####################################################################################################
# Combine MH- Central Node
setwd("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/")
source("run/mcmc_source_Tr.R")
combineMH(cycle=cycle,ens=ens,stage=1) # Combine MH
combineTotalParticles(cycle=cycle)# combine Total Particles for Covariance Matrix generation(proposal)
