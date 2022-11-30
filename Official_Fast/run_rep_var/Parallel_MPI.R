# # Copyright (C) 2022 Ben S. Lee
#email: slee287@gmu.edu
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
# Parallel_MPI.R: This file runs the FaMoS calibration. 
#######
# Initial Arguments passed in through PBS File
args = as.numeric(commandArgs(trailingOnly=TRUE))
print(args)
####################################################################################################
# Initialize Settings
cycle=args[1] # Cycle <- passed in through PBS file
ens <-args[2] # Number of particles
niter<-args[3] # Number of MCMC iterations for the mutation cycle
####################################################################################################
# Set the Input/Output Directories and load the helper function file mcmc_source_Tr.R
setwd("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/")
source("run_rep_var/mcmc_source_Tr.R")
inputDir<-"/glade/scratch/sanjib/runA_rep_var/input"
outputDir<-"/glade/scratch/sanjib/runA_rep_var/output"
####################################################################################
# Intiialize Parallelization via RMPI
library(snow);library(Rmpi);library(doParallel);library(foreach)
mp_type = "MPI"
cl <- parallel::makeCluster(spec = ens, type=mp_type)
doParallel::registerDoParallel(cl)
####################################################################################################
####################################################################################################
# Part 1A: Importance Sample 
if(cycle==1){ # Only for First Cycle
  source("run_rep_var/Initialize.R") # Initialize the particles according to Initialize.R
}else{# For subsequent cycles
  load(paste("output_rep_var/temperVal_",cycle-1,".RData",sep="")) # Load tempering values from previous cycle
  if(temperVal$cumulative>0.98){stop("Stopping Criterion Met")} # Stop if 98%+ of likelihood has been incorporated
  MCMCtemperVal<-temperVal$cumulative # Load in the cumulative tempering value 
  load(paste("output_rep_var/mhParameters_",cycle-1,".RData",sep="")) # Load in the mutated particles from the previous cycle
  # In the Following Step, we reweight the particles (from previous cycle) using the "Full Likelihood Value" 
  # We incorporate the full likelihood such that we can optimize the incorporate schedule with respect to the Effective Sample Size (ESS). 
  for (jobNum in 1:ens){
    source("run_rep_var/mcmc_source_Tr.R") # Load Helper File
    jobPar<-parMat[jobNum,] # Parameter Set
    llhd_t<-calcPF(cycle=cycle,jobNum=jobNum,llhdTemper=1, # we are tempering this by 1 to get the full likelihood
                   mcmcTemper=MCMCtemperVal, # This is the tempering value from the previous cycle's mutation stage
                   initResults=list(initResultsList[[1]][jobNum],initResultsList[[2]][[jobNum]]), 
                   priorPar=priorPar)
    # Save the file
    save(jobPar,llhd_t,file=paste(outputDir,"/PF_",cycle,"_",jobNum,".RData",sep="")) # Save file with the parameters and likelihood
  }
}

print("Stopped")
rm(list=setdiff(ls(), c("ens","cycle","niter","inputDir","outputDir"))) # Remove everything except the pertinent information
# ####################################################################################################
# ####################################################################################################
# Part 1B: Compute the optimal likelihood incorporation power (gamma) with respect to effective Sample Size. 
# In this study, we will set the gamma_min=0.1 and ESS_threshold=0.5
print("Central")
source("run_rep_var/mcmc_source_Tr.R") # Load Helper File
load(paste("output_rep_var/temperVal_",cycle-1,".RData",sep="")) # Load tempering value data
combineIS(cycle=cycle,cumulTemp=temperVal$cumulative,ens=ens, prop=0.5) # Compute the optimal incorporation power (gamma_t) and combine particles
print("Central Complete")
rm(list=setdiff(ls(), c("ens","cycle","niter","inputDir","outputDir")))

#################################################################################################################################################################################################
# Part 2A: Mutation Step
#################################################################################################################################################################################################
# Using the reweighted particles from Part 1B, we mutate each particles using the Metropolis-Hastings Kernel
# We employ RMPI and distributed computing for this part. 
foreach::foreach(jobNum=1:ens) %dopar% {
  setwd("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/") # Set directory
  source("run_rep_var/mcmc_source_Tr.R") # Load Helper File
  load(paste("output_rep_var/rsParameters_",cycle,".RData",sep="")) # Load data from previous step
  # Begin Mutation
  niter.mcmc = niter # Set the number of iterations/mutations per Markov Chain 
  par.init<-parMat[jobNum,] # Initial Parameters
  load(paste("output_rep_var/temperVal_",cycle,".RData",sep="")) # Load tempering file for the cumulative tempering value
  MCMCtemperVal<-temperVal$cumulative # Cumulative tempering value sum_{i=1}^{t}gamma_t
  temperVal<-temperVal$incremental # Incremental tempering value gamma_t

  # Generate a proposal covariance matrix for the MH updates
  if(cycle==1){
    CovMat<-genPropMat(cycle=cycle,scale=1) # Cycle 1
  }else{
    CovMat<-genPropMat(cycle=cycle,scale=0.2) # All subsequent cycles. We choose these values based on the overall performance
  }
  initResults<-list(initResultsList[[1]][jobNum],initResultsList[[2]][[jobNum]]) # Combine the output results from Part 1B. 
  # Particle Mutation: This performs an all-at-once update of the Metropolis-Hastings algorithm. 
  amcmc.out<-mcmcManual_tempered(iter=niter.mcmc, # Number of iterations
                                 init=par.init, # Initial Values
                                 propCov=CovMat, # Proposal Covariance Matrix
                                 inputDir = inputDir, # Input Directory
                                 outputDir = outputDir, # Output Directory
                                 obs = obs, # Observed Data
                                 priorPar = priorPar, # Prior Distribution Parameters (see mcmc_source_Tr.R)
                                 jobNum=jobNum, # Particle Index (1:ens)
                                 temper=MCMCtemperVal, # Cumulative tempering value sum_{i=1}^{t}gamma_t
                                 initResults=initResults, # Combine the output results from Part 1B. 
                                 parNames=parNames) # Parameter Names
  
  ##################
  save(amcmc.out,MCMCtemperVal,temperVal, # Save Data
       file=paste("/glade/scratch/sanjib/runA_rep_var/output/MCMC_",cycle,"_1_",jobNum,".RData",sep=""))
  ################################################################################
  rm(list=setdiff(ls(), c("ens","cycle","niter","inputDir","outputDir"))) # Remove unnecessary objects
}


####################################################################################################
####################################################################################################
#################################################################################################################################################################################################
# Part 2B: Combine Particles from Part 2A
#################################################################################################################################################################################################
setwd("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/") # Set working Directory
source("run_rep_var/mcmc_source_Tr.R") # Load Helper File
combineMH(cycle=cycle,ens=ens,stage=1) # Combine results from Part 2A. 
combineTotalParticles(cycle=cycle)# Combine Total Particles to generate the Proposal Covariance Matrix
