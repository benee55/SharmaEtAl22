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
args = as.numeric(commandArgs(trailingOnly=TRUE))
print(args)
####################################################################################################
# Intialize
cycle=args[1] # Cycle <- passed in through PBS file
niter<-args[2] # Number of MCMC iterations
ens <-mpi.universe.size() - 1 # Number of particles
cycle=1 ; niter=6 ; ens<-71
####################################################################################################
setwd("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/")
source("run/mcmc_source_Tr.R")
inputDir<-"/glade/scratch/sanjib/runA/input"
outputDir<-"/glade/scratch/sanjib/runA/output"
####################################################################################
# Parallelize
# library(snow);library(Rmpi);library(doParallel);library(foreach)
# mp_type = "MPI"
# cl <- parallel::makeCluster(spec = ens, type=mp_type)
# doParallel::registerDoParallel(cl)
####################################################################################
# foreach::foreach(jobNum=1:ens) %dopar% {
setwd("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/")
source("run/mcmc_source_Tr.R")
inputDir<-"/glade/scratch/sanjib/runA/input"
outputDir<-"/glade/scratch/sanjib/runA/output"
################################################################################
load(paste("output/rsParameters_",cycle,".RData",sep=""))
# MCMC
jobNum<-1
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
CovMat<-genPropMat(cycle=cycle,scale=1)   # Note that we use a different function. This finds a good proposal based on the sample cov of particles form current cycle.

initResults<-list(initResultsList[[1]][jobNum],initResultsList[[2]][[jobNum]])

# Test Finished
# curResults<-logPosterior_temper(par=par.init,
#                                 priorPar = priorPar,
#                                 obs = obs,
#                                 inputDir = inputDir,
#                                 outputDir = outputDir,
#                                 j=jobNum , 
#                                 temper=MCMCtemperVal)
# Test FInished
# output<-modelEval(par=par.init, j=jobNum , inputDir=inputDir , outputDir=outputDir) 


  # set.seed(jobNum*1234*cycle) #set seed
  ##################
  ##################
  ##################
  ##################
  ##################
                     # TO DO
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
  ##################
  ##################
  ##################
  save(amcmc.out,MCMCtemperVal,temperVal,
       file=paste("/glade/scratch/sanjib/runA/output/MCMC_",cycle,"_1_",jobNum,".RData",sep=""))
  ################################################################################
  rm(list=setdiff(ls(), c("ens","cycle","niter","inputDir","outputDir")))
}