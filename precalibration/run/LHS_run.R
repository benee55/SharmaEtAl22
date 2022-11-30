rm(list=ls())
# Load data + Parameter Values
setwd("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/LHSoutput/")
load("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/LHSoutput/input/LHS_Inputs.RData")
source("run/rWrapper.R")
##################################################################
# Observations
dat<-read.table("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/SBYP1_obs.txt")
obsFull<-dat[(30*365+10):(30*365+10+365-1),]
obs<-as.numeric(obsFull[,5])*0.0283 # Need to convert 

##################################################################
# Priors
priorPar<-rbind(c(0.02 , 0.02),  # sigma2
                c(0 , 0.3), # PCTIM
                c(0 , 0.5), # ADIMP
                c(-1.1 , -0.2), # UZTWM
                c(-2.5 , -0.25), # LZTWM
                c(-2.8 , -0.2), # LZFSM
                c(-3.2 , -0.3), # LZFPM
                c(-3.8 , -0.1), # LZSK
                c(-4.5 , -0.1), # UZFWM
                c(-3.5 , -0.1), # REXP
                c(-3.5 , -0.1), # UZK
                c(0.5,4.5), # rutpix_Q0CHN
                c(0.3,2.25)) # rutpix_QMCHN

##################################################################
##################################################################
inputDir<-"/gpfs/group/kzk10/default/private/hydrocalib/SGrove/LHSoutput/input"
outputDir<-'/gpfs/group/kzk10/default/private/hydrocalib/SGrove/LHSoutput/output'

# Parallelize
library(snow);library(doParallel);library(foreach); library(Rmpi)
nprocs <- 259
mp_type = "MPI" # PSOCK or MPI
cl <- parallel::makeCluster(nprocs, type=mp_type)
doParallel::registerDoParallel(cl)

foreach::foreach(j=1:1000,
                 .packages=c('invgamma')) %dopar% {
                   par<-c(runif(1,0.001,10),parMat[j,])
                   lpostVal<-logPosterior(par = par , priorPar = priorPar , 
                                          obs = obs, inputDir = inputDir , 
                                          outputDir = outputDir , j = j)
                   
                   save(lpostVal,par,
                        file=paste("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/LHSoutput/output/output",j,".RData",sep=""))
                   
                 }


