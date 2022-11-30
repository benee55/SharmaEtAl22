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
                c(0.3,3.4)) # rutpix_QMCHN

##################################################################
##################################################################
inputDir<-"/gpfs/group/kzk10/default/private/hydrocalib/SGrove/LHSoutput/input"
outputDir<-'/gpfs/group/kzk10/default/private/hydrocalib/SGrove/LHSoutput/output'

#Individual Components
j=5
# writeInput( par = parMat[1,] , j = j , dir = inputDir) # Function 1: Write Input
# writeOutput( j = j , dir = outputDir) # Function 2: Create Directory
# runHydroModel( j = j , dir = inputDir) # Function 3: Run Model + Save File
# dat<-readOutput( j = j, dir = outputDir)  # Function 4: Read Output File
# dat2<-runReadModel( j = j , inputDir =inputDir , outputDir = outputDir) # Function 5: Combined Run and Read
# Run posterior
par<-c(1,parMat[1,])
lpostVal<-logPosterior(par = par , priorPar = priorPar , 
                        obs = obs, inputDir = inputDir , 
                        outputDir = outputDir , j = j)
print(lpostVal)