rm(list=ls())
setwd("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/precalibration")
source("../run/rWrapper_Continuous.R")
source("../run/mcmc_source_Tr.R")

source("../run/rWrapper_Continuous.R")
source("../run/mcmc_source_Tr.R")
load("output/mhParameters_0.RData")
inputDir<-"/glade/scratch/sanjib/precalibration/input"
outputDir<-"/glade/scratch/sanjib/precalibration/output"
jobNum=1
jobPar<-parMat[jobNum,]
outputMat<- modelEval( par = jobPar , j = jobNum , inputDir =inputDir , outputDir = outputDir)
save(outputMat,file = "output/preCalibrationResults_single.RData")
