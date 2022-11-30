rm(list=ls())
library(snow);library(Rmpi);library(doParallel);library(foreach)
setwd("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/precalibration")
source("../run/rWrapper_Continuous.R")
source("../run/mcmc_source_Tr.R")

# Precalibration
ensembleN<-10000
parMat<-apply(boundMat,1,function(x,ens){runif(ens,min=x[1],max=x[2])},ens=ensembleN)
parMat<-cbind(rinvgamma(ensembleN,shape = priorPar[1,1], rate = priorPar[1,2]),parMat)
save(parMat, file="output/mhParameters_0.RData")

nprocs <-5*20-1
mp_type = "MPI" # PSOCK or MPI
cl <- parallel::makeCluster(spec = nprocs, type="MPI")
doParallel::registerDoParallel(cl)

outputMat<-foreach::foreach(jobNum=1:ensembleN , .combine = "cbind") %dopar% {
  source("../run/rWrapper_Continuous.R")
  source("../run/mcmc_source_Tr.R")
  load("output/mhParameters_0.RData")
  inputDir<-"/gpfs/scratch/skl5261/precalibration/input"
  outputDir<-"/gpfs/scratch/skl5261/precalibration/output"
  jobPar<-parMat[jobNum,]
  modelEval( par = jobPar , j = jobNum , inputDir =inputDir , outputDir = outputDir)
}
save(outputMat,file = "preCalibrationResults.RData")

