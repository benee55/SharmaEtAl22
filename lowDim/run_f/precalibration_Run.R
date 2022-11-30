jobIndex <- as.numeric(commandArgs(trailingOnly=TRUE))

library(snow);library(Rmpi);library(doParallel);library(foreach)
setwd("/glade/u/home/sanjib/FamosHydroModel/lowDim/precalibration")

# Parallelize
nprocs <-mpi.universe.size() - 1
print(nprocs)

mp_type = "MPI"
cl <- parallel::makeCluster(spec = nprocs, type=mp_type)
doParallel::registerDoParallel(cl)

# Values for runs
runIndex<-((jobIndex-1)*(nprocs))+(1:(nprocs))
outputMat<-foreach::foreach(jobNum=runIndex , .combine = "cbind" , .packages = c("mvtnorm","tmvtnorm","invgamma")) %dopar% {
  source("../run_f/rWrapper_Continuous.R")
  source("../run_f/mcmc_source_Tr.R")
  load(file="../precalibration/mhParameters_0.RData")
  inputDir<-"/glade/scratch/sanjib/precalibrationLowDim/input"
  outputDir<-"/glade/scratch/sanjib/precalibrationLowDim/output"
  jobPar<-parMat[jobNum,]
  modelEval(par = jobPar , j = jobNum , inputDir =inputDir , outputDir = outputDir)
}

save(outputMat,file = paste("preCalibrationResults",jobIndex,".RData",sep=""))

