jobIndex <- as.numeric(commandArgs(trailingOnly=TRUE))

library(snow);library(Rmpi);library(doParallel);library(foreach)
setwd("/glade/u/home/sanjib/FamosHydroModel/Official_Fast")

# Parallelize
nprocs <-mpi.universe.size() - 1
print(nprocs)

mp_type = "MPI"
cl <- parallel::makeCluster(spec = nprocs, type=mp_type)
doParallel::registerDoParallel(cl)

# Values for runs
runIndex<-((jobIndex-1)*(nprocs))+(1:(nprocs))
outputMat<-foreach::foreach(jobNum=runIndex , .combine = "cbind" , .packages = c("mvtnorm","tmvtnorm","invgamma")) %dopar% {
  source("run_validation/rWrapper.R")
  load(file="precalibration/output/mhParameters_0.RData")
  inputDir<-"/glade/scratch/sanjib/validation/input/precalibration"
  outputDir<-"/glade/scratch/sanjib/validation/output/precalibration"
  jobPar<-parMat[jobNum,]
  modelEval(par = jobPar , j = jobNum , inputDir =inputDir , outputDir = outputDir)
}

save(outputMat,file = paste("output_validation/preCalibrationResults",jobIndex,".RData",sep=""))

