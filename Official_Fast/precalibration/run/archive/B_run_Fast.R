library(snow);library(Rmpi);library(doParallel);library(foreach)
setwd("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/precalibration")

# Parallelize
nprocs <-mpi.universe.size() - 1
print(nprocs)

mp_type = "MPI"
cl <- parallel::makeCluster(spec = nprocs, type=mp_type)
doParallel::registerDoParallel(cl)

# Values for runs
runIndex<-1:(nprocs*2)
outputMat<-foreach::foreach(jobNum=runIndex , .combine = "cbind" , .packages = c("mvtnorm","tmvtnorm","invgamma")) %dopar% {
  source("../run/rWrapper_Continuous_Fast.R")
  source("../run/mcmc_source_Tr.R")
  load(file="output/mhParameters_0.RData")
  inputDir<-"/glade/scratch/sanjib/fastTest/input"
  outputDir<-"/glade/scratch/sanjib/fastTest/output"
  jobPar<-parMat[jobNum,]
  modelEval(par = jobPar , j = jobNum , inputDir =inputDir , outputDir = outputDir)
}

save(outputMat,file = "output/fast_preCalibrationResults.RData")

