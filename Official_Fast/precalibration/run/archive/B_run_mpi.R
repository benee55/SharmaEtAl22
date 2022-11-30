library(snow);library(Rmpi);library(doParallel);library(foreach)
setwd("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/precalibration")

# Parallelize
nprocs <-mpi.universe.size() - 1
print(nprocs)

mp_type = "MPI"
cl <- parallel::makeCluster(spec = nprocs, type=mp_type)
doParallel::registerDoParallel(cl)

# Values for runs
runIndex<-3080:(3080+nprocs-1)

outputMat<-foreach::foreach(jobNum=runIndex , .combine = "cbind" , .packages = c("mvtnorm","tmvtnorm","invgamma")) %dopar% {
  source("../run/rWrapper_Continuous.R")
  source("../run/mcmc_source_Tr.R")
  load(file="output/mhParameters_0.RData")
  inputDir<-"/glade/scratch/sanjib/precalibration/input"
  outputDir<-"/glade/scratch/sanjib/precalibration/output"
  jobPar<-parMat[jobNum,]
  modelEval(par = jobPar , j = jobNum , inputDir =inputDir , outputDir = outputDir)
}

save(outputMat,file ="output/preCalibrationResults_Additional.RData")

