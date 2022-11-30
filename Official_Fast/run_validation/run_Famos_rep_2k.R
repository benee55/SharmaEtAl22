
library(snow);library(Rmpi);library(doParallel);library(foreach)
setwd("/glade/u/home/sanjib/FamosHydroModel/Official_Fast")

# Parallelize
nprocs <-mpi.universe.size() - 1
print(nprocs)

mp_type = "MPI"
cl <- parallel::makeCluster(spec = nprocs, type=mp_type)
doParallel::registerDoParallel(cl)

# Values for runs
runIndex<-1:2015
outputMat<-foreach::foreach(jobNum=runIndex , .combine = "cbind" , .packages = c("mvtnorm","tmvtnorm","invgamma")) %dopar% {
  load(file="output_rep/mhParameters_4.RData")
  source("run_rep/mcmc_source_Tr.R")
  inputDir<-"/glade/scratch/sanjib/validation/input/Famos_rep"
  outputDir<-"/glade/scratch/sanjib/validation/output/Famos_rep"
  jobPar<-rep2orig(parMat[jobNum,])
  source("run_validation/rWrapper.R")
  modelEval(par = jobPar , j = jobNum , inputDir =inputDir , outputDir = outputDir)
}

save(outputMat,file ="output_validation/famosResults_rep_2k.RData")

