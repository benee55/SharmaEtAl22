
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
  load(file="output_rep_var/mhParameters_4.RData")
  source("run_rep_var/mcmc_source_Tr.R")
  inputDir<-"/glade/scratch/sanjib/validation/input"
  outputDir<-"/glade/scratch/sanjib/validation/output"
  jobPar<-rep2orig(parMat[jobNum,])
  source("run_validation/rWrapper.R")
  modelEval(par = jobPar , j = jobNum , inputDir =inputDir , outputDir = outputDir)
}

save(outputMat,file ="output_rep_var/famosResults_rep_var.RData")

