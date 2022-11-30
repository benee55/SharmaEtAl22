library(snow);library(doParallel);library(foreach)
library(Rmpi);
setwd("/glade/u/home/sanjib/FamosHydroModel/lowDim")

# Parallelize
nprocs <-mpi.universe.size() - 1
# nprocs <-4
print(nprocs)

mp_type = "MPI"
# mp_type = "PSOCK"
cl <- parallel::makeCluster(spec = nprocs, type=mp_type)
doParallel::registerDoParallel(cl)
load(file="input/design.RData")
ens<-nrow(parMat); rm(parMat)
# ens<-4
# Values for runs
outputMat<-foreach::foreach(jobNum=1:ens , .combine = "cbind") %dopar% {
  source("run/rWrapper_Continuous_projection.R")
  source("run/mcmc_source_Tr.R")
  load(file="input/design.RData")
  inputDir<-"/glade/scratch/sanjib/lowDim/input"
  outputDir<-"/glade/scratch/sanjib/lowDim/output"
  jobPar<-parMat[jobNum,]
  modelEval(par = jobPar , j = jobNum+2000 , inputDir =inputDir , outputDir = outputDir)
}

save(outputMat,file = "output/modelRuns_projection.RData")

