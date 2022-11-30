library(snow);library(Rmpi);library(doParallel);library(foreach)
setwd("/glade/u/home/sanjib/FamosHydroModel/Official_Fast")

# Parallelize
nprocs <-mpi.universe.size() - 1
print(nprocs)

mp_type = "MPI"
cl <- parallel::makeCluster(spec = nprocs, type=mp_type)
doParallel::registerDoParallel(cl)

# Values for runs
runIndex<-1:1007
outputMat<-foreach::foreach(jobNum=runIndex , .combine = "cbind" , .packages = c("mvtnorm","tmvtnorm","invgamma")) %dopar% {
  source("run_validation/rWrapper.R")
  load(file="output/mhParameters_6.RData")
  inputDir<-"/glade/scratch/sanjib/validation/input/Famos_variance"
  outputDir<-"/glade/scratch/sanjib/validation/output/Famos_variance"
  jobPar<-parMat[jobNum,]
  modelEval(par = jobPar , j = jobNum , inputDir =inputDir , outputDir = outputDir)
}

save(outputMat,file ="output_validation/famosResults_variance.RData")

