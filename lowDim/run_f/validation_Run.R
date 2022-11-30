jobIndex <- as.numeric(commandArgs(trailingOnly=TRUE))
library(snow);library(doParallel);library(foreach)
setwd("/glade/u/home/sanjib/FamosHydroModel/lowDim/output_f")

# Parallelize
nprocs <-35
print(nprocs)

mp_type = "PSOCK"
cl <- parallel::makeCluster(spec = nprocs, type=mp_type)
doParallel::registerDoParallel(cl)

# Values for runs
if(jobIndex==58){
  runIndex<-((jobIndex-1)*(nprocs))+(1:(nprocs))
  runIndex<-runIndex[runIndex<=2015]
}else{
  runIndex<-((jobIndex-1)*(nprocs))+(1:(nprocs))  
}

print("Begin")
outputMat<-foreach::foreach(jobNum=runIndex , .combine = "cbind" , .packages = c("mvtnorm","tmvtnorm","invgamma")) %dopar% {
  source("../run_f/rWrapper_validation.R")
  source("../run_f/mcmc_source_Tr.R")
  load(file="mhParameters_4.RData")
  inputDir<-"/glade/scratch/sanjib/validationLowDim/input"
  outputDir<-"/glade/scratch/sanjib/validationLowDim/output"
  jobPar<-rep2orig(parMat[jobNum,])
  modelEval(par = jobPar , j = jobNum , inputDir =inputDir , outputDir = outputDir)
}
print('finish')
save(outputMat,file = paste("validationResults",jobIndex,".RData",sep=""))

