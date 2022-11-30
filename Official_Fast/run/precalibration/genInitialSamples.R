rm(list=ls())
####################################################################################################
setwd("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/famos/Official_Fast/")
library(snow);library(Rmpi);library(doParallel);library(foreach)
ens<-nprocs <-39
mp_type = "MPI" # PSOCK or MPI
cl <- parallel::makeCluster(spec = nprocs, type="MPI")
doParallel::registerDoParallel(cl)

source("run/mcmc_source_Tr.R")
source("run/Initialize.R")

outputMat<-foreach::foreach(jobNum=1:ens) %dopar% {
  source("run/mcmc_source_Tr.R")
  source("run/rWrapper.R")
  load("output/mhParameters_0.RData")
  inputDir<-"/gpfs/group/kzk10/default/private/hydrocalib/SGrove/famos/Official_Fast/input"
  outputDir<-'/gpfs/group/kzk10/default/private/hydrocalib/SGrove/famos/Official_Fast/output'
  jobPar<-parMat[jobNum,]
  modelEval( par = jobPar , j = jobNum , inputDir =inputDir , outputDir = outputDir)
  
}

# rm(list=ls())
# setwd("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/famos/Official_Fast/")
# source("run/mcmc_source_Tr.R")
# source("run/rWrapper.R")
# load("output/mhParameters_0.RData")
# inputDir<-"/gpfs/group/kzk10/default/private/hydrocalib/SGrove/famos/Official_Fast/input"
# outputDir<-'/gpfs/group/kzk10/default/private/hydrocalib/SGrove/famos/Official_Fast/output'
# load("input/obsData.RData")
# ens=390
# outputMat<-matrix(nrow=ens,ncol=19)
# for(j in 1:ens){
#   print(j)
#   for(i in 1:5){
#     if(i==1){ # Read model output - First
#       output<-readOutput( j = j , interval=i, dir = outputDir)
#     }else{  # Read model output - Next
#       output<-c(output,readOutput( j = j , interval=i, dir = outputDir))
#     }
#   }
#   outputMat[j,]<-output
# }

save(outputMat, file="output/B_preCalibrationResults.RData")
