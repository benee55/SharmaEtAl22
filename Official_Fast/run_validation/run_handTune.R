setwd("/glade/u/home/sanjib/FamosHydroModel/Official_Fast")

source("run/mcmc_source_Tr.R")
inputDir<-"/glade/scratch/sanjib/validation/input/handTune"
outputDir<-"/glade/scratch/sanjib/validation/output/handTune"
jobPar<-c(2, # Variance - Not important in this evaluation 
          0 , # PCTIM 
              0.1 , # ADIMP 
              -0.7 , # UZTWM
              -0.7 , # LZTWM
              -1.3 , # LZFSM
              -1.2 , # LZFPM
              -2.3 , # LZSK
              1.0 , # snow_SCF
              -1.5 , # REXP
              -1.3 , # UZK
              3.0 , # rutpix_Q0CHN
              1.0) # rutpix_QMCHN

source("run_validation/rWrapper.R")
outputMat<-modelEval(par = jobPar , j = 1 , inputDir =inputDir , outputDir = outputDir)
save(outputMat , file = "output_validation/handTuneResults.RData")
source("run/rWrapper_Continuous.R")
modelRun<-modelEval(par = jobPar , j = 2 , inputDir =inputDir , outputDir = outputDir)
save(outputMat , modelRun,file = "output_validation/handTuneResults.RData")