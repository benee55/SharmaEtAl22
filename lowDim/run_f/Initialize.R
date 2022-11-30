# # Copyright (C) 2018 Ben S. Lee
#email: skl5261@psu.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

# Generate intial parameter set from pre-calibration
setwd("/glade/u/home/sanjib/FamosHydroModel/lowDim/")

for(i in 1:14){
  if(i==1){
    load(paste("precalibration/preCalibrationResults",i,".RData",sep=""))
    bar<-unlist(outputMat[2,])
    modelOutput<-matrix(bar, nrow=length(outputMat[2,]), ncol=length(outputMat[2,1][[1]]), byrow = TRUE)
  }else{
    load(paste("precalibration/preCalibrationResults",i,".RData",sep=""))
    bar<-unlist(outputMat[2,])
    foo<-matrix(bar, nrow=length(outputMat[2,]),ncol=length(outputMat[2,1][[1]]), byrow = TRUE)
    modelOutput<-rbind(modelOutput,foo)
  }
}


# Parameters
load("precalibration/mhParameters_0.RData")
parMat<-parMat[1:nrow(modelOutput),]
parMat[,1]<-1825.215^2 # Results from First Calibration
parMat<-t(apply(parMat, 1, orig2rep)) # Reparameterize
save(parMat,file="output_f/mhParameters_0.RData")

# Functions to compute posterior
logLikelihood_temper<-function(par, obs , temper , output, obsInd){
  sigma2<-par[1] # Variance Parameter
  extremeOutput<-output[obsInd]
  llhd<-temper*sum(dnorm(x=obs, mean = extremeOutput , sd= sqrt(sigma2), log = TRUE)) # COmpute Likelihood
  output<-list(extremeOutput,output)
  return(list(llhd,output))
}

for(jobNum in 1:nrow(parMat)){
  jobPar<-parMat[jobNum,]
  llhd_t<-logLikelihood_temper(par =jobPar, obs = obs ,  temper = 1 , output = modelOutput[jobNum,] , obsInd=obsInd)
  save(jobPar,llhd_t,file=paste("/glade/scratch/sanjib/run_f/output/PF_",cycle,"_",jobNum,".RData",sep=""))
}

# Tempering values
temperVal<-list()
temperVal$cumulative<-0
temperVal$incremental<-0
save(temperVal,file="output_f/temperVal_0.RData")

# Generate Covariance matrix for proposal
llhd<-apply(modelOutput, 2, function(x){sum(dnorm(x = obs, mean = x[obsInd] , sd = 1825.215, log = TRUE))})
# 
covIndex<-which(llhd>quantile(llhd,probs = 0.9))
keepParMat<-parMat[covIndex,-1]
# Use acceptance ratio from Rosenthal et al. 
CovMat<-cov(keepParMat)*((2.38^2)/ncol(keepParMat))
CovMat<-CovMat+diag(ncol(CovMat))*(0.001*diag(CovMat))
save(CovMat, file="output_f/BeginCovMat_Tr.RData")

