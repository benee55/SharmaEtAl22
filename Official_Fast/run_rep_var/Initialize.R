# # Copyright (C) 2022 Ben S. Lee
#email: slee287@gmu.edu
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

####################################################################################
####################################################################################
# Initialize Data for FaMoS Loading
####################################################################################
####################################################################################
# Set working Directory
setwd("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/")

# Add parameter settings from the precalibration runs (random sampling of parameters)
for(i in 1:14){
  load(paste("precalibration/output/preCalibrationResults",i,".RData",sep=""))
  if(i==1){
    modelOutput<-outputMat[-nrow(outputMat),]
  }else{
    modelOutput<-cbind(modelOutput,outputMat[-nrow(outputMat),])
  }
}

# Parameters from the initial sampling of parameters
load("precalibration/output/mhParameters_0.RData")
parMat<-parMat[1:ncol(modelOutput),]
parMat[,1]<-1825.215^2 # Results from First Calibration
parMat<-t(apply(parMat, 1, orig2rep)) # Transform the parameters to the reparameterized settings. 
# Reparameterization aids in the SMC
save(parMat,file="output_rep_var/mhParameters_0.RData") #Save Data

# Functions to compute likelihood function
logLikelihood_temper<-function(par, obs , temper , output, obsInd){
  sigma2<-par[1] # Variance Parameter
  extremeOutput<-output[obsInd]
  llhd<-temper*sum(dnorm(x=obs, mean = extremeOutput , sd= sqrt(sigma2), log = TRUE)) # Compute Likelihood
  output<-list(extremeOutput,output)
  return(list(llhd,output))
}

for(jobNum in 1:nrow(parMat)){
  jobPar<-parMat[jobNum,]
  # Compute the full likelihood (temper=1) to generated the weights
  llhd_t<-logLikelihood_temper(par =jobPar, obs = obs ,  temper = 1 , output = modelOutput[,jobNum] , obsInd=obsInd)
  save(jobPar,llhd_t,file=paste("/glade/scratch/sanjib/runA_rep_var/output/PF_",cycle,"_",jobNum,".RData",sep=""))
}

# Tempering values
temperVal<-list()
temperVal$cumulative<-0
temperVal$incremental<-0
save(temperVal,file="output_rep_var/temperVal_0.RData")

# Generate Covariance matrix for proposal
# 1. Compute Weights
llhd<-apply(modelOutput, 2, function(x){sum(dnorm(x = obs, mean = x[obsInd] , sd = 1825.215, log = TRUE))})
# 2. Choose the weights that are in the top 10%
covIndex<-which(llhd>quantile(llhd,probs = 0.9))
keepParMat<-parMat[covIndex,-1]
# Generate scaling matrix for the parameter set
CovMat<-cov(keepParMat)*((2.38^2)/ncol(keepParMat)) # Use acceptance ratio from Rosenthal et al. 
CovMat<-CovMat+diag(ncol(CovMat))*(0.001*diag(CovMat)) #Add white noise to ensure positive-definiteness
save(CovMat, file="output_rep_var/BeginCovMat_Tr.RData") # Save covariance matrix