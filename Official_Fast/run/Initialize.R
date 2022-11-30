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
setwd("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/")

for(i in 1:14){
  load(paste("precalibration/output/preCalibrationResults",i,".RData",sep=""))
  if(i==1){
    modelOutput<-outputMat[-nrow(outputMat),]
  }else{
    modelOutput<-cbind(modelOutput,outputMat[-nrow(outputMat),])
  }
}

# This provides all values in time series. Need to index
# sampleIndex<-sample(1:ncol(modelOutput), ens)
# modelOutput<-modelOutput[,sampleIndex]

# Parameters
load("precalibration/output/mhParameters_0.RData")
parMat<-parMat[1:ncol(modelOutput),]
save(parMat,file="output/mhParameters_0.RData")

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
  llhd_t<-logLikelihood_temper(par =jobPar, obs = obs ,  temper = 1 , output = modelOutput[,jobNum] , obsInd=obsInd)
  save(jobPar,llhd_t,file=paste("/glade/scratch/sanjib/runA/output/PF_",cycle,"_",jobNum,".RData",sep=""))
}

# Tempering values
temperVal<-list()
temperVal$cumulative<-0
temperVal$incremental<-0
save(temperVal,file="output/temperVal_0.RData")

# Generate Covariance matrix for proposal
if(FALSE){
llhd<-apply(modelOutput, 2, function(x){sum(dnorm(x = obs[1:2], mean = x[obsInd[1:2]] , sd = 1934.148, log = TRUE))})
}else{
  llhd<-apply(modelOutput, 2, function(x){sum(dnorm(x = obs, mean = x[obsInd] , sd = 2000, log = TRUE))})
}
# 
covIndex<-which(llhd>quantile(llhd,probs = 0.9))
keepParMat<-parMat[covIndex,-1]
# Use acceptance ratio from Rosenthal et al. 
CovMat<-cov(keepParMat)*((2.38^2)/ncol(keepParMat))
CovMat<-CovMat+diag(ncol(CovMat))*(0.001*diag(CovMat))
save(CovMat, file="output/BeginCovMat_Tr.RData")


# Save Bhattacharrya Distance for Exact Case
# ADD VALUES for Streamflow
# meanvar<-0.5379547;sdvar<-0.1259628 #Values from the test runs
# dat1<-rnorm(ens,mean=meanvar,sd=sdvar)
# empirBattaBS<-vector("numeric")
# breaks<-200
# for(j in 1:1000){
#   dat2<-rnorm(ens,mean=meanvar,sd=sdvar)
#   empirBattaBS[j]<-battaDistance(dat.current=dat1,dat.previous=dat2,breaks=breaks)
#   bdMCRange<-quantile(empirBattaBS,probs = c(0.025,0.5,0.975),na.rm = TRUE)
#   }
# save(bdMCRange,file="~/work/hydro/Official_Fast/output/Empirical_BD.RData")
