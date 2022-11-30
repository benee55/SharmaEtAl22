# # Copyright (C) 2020 Ben S. Lee
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
# We remove the variance parameter for the lGM

library(invgamma);
library(mvtnorm);
# Parameter Names
parNames<-c("PCTIM" , "ADIMP" , "Q0CHN" , "QMCHN")
boundMat<-rbind(c(0, 5) , # PCTIM 0.3=original maximum
                c(0 , 2), # ADIMP 0.5=original maximum
                c(0.5,4.5), # rutpix_Q0CHN
                c(0.3,1.9)) # rutpix_QMCHN Use Original: 3.4 ; BPrior: 2.25

rownames(boundMat)<-parNames
colnames(boundMat)<-c("lower","upper")  
priorPar<-rbind(c(0.01,0.01), # Inverse Gamma Hyperparameters
                boundMat)
##################################################################################################################
# Observations
# load("/glade/u/home/sanjib/FamosHydroModel/lowDim/input/obsData.RData")
##################################################################################################################
# Priors 
logPrior<-function(par , priorPar){
  lpriorVal<-sum(dunif(x = par[-1] , min = priorPar[-1,1] , max = priorPar[-1,2] , log = TRUE),
                 dinvgamma( x = par[1] , shape = priorPar[1,1] , rate = priorPar[1,2] , log = TRUE ))
}

# Likelihood
logLikelihood_temper<-function(par, obs , j , inputDir , outputDir , temper){
  # Add TRYCatch here 
  
  output<-modelEval(par=par, j=j , inputDir=inputDir , outputDir=outputDir) # Evaluate Model and Obtain output
  sigma2<-par[1] # Variance Parameter
  lenOut<-length(output[[1]])
  llhd<-temper*sum(dnorm(x=obs[1:lenOut], mean = output[[1]] , sd= sqrt(sigma2), log = TRUE)) # COmpute Likelihood
  return(list(llhd,output))
}
  
# Posterior 
logPosterior_temper<-function(par , priorPar , obs , inputDir , outputDir , j , temper){
  lPri<-logPrior( par = par , priorPar = priorPar )
  if(lPri==-Inf){
    return(-1e10)
  }else{
    llhd<-logLikelihood_temper(par = par, obs = obs , j = j , inputDir = inputDir , outputDir = outputDir , temper=temper)
    lpostVal<-lPri+llhd[[1]]
    output<-llhd[[2]]
    return(list(lpostVal,output))
  }
}
