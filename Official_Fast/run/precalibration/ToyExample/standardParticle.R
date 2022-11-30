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
####################################################################################################

rm(list=ls())
setwd("~/work/LeeEtal-PSUICE-3D/calibration/AOAS_code/ToyExample/")
source("source_IS.R")
load(file="mcmcGoldStandard.RData")
##########################################################
library(fields);library(classInt);library(mvtnorm);library(invgamma)
library(snow);library(doParallel);library(foreach);library(adaptMCMC)

# Step 1: Draw samples from Prior 
ens=2000
pfIter=200
type = "MPI" # PSOCK or MPI
cores =359

initPar<-cbind(rnorm(ens,mean=prior1[1],sd=sqrt(prior2[1])),
                    runif(ens,min=prior1[2],max=prior2[2]),
                    rinvgamma(ens,shape=prior1[3],scale=prior2[3]),
                    rinvgamma(ens,shape=prior1[4],scale=prior2[4]))

tempSched<-rep(0.1,10)
cumsum(tempSched)
testDat<-pSpider(ens=ens,
        initPar=initPar,
        cores=cores,
        type=type,
        obs=Z,
        Xvar=X,
        distMat=distMat,
        covFun=expCov,
        prior1=prior1,
        prior2=prior2,
        tempSched=tempSched,
        pfIter=pfIter,
        parTruth=parTruth,
        mcmcComp=chain1,
        TuningMat=TuningMat)






save.image(file="standardParticle.RData")



