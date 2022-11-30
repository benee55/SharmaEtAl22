# # Copyright (C) 2021 Ben S. Lee
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
rm(list=ls())
setwd("~/Dropbox/hydroFamos/run/precalibration")
source("test_Source.R")
inputDir<-"~/Dropbox/hydroFamos/run/precalibration/input"
outputDir<-"~/Dropbox/hydroFamos/run/precalibration/output"

#Initial Parameters
ensembleN<-100000
parMat<-apply(boundMat[-1,],1,function(x,ens){runif(ens,min=x[1],max=x[2])},ens=ensembleN)
outputMat<-t(apply(parMat,1,compModel, x=inputX))

SSE<-apply(outputMat,1,function(x,obs){sum((x-obs)^2)}, obs=obs)
LowSSE<-which(SSE<quantile(SSE,probs = 0.05))
CovMat<-cov(parMat[LowSSE,])*((2.38^2)/ncol(parMat))
# load("MCMCOutput.RData")
# CovMat<-cov(mcmcOutput[,1:2])*((2.38^2)/ncol(mcmcOutput))
save(CovMat, file="output/BeginCovMat_Tr.RData")
