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


#Initial Parameters
parMat<-apply(boundMat[-1,],1,function(x,ens){runif(ens,min=x[1],max=x[2])},ens=ensembleN)
parMat<-cbind(rinvgamma(ensembleN,shape = priorPar[1,1], rate = priorPar[1,2]),parMat)

save(parMat,file="output/mhParameters_0.RData")
temperVal<-list()
temperVal$cumulative<-0
temperVal$incremental<-0
save(temperVal,file="output/temperVal_0.RData")

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
# save(bdMCRange,file="output/Empirical_BD.RData")
