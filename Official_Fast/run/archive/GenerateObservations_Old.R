rm(list=ls())
# Observations
obs<-read.table("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/SBYP1_obs.txt" , stringsAsFactors = FALSE)
obs<-obs[-(1:2),]
obs<-t(apply(obs,1,as.numeric))
obsInd<-matrix(NA, nrow= 5, ncol=2)
obsInd[1,]<-c(which(obs[,2]%in%c(2004) & obs[,3]==8 & obs[,4]==15), # 2004/08/15-2004/09/21
              which(obs[,2]%in%c(2004) & obs[,3]==9 & obs[,4]==21))

obsInd[2,]<-c(which(obs[,2]%in%c(2004) & obs[,3]==12 & obs[,4]==15), # 2004/12/15-2005/04/07
              which(obs[,2]%in%c(2005) & obs[,3]==4 & obs[,4]==7))

obsInd[3,]<-c(which(obs[,2]%in%c(2006) & obs[,3]==5 & obs[,4]==15), # 2006/05/15-2006/06/30
              which(obs[,2]%in%c(2006) & obs[,3]==6 & obs[,4]==30))

obsInd[4,]<-c(which(obs[,2]%in%c(2007) & obs[,3]==2 & obs[,4]==15), # 2007/02/15-2007/03/18
              which(obs[,2]%in%c(2007) & obs[,3]==3 & obs[,4]==18))

obsInd[5,]<-c(which(obs[,2]%in%c(2008) & obs[,3]==1 & obs[,4]==10), # 2008/01/10-2008/03/11
              which(obs[,2]%in%c(2008) & obs[,3]==3 & obs[,4]==11))

fullObsInd<-c(obsInd[1,1]:obsInd[1,2],
              obsInd[2,1]:obsInd[2,2],
              obsInd[3,1]:obsInd[3,2],
              obsInd[4,1]:obsInd[4,2],
              obsInd[5,1]:obsInd[5,2])
obs<-obs[fullObsInd,]
obs<-as.numeric(obs[,5])*0.0283 # Need to convert

save(obs, file="/gpfs/group/kzk10/default/private/hydrocalib/SGrove/famos/Official_Fast/input/obsData.RData")

plot.ts(obs)