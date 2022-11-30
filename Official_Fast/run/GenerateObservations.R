rm(list=ls())
# Observations
obs<-read.table("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/input/SBYP1_obs.txt" , stringsAsFactors = FALSE)
obs<-obs[-(1:2),]
obs<-t(apply(obs,1,as.numeric))
# intervalMat<-c("20030601T00","20080331T23") # 2003/06/01-2008/03/31
startInd<-which(obs[,2]==2003 & obs[,3]==6 & obs[,4]==1)
endInd<-which(obs[,2]==2008 & obs[,3]==3 & obs[,4]==31)
obs<-obs[startInd:endInd,]
obs[,5]<-as.numeric(obs[,5])*0.0283 # Need to convert
# Subset Extremes
obsInd<-c(which(obs[,2]==2004 & obs[,3]==9 & obs[,4]%in%c(19,20)),
          which(obs[,2]==2005 & obs[,3]==1 & obs[,4]%in%c(15,16)),
          which(obs[,2]==2005 & obs[,3]==3 & obs[,4]%in%c(30:31)),
          which(obs[,2]==2005 & obs[,3]==4 & obs[,4]%in%c(3:6)),
          which(obs[,2]==2005 & obs[,3]==12 & obs[,4]%in%c(1)),
          which(obs[,2]==2006 & obs[,3]==6 & obs[,4]%in%c(28:30)),
          which(obs[,2]==2006 & obs[,3]==11 & obs[,4]%in%c(18)),
          which(obs[,2]==2007 & obs[,3]==3 & obs[,4]%in%c(16:17)),
          which(obs[,2]==2008 & obs[,3]==2 & obs[,4]%in%c(8)),
          which(obs[,2]==2008 & obs[,3]==3 & obs[,4]%in%c(6,9,10))
)
# Final reformatting
obs<-obs[,-1] #1766 total observations
fullObs<-obs[,4]
obs<-fullObs[obsInd]
save(obs, obsInd,fullObs, file="/glade/u/home/sanjib/FamosHydroModel/Official_Fast/input/obsData.RData")
