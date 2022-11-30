################################################################################################
rm(list=ls())
# Run this when the pre-calibration sample is finished 
# Observations
obs<-read.table("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/SBYP1_obs.txt" , stringsAsFactors = FALSE)
obs<-obs[-(1:2),]
fullDate<-(paste(obs$V2,obs$V3,obs$V4,sep=""))
newDate<-as.Date(fullDate,"%Y%m%d") # Convert to Date
useObs<-as.numeric(obs$V5)*0.0283


# Subset Observations to 2004-2008
obs<-t(apply(obs,1,as.numeric))
subset<-which(obs[,2]%in% c(2004:2008))

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
obs<-obs[obsInd,]
obs<-as.numeric(obs[,5])*0.0283 # Need to convert

# Time Intervals for each run
intervalMat<-rbind(c("20040810","20040921"), # 2004/08/10-2004/09/21
                   c("20041115","20050407"), # 2004/11/15-2005/04/07
                   c("20051101","20051202"), # 2005/11/01-2005/12/02
                   c("20060515","20060701"), # 2006/05/15-2006/07/01
                   c("20061015","20061119"), # 2006/10/15-2006/11/19
                   c("20070101","20070318"), # 2007/01/01-2007/03/18
                   c("20080101","20080311")) # 2008/01/10-2008/03/11
newIntervalDate<-t(apply(intervalMat,1,as.Date,format="%Y%m%d"))
newIntervalDate
# Plot of Observations from 2004-2008 with calibration observations (extremes) in red
par(mfrow=c(1,1))
plot(x=newDate[subset],y=useObs[subset] , typ="p" , ylab="Streamflow" , xlab="Date", pch=16, cex=0.5 , main="Selinsgrove Streamflow 2004-2009")
points(x=newDate[obsInd], y= obs, pch=16, col="red")
abline(h=4950.55, col="blue")
  legend("topright", legend=c("Observations", "Selected Events" , "Action Stage"), pch=c(16,16,NA), lty=c(NA,NA,1),  
         col=c("black","red","blue"))  




################################################################################################


# Examine the Model outputs 
setwd("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/famos/Official_Fast/")
source("run/mcmc_source_Tr.R")
source("run/rWrapper.R")
load("output/mhParameters_0.RData")
inputDir<-"/gpfs/group/kzk10/default/private/hydrocalib/SGrove/famos/Official_Fast/input"
outputDir<-'/gpfs/group/kzk10/default/private/hydrocalib/SGrove/famos/Official_Fast/output'
load("input/obsData.RData")

load("output/preCalibrationResults.RData") # Load randomly selected sample
ens<-length(outputMat)
useDate<-newDate[obsInd] # Use these dates for observations - 21 total
 
# Violin plots of the mode outputs
par(mfrow=c(5,5), mar=c(2,2,2,2))
outputMatOrig<-outputMat
outputMat<-matrix(NA,nrow=ens,ncol=21)
for(i in 1:ens){outputMat[i,]<-outputMatOrig[[i]][-22]}

MSE<-apply(outputMat,1,function(x){mean((x-obs)^2)})
for(i in 1:21){
  vioplot(outputMat[,i], ylim=range(outputMat[,i],obs[i]), main = useDate[i])
  points(x=1, y=obs[i], col="red" ,pch=16)
  abline(h=4950.55, col="blue" ,lwd=2) # ACtion Stage
}


dev.off()
outputMat<-data.frame(outputMat)
colnames(outputMat)<-useDate
par(mar=c(6,6,4,2))
boxplot(outputMat , las=2 , main="Model Output and observations")

for(i in 1:21){
  points(x=i, y=obs[i], col="red" ,pch=16)
}
abline(h=4950.55, col="blue" ,lwd=2) # ACtion Stage

par(mar=c(6,4,4,2))
vioplot(outputMat  , las=2 , main="Model Output and observations" , pch.col="blue")
for(i in 1:21){
  points(x=i, y=obs[i], col="red" ,pch=16 , cex=2)
}
abline(h=4950.55, col="blue" ,lwd=2) # ACtion Stage
legend("topright", legend=c("Observation" , "Action Stage" , "Ensemble Mean"), pch=c(16,16,16),  
       col=c("red","blue","white"))  
####################################################################################################
# How we can score each model run. 

# Log Likelihood function
llhd_function<-function(output,sigma2, obs){
  sum(dnorm(x = obs, mean = output , sd = sqrt(sigma2), log = TRUE))
}

# Calculate Weights for 0.1 tempering
llhd<-apply(outputMat, 1, function(x){sum(dnorm(x = obs, mean = x , sd = 1000, log = TRUE))})
fullWeight<- llhd# first weights are log likelihood
weights<-fullWeight*0.1
weights<-exp(weights-max(weights))
weights<-weights/sum(weights)

# Resample particles 
reSampID<-sample(x=1:ens,size = ens, prob= weights, replace = TRUE)
length(unique(reSampID)) # Number of unique particles 
keepParMat<-parMat[reSampID,] # Resampled particles


# Density plots of each parameter
par(mfrow=c(4,3), mar=c(2,2,2,2))
for(i in 2:13){
  plot(density(keepParMat[,i]) , xlim=range(boundMat[i-1,1:2]),  main=parNames[i])
  abline(v=boundMat[i-1,1:2], col="red")
}

# Violin plots of each model output
keepDat<-outputMat[reSampID,]
par(mfrow=c(4,5), mar=c(2,2,2,2))
for(i in 1:19){
  vioplot(keepDat[,i], ylim=range(keepDat[,i],obs[i],4950.55), main = useDate[i])
  points(x=1, y=obs[i], col="red" ,pch=16)
  abline(h=4950.55, col="blue" ,lwd=2) # ACtion Stage
}
