# Review LHS output's depth 
##==============================================================================
##
rm(list=ls())
library(hydroGOF)
## Script for converting streamflow data to Gage Height using 
## the USGS rating curve
main_path=setwd('/gpfs/group/kzk10/default/private/hydrocalib/SGrove/')
infile<-paste('/gpfs/group/kzk10/default/private/hydrocalib/SGrove/output/SBYP1_discharge_outlet.ts',sep='')
big<-scan(infile,skip=165,list(gauge_no=0,day=0,sim_flow=0))
sim_flow<-big$sim_flow # Observed Streamflow
# Rating Curve for Selinsgrove
ratingcurve<-paste('/gpfs/group/kzk10/default/private/hydrocalib/SGrove/ratingcurve_selinsgrove.txt',sep='')
big<-scan(ratingcurve,skip=1,list(height=0,flow=0))
flow<-big$flow*0.0283  # Streamflow
height<-big$height # Height 
# Approximate Guage Height for observed streamflow
estimated_gage_height=approx(flow,height,xout=sim_flow)$y # Approximate the gauge height
predictDepth<-function(flow,height,simFlow){approx(x=flow,y=height,xout=simFlow)$y} # Function to convert streamflow to depth
###########################################################################################
# Observation Data
load("LHSoutput/input/LHS_Inputs.RData")
dat<-read.table("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/SBYP1_obs.txt",stringsAsFactors = FALSE)
obsFull<-dat[(30*365+10):(30*365+10+365-1),]
obs<-as.numeric(obsFull[,5])*0.0283 # Observations
obs_height=approx(flow,height,xout=obs)$y # Approximate the gauge height
# LHS Input
lpostVect<-foo<-vector("numeric")
flowMat<-matrix(NA,nrow=365,ncol=1000)
for(i in 1:1000){
  print(i)
  load(paste("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/LHSoutput/output/output",i,".RData",sep=""))
  foo[i]<-1 ; 
  flowMat[,i]<-lpostVal[[2]][,2] ;
  lpostVect[i]<-lpostVal[[1]]
}
# Save Date, Parameters, streamflow, and depth
flowDepth<-list()
flowDepth$date<-lpostVal[[2]][,1]
flowDepth$parameters<-parMat
flowDepth$flow<-t(flowMat)
flowDepth$depth<-t(apply(flowDepth$flow,1,predictDepth, flow=flow , height=height))

save.image(file="SummaryData.RData")
########################################################################################################################################################
# Criteria to restrict Model Runs
maxID<-which.max(obs)
nseVect<-apply(flowDepth$depth,1,NSE , obs = estimated_gage_height , na.rm = TRUE) # NSE
kgeVect<-apply(flowDepth$depth,1,KGE , obs = estimated_gage_height , na.rm = TRUE) # KGE
# Criteria #1: Depth>20
# Criteria #2: NSE>0.5
# Criteria #3: KGE>0.5
keepID<-which(flowDepth$depth[,maxID]>20 & nseVect>0.5 & kgeVect>0.5)



# Visualize
par(mfrow=c(1,2), mar=c(2,2,2,2))
plot(x=1:365,y=obs,typ="l",lwd=1.5 , ylab="Streamflow" , xlab="Date (2011)" , 
     main="Pre-calibration Streamflow (2011)" , ylim=range(flowDepth$flow[keepID,],flowDepth$flow[visInd,],na.rm = TRUE))
for(i in keepID){
  lines(x=1:365,y=flowDepth$flow[i,],col="gray",lwd=1.5)
}
lines(x=1:365,y=obs,col="black",lwd=2)
# legend("topleft" , legend=c("Observations (1000 Total)" , "Model Output"),
       # lty=c(1,1) , col=c("black","gray"),lwd=rep(2,2),cex=1.5)

# everything
visInd<-round(seq(1,1000,length.out = 100))
plot(x=1:365,y=obs,typ="l",lwd=1.5 , ylab="Streamflow" , xlab="Date (2011)" , 
     main="All Streamflow (2011)" , ylim=range(flowDepth$flow[visInd,],na.rm = TRUE))
for(i in visInd){
  lines(x=1:365,y=flowDepth$flow[i,],col="gray",lwd=1.5)
}
lines(x=1:365,y=obs,col="black",lwd=2)
legend("topleft" , legend=c("Observations (1000 Total)" , "Model Output"),
       lty=c(1,1) , col=c("black","gray"),lwd=rep(2,2),cex=1.5)



##################################################################################################################################
group1<-visInd
d1<-density(flowDepth$depth[group1,maxID])
group2<-keepID
d2<-density(flowDepth$depth[group2,maxID])
# Narrow plot based on likelihood
dev.off()
plot(density(flowDepth$depth[,maxID]),ylim=range(d1$y,d2$y) , main="Depth at 9/9/11", col="red",
     xlab="depth")
lines(density(flowDepth$depth[group2,maxID]),col="blue")
abline(v=estimated_gage_height[maxID] , col="gray",lwd=2)
abline(v=obs_height[maxID] , col="black",lwd=2,lty=2)
legend("topleft" , legend=c("Prior Samples" , "Pre-calibration" ,  "Hand-tuned 'Best Guess'","Observed Depth") , 
       lty=c(1,1,1,2) , col=c("red","blue","gray","black"),cex=1.5 , lwd=2)



########################################################################
parNames<-c("PCTIM" ,"ADIMP" , "UZTWM" , "LZTWM" , 
            "LZFSM" , "LZFPM" , "LZSK" , "UZFWM" , 
            "REXP" , "UZK" , "Q0CHN" , "QMCHN")
priorPar<-rbind(c(0 , 0.3), # PCTIM
                c(0 , 0.5), # ADIMP
                c(-1.1 , -0.2), # UZTWM
                c(-2.5 , -0.25), # LZTWM
                c(-2.8 , -0.2), # LZFSM
                c(-3.2 , -0.3), # LZFPM
                c(-3.8 , -0.1), # LZSK
                c(-4.5 , -0.1), # UZFWM
                c(-3.5 , -0.1), # REXP
                c(-3.5 , -0.1), # UZK
                c(0.5,4.5), # rutpix_Q0CHN
                c(0.3,2.25)) # rutpix_QMCHN


par(mfrow=c(3,4),mar=c(2,2,2,2))
for(i in 1:nrow(priorPar)){
  plot(density(parMat[group2,i]),main=parNames[i],xlim=range(parMat[,i]))
  abline(v=priorPar[i,],col="red")
}

