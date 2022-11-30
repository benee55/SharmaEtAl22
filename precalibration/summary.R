# Code to combine LHS_output
##==============================================================================
##
rm(list=ls())
## Script for converting streamflow data to Gage Height using 
## the USGS rating curve
main_path=setwd('/gpfs/group/kzk10/default/private/hydrocalib/SGrove/')
infile<-paste('/gpfs/group/kzk10/default/private/hydrocalib/SGrove/output/SBYP1_discharge_outlet.ts',sep='')
big<-scan(infile,skip=165,list(gauge_no=0,day=0,sim_flow=0))
sim_flow<-big$sim_flow
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 203.4   368.6   844.4  1346.4  1826.2  8497.5 

ratingcurve<-paste('/gpfs/group/kzk10/default/private/hydrocalib/SGrove/ratingcurve_selinsgrove.txt',sep='')
big<-scan(ratingcurve,skip=1,list(height=0,flow=0))
height<-big$height
flow<-big$flow*0.0283
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 54.96  2289.18  6958.63  8329.58 13713.51 22063.92 


estimated_gage_height=approx(flow,height,xout=sim_flow)$y
plot(x=flow,y=height , pch=16 , col="red")
points(x=sim_flow,y=estimated_gage_height , pch=16)


# Function to convert streamflow to height
# approx(flow,height,xout=sim_flow)$y

###########################################################################################
# LHS Input
load("LHSoutput/input/LHS_Inputs.RData")
dat<-read.table("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/SBYP1_obs.txt",stringsAsFactors = FALSE)
obsFull<-dat[(30*365+10):(30*365+10+365-1),]


obs<-as.numeric(obsFull[,5])*0.0283
summary(obs)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 117.2   348.1   984.8  1478.8  1966.8 11404.9 
# Rating Curve
# 54.96  2289.18  6958.63  8329.58 13713.51 22063.92 

# LHS Input
lpostVect<-foo<-vector("numeric")
flowMat<-matrix(NA,nrow=365,ncol=1000)
for(i in 1:1000){
  print(i)
  load(paste("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/LHSoutput/output/output",i,".RData",sep=""))
  foo[i]<-1 ; 
  flowMat[,i]<-lpostVal[[2]][,2] ;
  lpostVect[i]<-lpostVal[[1]]
  # bar = try(load(paste("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/LHSoutput/output/output",i,".RData",sep="")), silent = TRUE)
  # if (inherits(bar, "try-error")){foo[i]<-0
  # }else{
  #   foo[i]<-1 ; flowMat[,i]<-lpostVal[[2]][,2] ;lpostVect[i]<-lpostVal[[1]]}
}


flowDepth<-list()
flowDepth$date<-lpostVal[[2]][,1]
flowDepth$parameters<-parMat
flowDepth$flow<-t(flowMat)

summary(as.numeric(flowDepth$flow))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 6.69   377.75   781.38  1424.45  1750.83 38302.12 


predictDepth<-function(flow,height,simFlow){approx(x=flow,y=height,xout=simFlow)$y}
bar<-apply(flowDepth$flow,1,predictDepth, flow=flow , height=height)

flowDepth$depth<-t(bar)
# save(flowDepth ,file = "/gpfs/group/kzk10/default/private/hydrocalib/SGrove/LHSoutput/LHS_depth.RData")

######################################

obsMat<-matrix(obs,nrow=nrow(flowDepth$flow),ncol=ncol(flowDepth$flow),byrow = TRUE)
diffMat<-(obsMat-flowDepth$flow)^2
mseMat<-apply(diffMat^2,1,mean)
quantileVect<-quantile(mseMat,probs = c(0.05,0.1,0.25,0.5))

plot(x=1:365,y=obs,typ="l",lwd=1.5 , ylab="Streamflow" , xlab="Date (2011)" , 
     main="Streamflow (2011)")
minID<-which.min(mseMat)
lines(x=1:365,y=flowDepth$flow[minID,],col="blue",lwd=1.5)
lines(x=1:365,y=sim_flow,col="red",lwd=1.5)
legend("topleft" , legend=c("Observations" , "Hand-tuned" , "Least Squares"),
       lty=c(1,1,1) , col=c("black","red" ,"blue"),lwd=rep(2,3) , cex=1.5)

visInd<-round(seq(1,1000,length.out = 100))
plot(x=1:365,y=obs,typ="l",lwd=1.5 , ylab="Streamflow" , xlab="Date (2011)" , 
     main="Streamflow (2011)" , ylim=range(flowDepth$flow[visInd,],na.rm = TRUE))
for(i in visInd){
  lines(x=1:365,y=flowDepth$flow[i,],col="gray",lwd=1.5)
}
lines(x=1:365,y=obs,col="black",lwd=2)
legend("topleft" , legend=c("Observations (1000 Total)" , "Model Output"),
       lty=c(1,1) , col=c("black","gray"),lwd=rep(2,2),cex=1.5)



maxID<-which.max(obs)
group1<-which(mseMat<quantileVect[1])
d1<-density(bar[maxID,group1])
group2<-which(mseMat<quantileVect[2])
d2<-density(bar[maxID,group2])
flowDepth10percent<-flowDepth
names(flowDepth10percent)
flowDepth10percent$parameters<-flowDepth10percent$parameters[group2,]
flowDepth10percent$flow<-flowDepth10percent$flow[group2,]
flowDepth10percent$depth<-flowDepth10percent$depth[group2,]
save(obs,flowDepth ,flowDepth10percent ,file = "/gpfs/group/kzk10/default/private/hydrocalib/SGrove/LHSoutput/LHS_depth.RData")

group3<-which(mseMat<quantileVect[3])
d3<-density(bar[maxID,group3])
group4<-which(mseMat<quantileVect[4])
d4<-density(bar[maxID,group4])

d0<-density(bar[maxID,])
flowDepth$date[maxID]
# Narrow plot based on likelihood
plot(density(bar[maxID,]),ylim=range(d0$y,d2$y) , main="Depth at 9/9/11", col="red",
     xlab="depth")
# lines(density(bar[maxID,group1]),col="blue")
lines(density(bar[maxID,group2]),col="blue")
abline(v=estimated_gage_height[maxID] , col="black",lwd=2)
legend("topleft" , legend=c("Prior Samples" , "Pre-calibration" , "Hand-tuned 'Best Guess'") , 
       lty=c(1,1,1) , col=c("red","blue","black"),cex=1.5 , lwd=2)

########################################################################
# PRe0calibration
load("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/LHSoutput/LHS_depth.RData")
#HydroGoaf
maxID<-which.max(obs)
dateMax<-flowDepth$date[maxID]; dateMax
max26Depth<-flowDepth$depth[,maxID] # Flow is greater than 26
over26<-which(max26Depth>26)
max26Depth<-max26Depth[over26]

NSE<-function(x,obs)({
  removeID<-which(is.na(x))
  meanObs<-mean(obs)
  1- (sum((x-obs)^2)/sum((x-meanObs)^2))
})

foo<-apply(flowDepth$flow,1, NSE, obs=obs)

thresh5<-which(foo>0.5)
thresh6<-which(foo>0.6)
# thresh7<-which(foo>0.7)
# 
# thresh6<-thresh6[which(flowDepth$depth[thresh6,maxID]>26)]
# thresh7<-thresh7[which(flowDepth$depth[thresh7,maxID]>26)]

par(mfrow=c(1,2), mar=c(2,2,2,2))
visInd<-round(seq(1,1000,length.out = 100))
plot(x=1:365,y=obs,typ="l",lwd=1.5 , ylab="Streamflow" , xlab="Date (2011)" , 
     main="Streamflow (2011)" , ylim=range(flowDepth$flow[visInd,],na.rm = TRUE))
for(i in thresh6){
  lines(x=1:365,y=flowDepth$flow[i,],col="gray",lwd=1.5)
}
lines(x=1:365,y=obs,col="black",lwd=2)
legend("topleft" , legend=c("Observations (1000 Total)" , "Model Output"),
       lty=c(1,1) , col=c("black","gray"),lwd=rep(2,2),cex=1.5)

# everything
visInd<-round(seq(1,1000,length.out = 100))
plot(x=1:365,y=obs,typ="l",lwd=1.5 , ylab="Streamflow" , xlab="Date (2011)" , 
     main="Streamflow (2011)" , ylim=range(flowDepth$flow[visInd,],na.rm = TRUE))
for(i in visInd){
  lines(x=1:365,y=flowDepth$flow[i,],col="gray",lwd=1.5)
}
lines(x=1:365,y=obs,col="black",lwd=2)
legend("topleft" , legend=c("Observations (1000 Total)" , "Model Output"),
       lty=c(1,1) , col=c("black","gray"),lwd=rep(2,2),cex=1.5)


##################
maxID<-which.max(obs)
group1<-thresh5
d1<-density(bar[maxID,group1])
group6<-thresh6
d6<-density(bar[maxID,group6])

dev.off()
# Narrow plot based on likelihood
plot(density(bar[maxID,]),ylim=range(d0$y,d1$y,d2$y) , main="Depth at 9/9/11", col="red",
     xlab="depth")
# lines(density(bar[maxID,group1]),col="blue")
lines(density(bar[maxID,group1]),col="blue")
lines(density(bar[maxID,group2]),col="green")
abline(v=estimated_gage_height[maxID] , col="black",lwd=2)
legend("topleft" , legend=c("Prior Samples" , "Pre-calibration 0.5" , "Pre-calibration 0.6", "Hand-tuned 'Best Guess'") , 
       lty=c(1,1,1,1) , col=c("red","blue","green","black"),cex=1.5 , lwd=2)


# NSE doesn't require tuning a variance parameter 
# NSE does not account for temporal correlations 




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
  plot(density(parMat[group1,i]),main=parNames[i],xlim=range(parMat[,i]))
  abline(v=priorPar[i,],col="red")
}



plot(x=1:365,y=obs,typ="l")
for(k in group1){lines(x=1:365 , y=flowDepth$flow[k,],col="gray")}
plot(density(parMat[group1,]))
priorPar





