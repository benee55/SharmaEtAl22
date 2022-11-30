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
obs<-as.numeric(dat[-(1:2),5])*0.0283 # Observations
date<-as.numeric(dat[-(1:2),2]) # Observations
obs_height=approx(flow,height,xout=obs)$y 
plot.ts(obs_height, main="Depth")
keepInd<-which(obs_height>=20 & date%in%c(2004:2013))
length(keepInd)
useDat<-dat[-(1:2),]
refDat<-cbind(useDat[keepInd,],obs_height[keepInd])
table(refDat$V2)

library(xtable)
print(xtable(refDat[,2:6]))
