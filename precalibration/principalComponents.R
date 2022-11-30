rm(list=ls())
setwd('/gpfs/group/kzk10/default/private/hydrocalib/SGrove/')
load("SummaryData.RData")
load("CalibrationData3D_Full.RData")
source("~/work/PICAR/rev_Code/comparison/sharedFunctions.R")
library(mlegp)

# Remove Mean
M<-flowDepth$flow # p x n
meanVect<-as.numeric(t(M)%*%rep(1/nrow(M),nrow(M))) # Mean
centerdM<-M - matrix(meanVect,nrow=nrow(M) , ncol= ncol(M),byrow=TRUE)

svdVal<-svd(centerdM)
Jy<-which.min(abs(cumsum(svdVal$d)/sum(svdVal$d)-0.6)) ; print(Jy)
K<-svdVal$v[,1:Jy]%*%diag(sqrt(svdVal$d[1:Jy])) # n x p
Y_i<-t(solve(t(K)%*%K)%*%t(K)%*%t(centerdM))

U<-svdVal$v[,1:Jy]
D_half<-diag(sqrt(svdVal$d[1:Jy]))

Z_proj<-(U%*%D_half)%*%t(Y_i)

p1<-Z_proj[,2]+meanVect
p2<-centerdM[2,]+meanVect
err<-p1-p2
sum(err^2)/365


plot(x=1:365, y=p1 , typ="l")
lines(x=1:365, y=p2 , col="red")



# Observations
svdVal<-svd(centerdM)
Jy<-which.min(abs(cumsum(svdVal$d)/sum(svdVal$d)-0.65)) ; print(Jy)
K<-svdVal$v[,1:Jy]%*%diag(sqrt(svdVal$d[1:Jy])) # n x p
U<-svdVal$v[,1:Jy]
D_half<-diag(sqrt(svdVal$d[1:Jy]))
Y_i<-t(solve(t(K)%*%K)%*%t(K)%*%t(centerdM))

projZ<-t(solve(t(K)%*%K)%*%t(K)%*%(obs-meanVect))
p1<-obs
p2<-(U%*%D_half)%*%t(projZ)+meanVect
err<-p1-p2
sum(err^2)/365
plot(x=1:365, y=p1 , typ="l")
lines(x=1:365, y=p2 , col="red")


plot(x=1:365, y=err , typ="l")



######## WHICH REDUCED SPACE COEEFS

redY<-apply(Y_i,1,function(x)sum((x-projZ)^2))
closeInd<-which.min(redY)

plot(x=1:365, y=obs , typ="l" ,ylim=range(flowDepth$flow[closeInd,],obs))
lines(x=1:365, y=flowDepth$flow[closeInd,] , col="red")
sum((obs-flowDepth$flow[closeInd,])^2)/365


basis<-(U%*%D_half)
plot.ts(basis[,6])
