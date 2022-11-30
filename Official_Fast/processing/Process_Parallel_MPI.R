# CHeck to make sure everything passess through windows

setwd("~/Box Sync/LeeEtal-PSUICE-3D/calibration/MPITest/discr_output/InitialRun/")
rm(list = ls())

source("../../run_discr/mcmc_source.R")
m.files<-list.files(pattern="mhParameters")
c.files<-c("rsParameters_1.RData",list.files(pattern="contribution"))
parList<-list()
contributionList<-list();for(i in 1:10){contributionList[[i]]<-list()}

parNames<-c( "OCFACMULT",
             "OCFACMULTASE",
             "CALVNICK",
             "CRHSHELF",
             "TAUASTH",
             "CALVLIQ",
             "CLIFFVMAX",
             "FACEMELTRATE",
             "ENHANCESHEET",
             "ENHANCESHELF",
             "crhfac")
baseVect<-rep(NA,11)
logVect<-c(1,2,4,9,10,11)
baseVect[logVect]<-c(10,10,10,10,3,10)
sle<-57.554
for(i in 1:length(m.files)){
  load(m.files[i])
  
  parMat[,10]<-log(parMat[,10],base = 3)
  parMat[,logVect[-5]]<-log(parMat[,logVect[-5]],base = 10)
  parList[[i]]<-parMat
  load(c.files[i])
  contributionList[[1]][[i]]<-sle-contributionMat[,1]
  contributionList[[2]][[i]]<-sle-contributionMat[,2]
  contributionList[[3]][[i]]<-contributionMat[,6]-contributionMat[,3]
  contributionList[[4]][[i]]<-contributionMat[,4]/1e+15
  contributionList[[5]][[i]]<-contributionMat[,5]/1e+12
  contributionList[[6]][[i]]<-contributionMat[,6]-contributionMat[,8]
  contributionList[[7]][[i]]<-contributionMat[,6]-contributionMat[,9]
  contributionList[[8]][[i]]<-contributionMat[,6]-contributionMat[,10]
  contributionList[[9]][[i]]<-contributionMat[,6]-contributionMat[,11]
  contributionList[[10]][[i]]<-contributionMat[,6]-contributionMat[,12]
  }


LogparPollard<-parPollard
LogparPollard[logVect]<-log(LogparPollard[logVect],base=baseVect[logVect])

pdf(file="FullRuns.pdf",height = 8.5,width=11)
par(mfrow=c(4,3),mar=c(2,2,2,2))
for(k in 1:11){
for(j in 1:length(m.files)){
  if(j==1){
    plot(density(parList[[j]][,k]),col=j,
         ylim=range(density(parList[[1]][,k])$y,
                    density(parList[[2]][,k])$y,
                    density(parList[[3]][,k])$y,
                    density(parList[[4]][,k])$y,
                    density(parList[[5]][,k])$y,
                    density(parList[[6]][,k])$y,
                    density(parList[[7]][,k])$y,
                    density(parList[[8]][,k])$y,
                    density(parList[[9]][,k])$y,
                    density(parList[[10]][,k])$y,
                    density(parList[[11]][,k])$y),
         xlim=range(density(parList[[1]][,k])$x,
                    density(parList[[2]][,k])$x,
                    density(parList[[3]][,k])$x,
                    density(parList[[4]][,k])$x,
                    density(parList[[5]][,k])$x,
                    density(parList[[6]][,k])$x,
                    density(parList[[7]][,k])$x,
                    density(parList[[8]][,k])$x,
                    density(parList[[9]][,k])$x,
                    density(parList[[10]][,k])$x,
                    density(parList[[11]][,k])$x),
         main=parNames[k])
  }
  lines(density(parList[[j]][,k]),col=j)
  abline(v=LogparPollard[k],lwd=2,col="red",lty=2)
}
}

###########################################################################
obsNames<-c("Pliocene","LIG","LGM","Modern Volume","Modern Area","SLR2100","SLR2200","SLR2300","SLR2400","SLR2500")
lower.window<-c(5,3,-20,24.5,10.8)
upper.window<-c(30,8,-4,29.5,13.8)
load("pollardResults.RData")
sle<-57.554
pollardVal<-c(sle-resultsPollard[[2]][1,3],
              sle-resultsPollard[[2]][2,3],
              resultsPollard[[2]][4,3]-resultsPollard[[2]][3,3],
              resultsPollard[[2]][4,1]/1e+15,
              resultsPollard[[2]][4,2]/1e+12,
              resultsPollard[[2]][4,3]-resultsPollard[[2]][6,3],
              resultsPollard[[2]][4,3]-resultsPollard[[2]][7,3],
              resultsPollard[[2]][4,3]-resultsPollard[[2]][8,3],
              resultsPollard[[2]][4,3]-resultsPollard[[2]][9,3],
              resultsPollard[[2]][4,3]-resultsPollard[[2]][10,3])


windows<-matrix(c(lower.window,upper.window),nrow=2,ncol=5,byrow = TRUE)
par(mfrow=c(4,3),mar=c(2,2,2,2))
for(k in 1:length(contributionList)){
  for(j in 1:length(m.files)){
    if(j==1){
      plot(density(contributionList[[k]][[j]],na.rm=TRUE),col=j,lwd=2,
           ylim=range(density(contributionList[[k]][[1]],na.rm=TRUE)$y,
                      density(contributionList[[k]][[2]],na.rm=TRUE)$y,
                      density(contributionList[[k]][[3]],na.rm=TRUE)$y,
                      density(contributionList[[k]][[4]],na.rm=TRUE)$y,
                      density(contributionList[[k]][[5]],na.rm=TRUE)$y,
                      density(contributionList[[k]][[6]],na.rm=TRUE)$y,
                      density(contributionList[[k]][[7]],na.rm=TRUE)$y,
                      density(contributionList[[k]][[8]],na.rm=TRUE)$y,
                      density(contributionList[[k]][[9]],na.rm=TRUE)$y,
                      density(contributionList[[k]][[10]],na.rm=TRUE)$y,
                      density(contributionList[[k]][[11]],na.rm=TRUE)$y),
           xlim=range(density(contributionList[[k]][[1]],na.rm=TRUE)$x,
                      density(contributionList[[k]][[2]],na.rm=TRUE)$x,
                      density(contributionList[[k]][[3]],na.rm=TRUE)$x,
                      density(contributionList[[k]][[4]],na.rm=TRUE)$x,
                      density(contributionList[[k]][[5]],na.rm=TRUE)$x,
                      density(contributionList[[k]][[6]],na.rm=TRUE)$x,
                      density(contributionList[[k]][[7]],na.rm=TRUE)$x,
                      density(contributionList[[k]][[8]],na.rm=TRUE)$x,
                      density(contributionList[[k]][[9]],na.rm=TRUE)$x,
                      density(contributionList[[k]][[10]],na.rm=TRUE)$x,
                      density(contributionList[[k]][[11]],na.rm=TRUE)$x,
                      cbind(windows,rep(0.1,2),rep(0.1,2),rep(0.1,2),rep(0.1,2),rep(0.1,2))[,k]),
           main=obsNames[k])
    }
    lines(density(contributionList[[k]][[j]],na.rm=TRUE),col=j,lwd=2)
  }
  if(k%in%1:5){abline(v=windows[,k],col="blue",lwd=2,lty=2)}
  abline(v=pollardVal[k],lty=1,col="red",lwd=2)
  
}
dev.off()
#check if they pass through windows
minMaxVal<-matrix(NA,nrow=2,ncol=6)
for(k in 1:6){
  minMaxVal[,k]<-c(min(contributionList[[k]][[2]],
      contributionList[[k]][[3]],
      contributionList[[k]][[4]],
      contributionList[[k]][[5]],
      contributionList[[k]][[6]],
      contributionList[[k]][[7]],
      contributionList[[k]][[8]],
      contributionList[[k]][[9]],
      contributionList[[k]][[10]],
      contributionList[[k]][[11]]),
      max(contributionList[[k]][[2]],
           contributionList[[k]][[3]],
           contributionList[[k]][[4]],
           contributionList[[k]][[5]],
           contributionList[[k]][[6]],
           contributionList[[k]][[7]],
          contributionList[[k]][[8]],
          contributionList[[k]][[9]],
          contributionList[[k]][[10]],
          contributionList[[k]][[11]]))
      
  
}
#Check if all values are within the windows
print(minMaxVal)
minMaxVal[1,-6]>windows[1,]
minMaxVal[2,-6]<windows[2,]
