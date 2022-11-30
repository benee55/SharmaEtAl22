# CHeck to make sure everything passess through windows


rm(list = ls())

source("../../run_discr/mcmc_source.R")


contributionList<-list()
for(j in 1:3){
  contributionList[[j]]<-list()
  for(i in 1:10){
    contributionList[[j]][[i]]<-list()
  } 
  }
parList<-list();for(i in 1:3){parList[[i]]<-list()} 
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

verFileList<-c("InitialRun","Speed10","NoSpeed")
for (j in 1:3){
  setwd(paste("~/Box Sync/LeeEtal-PSUICE-3D/calibration/MPITest/discr_output/",verFileList[j],"/",sep=""))
  m.files<-list.files(pattern="mhParameters")[c(1,2,5:12,3,4)]
  c.files<-c("rsParameters_1.RData",list.files(pattern="contribution")[c(1,4:11,2,3)])
  
for(i in 1:length(m.files)){
  load(m.files[i])
  
  parMat[,10]<-log(parMat[,10],base = 3)
  parMat[,logVect[-5]]<-log(parMat[,logVect[-5]],base = 10)
  parList[[j]][[i]]<-parMat
  load(c.files[i])
  contributionList[[j]][[1]][[i]]<-sle-contributionMat[,1]
  contributionList[[j]][[2]][[i]]<-sle-contributionMat[,2]
  contributionList[[j]][[3]][[i]]<-contributionMat[,6]-contributionMat[,3]
  contributionList[[j]][[4]][[i]]<-contributionMat[,4]/1e+15
  contributionList[[j]][[5]][[i]]<-contributionMat[,5]/1e+12
  contributionList[[j]][[6]][[i]]<-contributionMat[,6]-contributionMat[,8]
  contributionList[[j]][[7]][[i]]<-contributionMat[,6]-contributionMat[,9]
  contributionList[[j]][[8]][[i]]<-contributionMat[,6]-contributionMat[,10]
  contributionList[[j]][[9]][[i]]<-contributionMat[,6]-contributionMat[,11]
  contributionList[[j]][[10]][[i]]<-contributionMat[,6]-contributionMat[,12]
}
}

LogparPollard<-parPollard
LogparPollard[logVect]<-log(LogparPollard[logVect],base=baseVect[logVect])

pdf(file="../FinalRunComparison.pdf",height = 8.5,width=11)
par(mfrow=c(4,3),mar=c(2,2,2,2))


for(k in 1:11){#Parameters
      plot(density(parList[[1]][[1]][,k]),col="black",
           ylim=range(density(parList[[1]][[1]][,k])$y,
                      density(parList[[1]][[12]][,k])$y,
                      density(parList[[2]][[12]][,k])$y,
                      density(parList[[3]][[12]][,k])$y),
           xlim=range(density(parList[[1]][[1]][,k])$x,
                      density(parList[[1]][[12]][,k])$x,
                      density(parList[[2]][[12]][,k])$x,
                      density(parList[[3]][[12]][,k])$x),
           main=parNames[k])
      if(k==1){
        legend("topleft",legend=c("prior","SpeedLim=5","SpeedLim=10","No Speed Lim"),col=c("black","blue","red","darkgreen"),
               lty=c(1,1,1,1),cex=1)  
      }
lines(density(parList[[1]][[12]][,k]),col="blue") 
lines(density(parList[[2]][[12]][,k]),col="red") 
lines(density(parList[[3]][[12]][,k]),col="darkgreen") 
abline(v=LogparPollard[k],lwd=2,col="black",lty=2)
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
for(k in 1:length(contributionList[[1]])){
  
      plot(density(contributionList[[1]][[k]][[1]],na.rm=TRUE),col="black",lwd=2,
           ylim=range(density(contributionList[[1]][[k]][[1]],na.rm=TRUE)$y,
                      density(contributionList[[1]][[k]][[12]],na.rm=TRUE)$y,
                      density(contributionList[[2]][[k]][[12]],na.rm=TRUE)$y,
                      density(contributionList[[3]][[k]][[12]],na.rm=TRUE)$y),
           xlim=range(density(contributionList[[1]][[k]][[1]],na.rm=TRUE)$x,
                      density(contributionList[[1]][[k]][[12]],na.rm=TRUE)$x,
                      density(contributionList[[2]][[k]][[12]],na.rm=TRUE)$x,
                      density(contributionList[[3]][[k]][[12]],na.rm=TRUE)$x,
                      cbind(windows,rep(0.1,2),rep(0.1,2),rep(0.1,2),rep(0.1,2),rep(0.1,2))[,k]),
           main=obsNames[k])
      if(k==1){
        legend("topleft",legend=c("prior","SpeedLim=5","SpeedLim=10","No Speed Lim"),col=c("black","blue","red","darkgreen"),
               lty=c(1,1,1,1),cex=1)  
      }
  if(k%in%1:5){abline(v=windows[,k],col="darkgray",lwd=2,lty=2)}
  
      lines(density(contributionList[[1]][[k]][[12]],na.rm=TRUE),col="blue",lwd=2)    
      lines(density(contributionList[[2]][[k]][[12]],na.rm=TRUE),col="red",lwd=2)    
      lines(density(contributionList[[3]][[k]][[12]],na.rm=TRUE),col="darkgreen",lwd=2)    
  
  
  abline(v=pollardVal[k],lty=2,col="black",lwd=2)
  
}
dev.off()

#check if they pass through windows
minMaxVal<-matrix(NA,nrow=2,ncol=6)
for(j in 1:3){
for(k in 1:6){
  minMaxVal[,k]<-c(min(contributionList[[j]][[k]][[2]],
                       contributionList[[j]][[k]][[3]],
                       contributionList[[j]][[k]][[4]],
                       contributionList[[j]][[k]][[5]],
                       contributionList[[j]][[k]][[6]],
                       contributionList[[j]][[k]][[7]],
                       contributionList[[j]][[k]][[8]],
                       contributionList[[j]][[k]][[9]],
                       contributionList[[j]][[k]][[10]],
                       contributionList[[j]][[k]][[11]],
                       contributionList[[j]][[k]][[12]]),
                   max(contributionList[[j]][[k]][[2]],
                       contributionList[[j]][[k]][[3]],
                       contributionList[[j]][[k]][[4]],
                       contributionList[[j]][[k]][[5]],
                       contributionList[[j]][[k]][[6]],
                       contributionList[[j]][[k]][[7]],
                       contributionList[[j]][[k]][[8]],
                       contributionList[[j]][[k]][[9]],
                       contributionList[[j]][[k]][[10]],
                       contributionList[[j]][[k]][[11]],
                       contributionList[[j]][[k]][[12]]))
  
  
}
  print(minMaxVal)
  print(rbind(minMaxVal[1,-6]>windows[1,],minMaxVal[2,-6]<windows[2,]))
}
#Check if all values are within the windows
print(minMaxVal)
minMaxVal[1,-6]>windows[1,]
minMaxVal[2,-6]<windows[2,]
