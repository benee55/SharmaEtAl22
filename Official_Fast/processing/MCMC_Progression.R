rm(list=ls())
setwd("~/work/LeeEtal-PSUICE-3D/calibration/Official_Fast/output/")
source("../run/mcmc_source_Tr.R")

# Check Number of cycles
output2Contribution<-function(x){
  sle<-57.554
  output<-c(sle-x[[2]]["plio",  "sle"],
            sle-x[[2]]["lig",  "sle"],
            x[[2]]["mod", "sle"]- x[[2]]["lgm", "sle"],
            x[[2]]["mod", "total"]/1e+15,
            x[[2]]["mod", "grounded"]/1e+12,
            x[[2]]["mod", "sle"]- x[[2]]["2100", "sle"],
            x[[2]]["mod", "sle"]- x[[2]]["2200", "sle"],
            x[[2]]["mod", "sle"]- x[[2]]["2300", "sle"],
            x[[2]]["mod", "sle"]- x[[2]]["2400", "sle"],
            x[[2]]["mod", "sle"]- x[[2]]["2500", "sle"])
}
# combine MCMC chains

mhProgression<-function(cycle){
  
  fileDirLoad<-list.files(pattern =  paste("MCMC_",cycle,"_",1,"_",sep=""))
  ens<-length(fileDirLoad)
  orderNum<-sapply(fileDirLoad,function(x) as.numeric(strsplit(x,split="_|[/.]")[[1]][4])) # need to order it by numerical JobNum
  fileDirLoad<-fileDirLoad[order(orderNum)]
  
  parList<-list()
  resultsList<-list()
  
for(i in 1:length(fileDirLoad)){
  load(fileDirLoad[i])
  parList[[i]]<-amcmc.out[[1]]
  resultsList[[i]]<-t(sapply(amcmc.out[[5]],output2Contribution))
}
 
  summaryParList<-list()
  summaryResultsList<-list()
  for(h in 1:ncol(parList[[1]])){
    summaryParList[[h]]<-t(sapply(parList,function(x) x[,h]))
    if(h!=11){
      summaryResultsList[[h]]<-t(sapply(resultsList,function(x)x[,h]))  
    }
    
  }
  
  return(list(parameters=summaryParList,results=summaryResultsList))
}

if(FALSE){
  # Save
  mhProgressionList<-list()
  for(cycle in 1:2){
    mhProgressionList[[cycle]]<-mhProgression(cycle=cycle)
  }
  save(mhProgressionList,file="mhProgression.RData")
  
}



# Processing
load("mhProgression.RData")
sampInd<-round(seq(1,15,length.out = 5))
par(mfrow=c(2,3),mar=c(2,2,2,2))
for(i in 6:10){
  d1<-density(mhProgressionList[[1]][[2]][[i]][,sampInd[5]])
  d2<-density(mhProgressionList[[2]][[2]][[i]][,sampInd[5]])
  
  plot(d1,main=2000+(i-5)*100,col=1,ylim=range(d1$y,d2$y),xlim=range(d1$x,d2$x))
  lines(d2,col=2);
}


################################################################################################
hpd <- function(samp,p=0.05)
{
  ## to find an approximate (1-p)*100% HPD interval from a
  ## given posterior sample vector samp
  
  r <- length(samp)
  samp <- sort(samp)
  rang <- matrix(0,nrow=trunc(p*r),ncol=3)
  dimnames(rang) <- list(NULL,c("low","high","range"))
  for (i in 1:trunc(p*r)) {
    rang[i,1] <- samp[i]
    rang[i,2] <- samp[i+(1-p)*r]
    rang[i,3] <- rang[i,2]-rang[i,1]
  }
  hpd <- rang[order(rang[,3])[1],]
  hpd<-c(hpd[1],mean(samp),hpd[2])
  names(hpd)<-c("low","mean","high")
  return(hpd)
}

################################################################################################
FinalParMatList<-list()
  for(k in 1:11){
    for(i in 1:2){
      if(i==1){
        FinalParMatList[[k]]<-apply(mhProgressionList[[i]][[1]][[k]],2,hpd,p=0.05)    
      }else{
        FinalParMatList[[k]]<-cbind(FinalParMatList[[k]],apply(mhProgressionList[[i]][[1]][[k]],2,hpd,p=0.05)  )
      }
    }
    if(k%in%c(1,2,4,9,11)){FinalParMatList[[k]]<-log(FinalParMatList[[k]],10)}
    if(k==10){FinalParMatList[[k]]<-log(FinalParMatList[[k]]/0.3,10)}
  }
  

FinalProjectionList<-list()
for(k in 1:10){
  for(i in 1:2){
    if(i==1){
      FinalProjectionList[[k]]<-apply(mhProgressionList[[i]][[2]][[k]],2,hpd,p=0.05)    
    }else{
      FinalProjectionList[[k]]<-cbind(FinalProjectionList[[k]],apply(mhProgressionList[[i]][[2]][[k]],2,hpd,p=0.05)  )
    }
  }
}


pdf(file="../processing/ProjectionsParameterHPD.pdf",height=8.5,width=11)
#PLot Projections
lengthVect<-c(15,15)
endVect<-cumsum(lengthVect)
startVect<-endVect-(lengthVect-1)
#Parameters
par(mfrow=c(3,4),mar=c(4,2,2,2))
for(i in 1:length(FinalParMatList)){
  x<-FinalParMatList[[i]]
  plot(x=1:endVect[2],y=x[2,],typ="n",main=parNames[i],ylim=range(x),
       xlab="Iterations",ylab="Value")
  for(k in 1:2){
    xPlot<-startVect[k]:endVect[k]
    lines(x=xPlot,y=x[2,xPlot],col=k,lwd=2)
    lines(x=xPlot,y=x[1,xPlot],col=k,lwd=1,lty=1)
    lines(x=xPlot,y=x[3,xPlot],col=k,lwd=1,lty=1)
    abline(v=startVect[k],col=k)
  }
}

#Projections
projNames<-c("Pliocene","LIG","LGM","Modern Volume","Modern Area","2100","2200","2300","2400","2500")
par(mfrow=c(2,3),mar=c(4,2,2,2))
for(i in 1:length(FinalProjectionList)){
  if(i==6){par(mfrow=c(2,3),mar=c(4,2,2,2))}
  x<-FinalProjectionList[[i]]
  plot(x=1:endVect[2],y=x[2,],typ="n",main=projNames[i],ylim=range(x),
       xlab="Iterations",ylab="Value")
  for(k in 1:2){
    xPlot<-startVect[k]:endVect[k]
    lines(x=xPlot,y=x[2,xPlot],col=k,lwd=2)
    lines(x=xPlot,y=x[1,xPlot],col=k,lwd=1,lty=1)
    lines(x=xPlot,y=x[3,xPlot],col=k,lwd=1,lty=1)
    abline(v=startVect[k],col=k)
  }
}
dev.off()  

##############################################
FullParMatList<-list()
for(k in 1:11){
  for(i in 1:2){
    if(i==1){
      FullParMatList[[k]]<-mhProgressionList[[i]][[1]][[k]]
    }else{
      FullParMatList[[k]]<-cbind(FullParMatList[[k]],mhProgressionList[[i]][[1]][[k]])
    }
  }
  if(k%in%c(1,2,4,9,11)){FullParMatList[[k]]<-log(FullParMatList[[k]],10)}
  if(k==10){FullParMatList[[k]]<-log(FullParMatList[[k]]/0.3,10)}
}

FullProjectionMatList<-list()
for(k in 1:10){
  for(i in 1:2){
    if(i==1){
      FullProjectionMatList[[k]]<-mhProgressionList[[i]][[2]][[k]]
    }else{
      FullProjectionMatList[[k]]<-cbind(FullProjectionMatList[[k]],mhProgressionList[[i]][[2]][[k]])
    }
  }
}

# Function to evaluate Battacharrya Distance 
library(snow);library(doParallel);library(foreach);
mp_type = "PSOCK" # PSOCK or MPI
nprocs =3
cl <- parallel::makeCluster(nprocs, type=mp_type)
doParallel::registerDoParallel(cl)
meanTest<-0
sdTest<-1
useBreaks<-200
test<-rnorm(1007,mean=meanTest,sd=sdTest)
pt<-proc.time()
testMat  <- foreach::foreach(ens.member=1:10000,
                                .combine='c') %dopar% {
                                  battaDistance(dat.current=rnorm(1007,mean=meanTest,sd=sdTest),dat.previous=test,breaks=useBreaks)
                                }
proc.time()-pt
MCErr<-c(quantile(testMat,probs = 0.025),mean(testMat),quantile(testMat,probs = 0.975))
MCErr
########################################
bdMat<-list()
cl <- parallel::makeCluster(nprocs, type=mp_type)
doParallel::registerDoParallel(cl)
pt<-proc.time()
bdMat[[1]]  <- foreach::foreach(i=1:length(FullParMatList),
                                .combine='rbind') %dopar% {
                                  if(i%in%c(1,2,4,9,11)){
                                    finalSamp<-log(mhProgressionList[[2]][[1]][[i]][,15],base = 10)
                                  }else if(i==10){
                                    finalSamp<-log(mhProgressionList[[2]][[1]][[i]][,15]/0.3,base = 10)
                                  }else{
                                    finalSamp<-mhProgressionList[[2]][[1]][[i]][,15] 
                                    }
                                   
                                  apply(FullParMatList[[i]],2,battaDistance,dat.previous=finalSamp,breaks=useBreaks)
                                }

bdMat[[2]]  <- foreach::foreach(i=1:length(FullProjectionMatList),
                                .combine='rbind') %dopar% {
  finalSamp<-mhProgressionList[[2]][[2]][[i]][,15]  
  apply(FullProjectionMatList[[i]],2,battaDistance,dat.previous=finalSamp,breaks=useBreaks)
                             }

proc.time()-pt



#PLot Projections
lengthVect<-c(15,15)
endVect<-cumsum(lengthVect)
startVect<-endVect-(lengthVect-1)

#Parameters
pdf(file="../processing/fullBD.pdf",height=8.5,width=11)
par(mfrow=c(3,4),mar=c(4,4,2,2))
for(i in 1:nrow(bdMat[[1]])){
  x<-bdMat[[1]][i,]
  plot(x=1:endVect[2],y=x,typ="n",main=parNames[i],ylim=range(x),
       ylab="Bhattacharyya Distance", xlab="Iteration")
  for(k in 1:2){
    xPlot<-startVect[k]:endVect[k]
    lines(x=xPlot,y=x[xPlot],col=k,lwd=2)
    lines(x=xPlot,y=x[xPlot],col=k,lwd=2,lty=3)
    lines(x=xPlot,y=x[xPlot],col=k,lwd=2,lty=3)
    abline(v=startVect[k],col=k)
    
  }
  abline(h=MCErr[3],lty=2,col="black",lwd=2)
}


#Projections
projNames<-c("Pliocene","LIG","LGM","Modern Volume","Modern Area","2100","2200","2300","2400","2500")
par(mfrow=c(2,3),mar=c(4,4,2,2))
for(i in 1:nrow(bdMat[[2]])){
  if(i==6){par(mfrow=c(2,3),mar=c(4,4,2,2))}
  x<-bdMat[[2]][i,]
  plot(x=1:endVect[2],y=x,typ="n",main=projNames[i],ylim=range(x),
       ylab="Bhattacharyya Distance", xlab="Iteration")
  for(k in 1:2){
    xPlot<-startVect[k]:endVect[k]
    lines(x=xPlot,y=x[xPlot],col=k,lwd=2)
    lines(x=xPlot,y=x[xPlot],col=k,lwd=2,lty=3)
    lines(x=xPlot,y=x[xPlot],col=k,lwd=2,lty=3)
    abline(v=startVect[k],col=k)
    
  }
  abline(h=MCErr[3],lty=2,col="black",lwd=2)
}


dev.off()
# Plot Density
projNames<-c("Pliocene","LIG","LGM","Modern Volume","Modern Area","2100","2200","2300","2400","2500")
par(mfrow=c(2,5),mar=c(2,2,2,2))
for(i in 1:length(FullProjectionMatList)){
  x<-FullProjectionMatList[[i]]
  d1<-density(x[,endVect[2]])
  d2<-density(x[,startVect[2]+5])
  d3<-density(x[,startVect[2]+10])
  plot(d1,col="black",main=projNames[i],
       ylim=range(d1$y,d2$y,d3$y),
       xlim=range(d1$x,d2$x,d3$x),lwd=2)
  lines(d2,col="red",lwd=2)
  lines(d3,col="blue",lwd=2)
  abline(v=hpd(x[,endVect[2]])[-2],lty=2)
  abline(v=hpd(x[,startVect[2]+5])[-2],lty=2,col="red")
  abline(v=hpd(x[,startVect[2]+10])[-2],lty=2,col="blue")
}

##############################################################################################################################
# Check Bhattacharrya for every last in each cycle vs. Rest 
#Cycle 1 
lengthVect<-c(15,15)
endVect<-cumsum(lengthVect)
startVect<-endVect-(lengthVect-1)
startVect<-startVect

pdf(file="../processing/BD_Projections.pdf",height=8.5,width=11)
for(i in 1:10){
  par(mfrow=c(1,2))
  for(k in 1:2){
  print(i)
  finalSamp<-FullProjectionMatList[[i]][,endVect[k]]
  bd<-apply(FullProjectionMatList[[i]][,startVect[k]:endVect[k]],2,battaDistance,dat.previous=finalSamp,breaks=useBreaks)
  plot.ts(bd,main=paste(projNames[i],"cycle",k),xlab="Iteration",ylab="Bhattacharyya Distance",ylim=range(bd,0,MCErr[3]))
  abline(h=MCErr[3],col="Red")
}
}
dev.off()


pdf(file="../processing/BD_Parameters.pdf",height=8.5,width=11)
for(i in 1:11){
  print(i)
  par(mfrow=c(1,2))
  for(k in 1:2){
    finalSamp<-FullParMatList[[i]][,endVect[k]]
    bd<-apply(FullParMatList[[i]][,startVect[k]:endVect[k]],2,battaDistance,dat.previous=finalSamp,breaks=useBreaks)
    plot.ts(bd,main=paste(parNames[i],"cycle",k),xlab="Iteration",ylab="Bhattacharyya Distance",ylim=range(bd,0,MCErr[3]))
    abline(h=MCErr[3],col="Red")
  }
}
dev.off()

########################################################################################################################
# Number of unique particles 
pdf(file="../processing/UniqueParticles.pdf",height=8.5,width=11)
foo<-apply(FullParMatList[[1]],2,function(x)length(unique(x)))
stopPoint<-14
plot(1:endVect[2],y=foo,typ="n",ylab="Unique Particles",xlab="Iteration",
     main=paste("Unique Particles vs. Iteration \n After Iteration",stopPoint+1),
     ylim=c(1000,2015))
for(i in 1:2){
  xPlot<-startVect[i]:endVect[i]
  lines(x=xPlot,y=foo[xPlot],col=i)
  abline(v=startVect[i]+stopPoint,lty=2,col=i)
  text(x=startVect[i]+stopPoint,y=foo[startVect[i]+stopPoint],col=i,labels=foo[startVect[i]+stopPoint])
}
dev.off()
