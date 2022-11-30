# # Copyright (C) 2020 Ben S. Lee
#email: slee287@gmu.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# We remove the variance parameter for the lGM

library(invgamma);library(lhs);library(mvtnorm);library(tmvtnorm)
# Parameter Names
parNames<-c("PCTIM" , "ADIMP" , "UZTWM" ,"LZTWM" , "LZFSM" , "LZFPM" , "LZSK" , "UZFWM" , "REXP" , "UZK" , "Q0CHN" , "QMCHN")
boundMat<-rbind(c(0 , 0.3), # PCTIM
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
                c(0.3,2.25)) # rutpix_QMCHN Use 2.25 instead of 3.4

rownames(boundMat)<-parNames
colnames(boundMat)<-c("lower","upper")  
parNames<-c("S2","PCTIM" , "ADIMP" , "UZTWM" ,"LZTWM" , "LZFSM" , "LZFPM" , "LZSK" , "UZFWM" , "REXP" , "UZK" , "Q0CHN" , "QMCHN")

priorPar<-rbind(c(0.001,0.001), # Inverse Gamma Hyperparameters
                boundMat)
##################################################################################################################
# Observations
load("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/famos/Official_Fast/input/obsData.RData")
##################################################################################################################
# Load Source for RWrapper and Forward Model 
source("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/famos/Official_Fast/run/rWrapper.R")
##################################################################################################################
# Priors 
logPrior<-function(par , priorPar){
  lpriorVal<-sum(dunif(x = par[-1] , min = priorPar[-1,1] , max = priorPar[-1,2] , log = TRUE),
                 dinvgamma( x = par[1] , shape = priorPar[1,1] , rate = priorPar[1,2] , log = TRUE ))
}

# Likelihood
logLikelihood_temper<-function(par, obs , j , inputDir , outputDir , temper){
  writeInput( par = par[-1] , j = j , dir = inputDir) # Function 1: Write Input 
  writeOutput( j = j , dir = outputDir) # Function 2: Create Directory
  output<-runReadModel( j = j, inputDir = inputDir , outputDir = outputDir )
  sigma2<-par[1]
  llhd<-temper*sum(dnorm(x=obs, mean = output[,2] , sd= sqrt(sigma2), log = TRUE))
  return(list(llhd,output))
}

# Posterior 
logPosterior_temper<-function(par , priorPar , obs , inputDir , outputDir , j , temper){
  lPri<-logPrior( par = par , priorPar = priorPar )
  if(lPri==-Inf){
    return(-1e10)
  }else{
    llhd<-logLikelihood_temper(par = par, obs = obs , j = j , inputDir = inputDir , outputDir = outputDir , temper=temper)
    lpostVal<-lPri+llhd[[1]]
    output<-llhd[[2]]
    return(list(lpostVal,output))
  }
}
#######################################################################################################
# Function to Calculate New Likelihood with respect to tempering 
calcPF<-function(cycle,jobNum,llhdTemper,mcmcTemper,initResults){
  par<-parMat[jobNum,]
  ################################################################################################
  names(par)<-parNames
  ################################################################################################
  lprior<-logPrior(par=par)
  llhd<-((initResults[[1]]-lprior)/mcmcTemper)*llhdTemper
  results<-initResults[[2]]
  return(list(llhd,results))
}
#############################################################################

mcmcManual_tempered<-function(iter,
                              init,
                              propCov,
                              inputDir , 
                              outputDir ,
                              obs, 
                              priorPar,
                              jobNum,
                              temper,
                              llhdTemper,
                              initResults,
                              parNames,
                              additional=FALSE){

  # Initialize Values
  n=length(obs)
  names(init)<-parNames
  # Intialize containers
  alphaMat<-parMat<-candMat<-matrix(NA,nrow=iter,ncol=length(init))
  candMat[1,]<-parMat[1,]<-init ; alphaMat[1,]<-1
  colnames(parMat)<-colnames(candMat)<-colnames(alphaMat)<-parNames
  # Add posterior
  lpost<-vector("numeric");
  if(additional==FALSE){
    lpost[1]<-temper*(initResults[[1]]/llhdTemper)+logPrior(par=init) # Calculate Log Posterior from previous run  - IS step
  }else{
    lpost[1]<-initResults[[1]] # Use log posterior from previous run
  }
  # Note that we are doing this (MCMC Weight- Cumulative) x (Prev)/likelihoodweight + Log Prior, where Prev=(Log Likelihood *likelihoodweight)
  resultsList<-list();resultsList[[1]]<-initResults[[2]]
  curResultsList<-list();curResultsList[[1]]<-initResults[[2]]
  #Begin MCMC
  for(i in 2:iter){
   # Propose on the log scale from the mutlivariate truncated normal
    candMat[i,]<-c(parMat[i-1,1] , rtmvnorm(n = 1,mean = parMat[i-1,-1],
                                            sigma = propCov,
                                            lower = boundMat[1,],upper = boundMat[3,]) )

      
    curResults<-logPosterior_temper(par=as.numeric(candMat[i,]),
                                    priorPar = priorPar,
                                    obs = obs,
                                    inputDir = inputDir,
                                    outputDir = outputDir,
                                    j=jobNum , 
                                    temper=temper)
        ############################################################################################################################
    # account for asymmetric proposal
    # Do some simulations to make sure this working properly
    boundMat<-priorPar[-1,]
    adj1<-log(pmvnorm(upper= boundMat[2,],lower=boundMat[1,],mean=as.numeric(parMat[i-1,-1]),sigma=propCov))[1]#previous
    adj2<-log(pmvnorm(upper= boundMat[2,],lower=boundMat[1,],mean=as.numeric(candMat[i,-1]),sigma=propCov))[1]#current 
    alpha<-min(curResults[[1]]+adj1-lpost[i-1]-adj2,log(1)) # Asymmetric Proposal. Follow Example. Add and subtract constants. 
    #https://journal.r-project.org/archive/2010-1/RJournal_2010-1_Wilhelm+Manjunath.pdf
    ############################################################################################################################
    
    alphaMat[i,]<-exp(alpha)
    
    # Accept/Reject
    if(log(runif(n = 1,min = 0,max = 1))<alpha){ # Accept
      parMat[i,]<-candMat[i,];
      lpost[i]<-curResults[[1]]
      resultsList[[i]]<-curResults[[2]]
    }else{ # REject
      parMat[i,]<-parMat[i-1,]
      lpost[i]<-lpost[i-1]
      resultsList[[i]]<-resultsList[[i-1]]
    }
    
    curResultsList[[i]]<-curResults[[2]]
    
    # Update Gamma Parameter - Gibbs Update
    alphaPar<-priorPar[1,1]+0.5*n 
    betaPar<-(sum((resultsList[[i]][,2]-obs)^2)+2*priorPar[1,2])/2
    parMat[i,1]<-rinvgamma(n=1, shape = alphaPar, rate=betaPar)
    
  }
  return(list(parMat,candMat,alphaMat,lpost,resultsList,curResultsList))
}



accRateFunc<-function(x){
  accRate<-(length(unique(x))-1)/(length(x)-1)
  return(accRate)
}
########################################################################################################
# Turn Full WEights to tempered weights and ESS
fullToESS<-function(fullWeight,tempSeq,showESS=FALSE){
  weights<-fullWeight*tempSeq
  weights<-exp(weights-max(weights))
  weights<-weights/sum(weights)
  ESS<-1/sum((weights^2))
  
  if(showESS==TRUE){
    return(ESS)
  }else{
    return(list(weights=weights,ESS=ESS))
  }
}
#Maximize ESS
maxESS<-function(fullWeight,tempSeq,ens){
  ESS<-fullToESS(fullWeight,tempSeq,showESS=TRUE)
  (ESS-0.5*ens)^2
}
# Optimize Weights to achieve criteria
optimizeWeights<-function(weightVect,cumulTemp,ens){
  optimTemp<-optimize(f = maxESS,fullWeight=weightVect,
                      ens=ens,interval = c(0.1,1-cumulTemp))
  useTemp<-optimTemp$minimum;
  useESSWeights<-fullToESS(fullWeight = weightVect,tempSeq = useTemp,showESS = FALSE)
  return(list(useTemp,useESSWeights))
}
########################################################################
# Function to calculate weights + resample
combineIS<-function(cycle,ens,cumulTemp){
  setwd("output/")
  fileDirLoad<-list.files(pattern = paste("PF_",cycle,sep=""))
  orderNum<-sapply(fileDirLoad,function(x) as.numeric(strsplit(x,split="_|[/.]")[[1]][3])) # Need to order it by numerical JobNum
  fileDirLoad<-fileDirLoad[order(orderNum)]
  
  ens<-length(fileDirLoad)
  weightVect<-vector("numeric")
  resultsList<-list()
  # Load Files
  for(i in 1:length(fileDirLoad)){
    load(fileDirLoad[i])
    resultsList[[i]]<-llhd_t[[2]]
    if(i==1){
      k<-length(jobPar)
      parWeightMat<-matrix(NA,nrow=ens,ncol=k)
      parWeightMat[1,]<-jobPar
    }else{parWeightMat[i,]<-jobPar}
    weightVect[i]<-llhd_t[[1]]
  }
  #Optimize
  optimList<-optimizeWeights(weightVect=weightVect,cumulTemp=cumulTemp,ens=ens)
  temperVal<-list();
  temperVal$cumulative<-cumulTemp+optimList[[1]]
  temperVal$incremental<-optimList[[1]]
  save(temperVal,file=paste("temperVal_",cycle,".RData",sep=""))
  weights<-optimList[[2]]$weights
  # Weight 
  reSampleInd<-sample(x=1:ens,size = ens,replace = TRUE,prob = weights)
  # Resample
  parMat<-parWeightMat[reSampleInd,]
  initResultsList<-list(weightVect[reSampleInd],resultsList[reSampleInd])
  #Transformed COntribution Mat
  # Save final files
  save(parMat,parWeightMat,weightVect,weights,reSampleInd,
       initResultsList,optimList,
       file=paste("rsParameters_",cycle,".RData",sep=""))
}
########################################################################

########################################################################
# Function to combine MH output 

combineMH<-function(cycle,ens,stage){

fileDirLoad<-list.files(pattern =  paste("MCMC_",cycle,"_",stage,"_",sep=""))
ens<-length(fileDirLoad)
orderNum<-sapply(fileDirLoad,function(x) as.numeric(strsplit(x,split="_|[/.]")[[1]][4])) # need to order it by numerical JobNum
fileDirLoad<-fileDirLoad[order(orderNum)]

weightVect<-vector("numeric")
resultsList<-list()

# Load Files
for(i in 1:length(fileDirLoad)){
  load(fileDirLoad[i])

  if(i==1){
    parMat<-matrix(NA,nrow=ens,ncol=ncol(amcmc.out[[1]]))
    parMat[1,]<-amcmc.out[[1]][nrow(amcmc.out[[1]]),]
    
    # Total particles
    if(stage==1){ # If first stage create new matrix
      TotalParMat<-amcmc.out[[1]]
    }else{ # just append if later stages
      load(paste("output/totalParticles_",cycle,".RData",sep="")) 
      TotalParMat<-rbind(TotalParMat,amcmc.out[[1]])
    }
  
    acceptVect<-vector("numeric")
    acceptVect[i]<-accRateFunc(amcmc.out[[1]][,1])
    kMCMC<-length(amcmc.out[[5]])
    
  }else{
    parMat[i,]<-amcmc.out[[1]][nrow(amcmc.out[[1]]),]
    TotalParMat<-rbind(TotalParMat,amcmc.out[[1]])
    acceptVect[i]<-accRateFunc(amcmc.out[[1]][,1])
    kMCMC<-length(amcmc.out[[5]])
    
  }
  resultsList[[i]]<-amcmc.out[[5]][[kMCMC]]
  weightVect[i]<-amcmc.out[[4]][[kMCMC]]
}

#Initial Results for next run
initResultsList<-list(weightVect,resultsList)


# Save 
# total particles
save(TotalParMat,
     file=paste("output/totalParticles_",cycle,".RData",sep=""))

# Final Particles
save(parMat,acceptVect,initResultsList,
     file=paste("output/mhParameters_",cycle,".RData",sep=""))
}
########################################################################################################################
# Combine Total Particles Function
combineTotalParticles<-function(cycle){

  if(cycle==1){
    load(paste("output/totalParticles_",cycle,".RData",sep=""))
    masterTotalParticles<-TotalParMat
  }else{
    load("output/totalParticles_1.RData")
    masterTotalParticles<-TotalParMat
    for(i in 2:cycle){# Loop Through all Total particles files
      load(paste("output/totalParticles_",i,".RData",sep=""))
      masterTotalParticles<-rbind(masterTotalParticles,TotalParMat)  
    }
    
  }
  save(masterTotalParticles,file="output/masterTotalParticles.RData")
}

########################################################################
# Bhattacharryya Distance
# Count Table Function
countTable<-function(x,dat1,dat2){
  d1<-sum(dat1>=x[1]&dat1<x[2],na.rm = TRUE)/length(dat1)
  d2<-sum(dat2>=x[1]&dat2<x[2],na.rm = TRUE)/length(dat2)
  BC<-(sqrt(d1*d2))
}
# Function to evaluate Battacharrya Distance 
battaDistance<-function(dat.current,dat.previous,breaks){
  
  penEnd<-quantile(dat.previous,probs = c(0.005,0.995),na.rm = TRUE)
  partition<-c(-Inf,seq(penEnd[1],penEnd[2],length.out = breaks-2),Inf)
  partitionTable<-rbind(partition[-length(partition)],partition[-1])
  
  sampBD<--log(sum(apply(partitionTable,2,countTable,
                            dat1=dat.previous,
                            dat2=dat.current)))
  return(sampBD)
}
########################################################################
# Function to generate covariance matrix
genPropMat<-function(cycle,scale){
  load(paste("output/rsParameters_",cycle,".RData",sep=""))
  
  # Covariance Matrix
  if(cycle==1){
    load("output/BeginCovMat_Tr.RData")
    CovMat<-CovMat*scale # Optimal proposal from Rosenthal et al. (2008)
  }else{
    load("output/masterTotalParticles.RData")
    TotalParMat<-rbind(masterTotalParticles,parMat)
    uniqueID<-!duplicated(TotalParMat[,1])
    TotalParMat<-TotalParMat[uniqueID,]
    colnames(TotalParMat)<-parNames
    #####################################################################################################################
    CovMat<-cov(TotalParMat)*((2.38^2)/ncol(TotalParMat))*scale # Optimal proposal from Rosenthal et al. (2008)
  }
   return(CovMat)
}




