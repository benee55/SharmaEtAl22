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

# Libraries ---------------------------------------------------------------
library(invgamma);
library(mvtnorm);
library(tmvtnorm)

# Parameters --------------------------------------------------------------
parNames<-c("PCTIM" , "ADIMP" , "UZTWM" ,"LZTWM" , 
            "LZFSM" , "LZFPM" , "LZSK" , "snow_SCF" ,
            "REXP" , "UZK" , "Q0CHN" , "QMCHN")

boundMat<-rbind(c(0, 5) , # PCTIM 0.3=original maximum
                c(0 , 2), # ADIMP 0.5=original maximum
                c(-50 , -0.1), # UZTWM
                c(-70 , -0.1), # LZTWM
                c(-100 , -0.1), # LZFSM
                c(-100 , -0.1), # LZFPM Old -120
                c(-3.8 , -0.1), # LZSK
                c(0.5 , 1.5), # snow_SCF
                c(-3.5 , -0.1), # REXP
                c(-3.5 , -0.1), # UZK
                c(0.5,4.5), # rutpix_Q0CHN
                c(0.3,1.9)) # rutpix_QMCHN Use Original: 3.4 ; BPrior: 2.25 

rownames(boundMat)<-parNames
colnames(boundMat)<-c("lower","upper")  
parNames<-c("S2",
            "PCTIM" , "ADIMP" , "UZTWM" ,"LZTWM" , 
            "LZFSM" , "LZFPM" , "LZSK" , "snow_SCF" , 
            "REXP" , "UZK" , "Q0CHN" , "QMCHN")

originalBoundMat<-boundMat
rep2orig<-function(par){ # Convert from reparameterized to original parameters
  c(par[1],par[-1]*(originalBoundMat[,2]-originalBoundMat[,1])/10+originalBoundMat[,1])
}

orig2rep<-function(par){ # Convert from original parameters to reparameterized
  c(par[1],10*(par[-1]-originalBoundMat[,1])/(originalBoundMat[,2]-originalBoundMat[,1]))
}

priorPar<-rbind(c(0.01,0.01), # Inverse Gamma Hyperparameters
                cbind(rep(0,12),rep(10,12)))


##################################################################################################################
# Observations
load("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/input/obsData.RData")
##################################################################################################################
# Load Source for RWrapper and Forward Model 
if(FALSE){
  source("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/run_rep/rWrapper_Continuous_Fast.R")  
}else{
  source("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/run_rep/rWrapper_Continuous.R")  
}
##################################################################################################################
# Priors 
logPrior<-function(par , priorPar){
  lpriorVal<-sum(dunif(x = par[-1] , min = priorPar[-1,1] , max = priorPar[-1,2] , log = TRUE),
                 dinvgamma( x = par[1] , shape = priorPar[1,1] , rate = priorPar[1,2] , log = TRUE ))
}

# Likelihood
logLikelihood_temper<-function(par, obs , j , inputDir , outputDir , temper){
  # Add TRYCatch here 
  
  output<-modelEval(par=par, j=j , inputDir=inputDir , outputDir=outputDir) # Evaluate Model and Obtain output
  sigma2<-par[1] # Variance Parameter
  lenOut<-length(output[[1]])
  llhd<-temper*sum(dnorm(x=obs[1:lenOut], mean = output[[1]] , sd= sqrt(sigma2), log = TRUE)) # COmpute Likelihood
  return(list(llhd,output))
}
  
# Posterior 
logPosterior_temper<-function(par , priorPar , obs , inputDir , outputDir , j , temper){
  
  lPri<-logPrior( par = par , priorPar = priorPar )
  if(lPri==-Inf){
    return(-1e10)
  }else{
    usePar<-rep2orig(par) # Reparameterize
    llhd<-logLikelihood_temper(par = usePar, obs = obs , j = j , inputDir = inputDir , outputDir = outputDir , temper=temper)
    lpostVal<-lPri+llhd[[1]]
    output<-llhd[[2]]
    return(list(lpostVal,output))
  }
}
#######################################################################################################
# Function to Calculate New Likelihood with respect to tempering 
calcPF<-function(cycle,jobNum,llhdTemper,mcmcTemper,initResults, priorPar){
  par<-parMat[jobNum,]
  ################################################################################################
  names(par)<-parNames
  ################################################################################################
  lprior<-logPrior(par=par, priorPar=priorPar)
  llhd<-((initResults[[1]]-lprior)/mcmcTemper)*llhdTemper
  results<-initResults[[2]]
  return(list(llhd,results))
}

#############################################################################
# MCMC 
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
                              initResults,
                              parNames,
                              additional=FALSE){

  # Initialize Values
  n=length(obs)
  names(init)<-parNames
  boundMat<-priorPar[-1,]
  # Intialize containers
  alphaMat<-parMat<-candMat<-matrix(NA,nrow=iter,ncol=length(init))
  candMat[1,]<-parMat[1,]<-init ; alphaMat[1,]<-1
  colnames(parMat)<-colnames(candMat)<-colnames(alphaMat)<-parNames
  # Add posterior
  lpost<-vector("numeric");
  if(additional==FALSE){
    lpost[1]<-temper*(initResults[[1]])+logPrior(par=init, priorPar=priorPar) # Calculate Log Posterior from previous run  - IS step
  }else{
    lpost[1]<-initResults[[1]] # Use log posterior from previous run
  }
  # Note that we are doing this (MCMC Weight- Cumulative) x (Prev)/likelihoodweight + Log Prior, where Prev=(Log Likelihood *likelihoodweight)
  resultsList<-list();resultsList[[1]]<-initResults[[2]]
  curResultsList<-list();curResultsList[[1]]<-initResults[[2]]
  #Begin MCMC
  for(i in 2:iter){

   # Propose on the log scale from the multivariate truncated normal
    candMat[i,]<-c(parMat[i-1,1] , rtmvnorm(n = 1,mean = parMat[i-1,-1],
                                            sigma = propCov,
                                            lower = boundMat[,1],upper = boundMat[,2]) )

      
    curResults<-try({logPosterior_temper(par=as.numeric(candMat[i,]),
                                    priorPar = priorPar,
                                    obs = obs,
                                    inputDir = inputDir,
                                    outputDir = outputDir,
                                    j=jobNum , 
                                    temper=temper)},silent = TRUE)
    if (inherits(curResults, "try-error")){curResults<-list(-1e10,NA)}
        ############################################################################################################################
    # account for asymmetric proposal
    # Do some simulations to make sure this working properly
    adj1<-log(pmvnorm(upper= boundMat[,2],lower=boundMat[,1],mean=as.numeric(parMat[i-1,-1]),sigma=propCov))[1]#previous
    adj2<-log(pmvnorm(upper= boundMat[,2],lower=boundMat[,1],mean=as.numeric(candMat[i,-1]),sigma=propCov))[1]#current 
    alpha<-try(min(curResults[[1]]+adj1-lpost[i-1]-adj2,log(1)),silent = TRUE) # Asymmetric Proposal. Follow Example. Add and subtract constants. 
    if (inherits(alpha, "try-error")){alpha<--1e10}
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
    
    curResultsList[[i]]<-curResults[[2]] # Candidate Results -> Many are rejected in MH step
    
    # Update Gamma Parameter - Gibbs Update
    # lenOut<-length(resultsList[[i]][[1]])
    # alphaPar<-priorPar[1,1]+0.5*lenOut 
    # betaPar<-(sum((resultsList[[i]][[1]]-obs[1:lenOut])^2)+2*priorPar[1,2])/2
    parMat[i,1]<-candMat[i,1]<-parMat[i-1,1]
    
  }
  return(list(parMat,candMat,alphaMat,lpost,resultsList,curResultsList))
}

########################################################################################################
# Acceptance Rate
########################################################################################################
accRateFunc<-function(x){
  accRate<-(length(unique(x))-1)/(length(x)-1)
  return(accRate)
}
########################################################################################################
# Function to transform the log-likelihood to the tempered weights and ESS
########################################################################################################
fullToESS<-function(fullWeight,tempSeq,showESS=FALSE){
  weights<-fullWeight*tempSeq # Tempered Weights (from log likelihood)
  weights<-exp(weights-max(weights)) # Subtract max weight and exponentiate
  weights<-weights/sum(weights) # Normalize
  ESS<-1/sum((weights^2)) # ESS
  if(showESS==TRUE){
    return(ESS)
  }else{
    return(list(weights=weights,ESS=ESS))
  }
}
########################################################################################################
#Function that maximizes the ESS to a desired proportion of the ensemble size
########################################################################################################
maxESS<-function(fullWeight,tempSeq,ens ,prop){
  ESS<-fullToESS(fullWeight,tempSeq,showESS=TRUE)
  (ESS-prop*ens)^2
}
########################################################################################################
# Function that optimizes the weights according to some criteria
# Note that this relies on the maxESS() and fullToESS() functions
########################################################################################################
optimizeWeights<-function(weightVect,cumulTemp,ens,prop){
  if(1-cumulTemp<0.1){
    useTemp<-1-cumulTemp
  }else{
    optimTemp<-optimize(f = maxESS,fullWeight=weightVect,prop=prop,
                        ens=ens,interval = c(0.1,1-cumulTemp))
    useTemp<-optimTemp$minimum; 
  }
  
  
  useESSWeights<-fullToESS(fullWeight = weightVect,tempSeq = useTemp,showESS = FALSE)
  return(list(useTemp,useESSWeights))
}
########################################################################
# Function to calculate weights + resample
combineIS<-function(cycle,cumulTemp,prop=0.5,ens){
  setwd("/glade/scratch/sanjib/runA_rep/output/")
  fileDirLoad<-list.files(pattern = paste("PF_",cycle,sep=""))
  orderNum<-sapply(fileDirLoad,function(x) as.numeric(strsplit(x,split="_|[/.]")[[1]][3])) # Need to order it by numerical JobNum
  fileDirLoad<-fileDirLoad[order(orderNum)]
  
  ensFiles<-length(fileDirLoad) # Size of Samples from Importance Function
  weightVect<-vector("numeric")
  resultsList<-list()
  # Load Files
  for(i in 1:length(fileDirLoad)){
    load(fileDirLoad[i])
    resultsList[[i]]<-llhd_t[[2]]
    if(i==1){
      k<-length(jobPar)
      parWeightMat<-matrix(NA,nrow=ensFiles,ncol=k)
      parWeightMat[1,]<-jobPar
    }else{parWeightMat[i,]<-jobPar}
    weightVect[i]<-llhd_t[[1]] # Full Log-LIkelihood
  }
  #Optimize
  optimList<-optimizeWeights(weightVect=weightVect,cumulTemp=cumulTemp,ens=ensFiles,prop=prop)
  temperVal<-list();
  temperVal$cumulative<-cumulTemp+optimList[[1]]
  temperVal$incremental<-optimList[[1]]
  save(temperVal,file=paste("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/output_rep/temperVal_",cycle,".RData",sep=""))
  weights<-optimList[[2]]$weights
  # Weight 
  reSampleInd<-sample(x=1:ensFiles,size = ens,replace = TRUE,prob = weights) # Note that this allows for subsampling
  # Resample
  parMat<-parWeightMat[reSampleInd,]
  initResultsList<-list(weightVect[reSampleInd],resultsList[reSampleInd]) # Provides full log likelihood and results
  #Transformed COntribution Mat
  # Save final files
  save(parMat,parWeightMat,weightVect,weights,reSampleInd,
       initResultsList,optimList,
       file=paste("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/output_rep/rsParameters_",cycle,".RData",sep=""))
}
########################################################################
########################################################################
# Function to combine Mutation output 

combineMH<-function(cycle,ens,stage){
  setwd("/glade/scratch/sanjib/runA_rep/output/")
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
      load(paste("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/output_rep/totalParticles_",cycle,".RData",sep="")) 
      TotalParMat<-rbind(TotalParMat,amcmc.out[[1]])
    }
  
    acceptVect<-vector("numeric")
    acceptVect[i]<-accRateFunc(amcmc.out[[1]][,2])
    kMCMC<-length(amcmc.out[[5]])
    
  }else{
    parMat[i,]<-amcmc.out[[1]][nrow(amcmc.out[[1]]),]
    TotalParMat<-rbind(TotalParMat,amcmc.out[[1]])
    acceptVect[i]<-accRateFunc(amcmc.out[[1]][,2])
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
     file=paste("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/output_rep/totalParticles_",cycle,".RData",sep=""))

# Final Particles
save(parMat,acceptVect,initResultsList,
     file=paste("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/output_rep/mhParameters_",cycle,".RData",sep=""))
}
########################################################################################################################
# Combine Total Particles Function
combineTotalParticles<-function(cycle){

  if(cycle==1){
    load(paste("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/output_rep/totalParticles_",cycle,".RData",sep=""))
    masterTotalParticles<-TotalParMat
  }else{
    load("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/output_rep/totalParticles_1.RData")
    masterTotalParticles<-TotalParMat
    for(i in 2:cycle){# Loop Through all Total particles files
      load(paste("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/output_rep/totalParticles_",i,".RData",sep=""))
      masterTotalParticles<-rbind(masterTotalParticles,TotalParMat)  
    }
    
  }
  save(masterTotalParticles,file="/glade/u/home/sanjib/FamosHydroModel/Official_Fast/output_rep/masterTotalParticles.RData")
}

########################################################################
# Function to generate covariance matrix using all samples
# This uses the scaling factor from Rosenthal et. al (2008)
########################################################################
genPropMat<-function(cycle,scale){
  load(paste("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/output_rep/rsParameters_",cycle,".RData",sep=""))
  
  # Covariance Matrix
  if(cycle==1){
    load("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/output_rep/BeginCovMat_Tr.RData")
    CovMat<-CovMat*scale # Optimal proposal from Rosenthal et al. (2008)
  }else{
    load("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/output_rep/masterTotalParticles.RData")
    TotalParMat<-rbind(masterTotalParticles,parMat)
    uniqueID<-!duplicated(TotalParMat[,2])
    TotalParMat<-TotalParMat[uniqueID,]
    colnames(TotalParMat)<-parNames
    #####################################################################################################################
    CovMat<-cov(TotalParMat[,-1])*((2.38^2)/ncol(TotalParMat[,-1]))*scale # Optimal proposal from Rosenthal et al. (2008)
  }
   return(CovMat)
}


########################################################################
# Functions to compute Bhattacharryya Distance
########################################################################
# Count Table Function
########################################################################
countTable<-function(x,dat1,dat2){ # breaks, dataset 1 , dataset 2
  d1<-sum(dat1>=x[1]&dat1<x[2],na.rm = TRUE)/length(dat1)
  d2<-sum(dat2>=x[1]&dat2<x[2],na.rm = TRUE)/length(dat2)
  BC<-(sqrt(d1*d2)) # Bhattacharrya coefficient for break x
}
########################################################################
# Function to evaluate Bhattacharrya Distance 
########################################################################
battaDistance<-function(dat.current,dat.previous,breaks){
  penEnd<-quantile(dat.previous,probs = c(0.005,0.995),na.rm = TRUE) # Endpoints set at 0.005 and 0.995 quantiles
  partition<-c(-Inf,seq(penEnd[1],penEnd[2],length.out = breaks-2),Inf) # Sequence of cutoff points
  partitionTable<-rbind(partition[-length(partition)],partition[-1]) # Set cutoff points for each partition
  sampBD<--log(sum(apply(partitionTable,2,countTable,
                         dat1=dat.previous,
                         dat2=dat.current))) # Compute Empirical Bhattacharrya Distance
  return(sampBD)
}

