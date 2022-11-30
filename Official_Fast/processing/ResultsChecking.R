####################################################################################
setwd("../output/")

####################################################################################
#Combined Results
#####################################################################################
# Select Random Samples
#sampID<-sample(1:71,5)
CheckEntireData<-function(parInd){
ens<-length(list.files(pattern="PF_1_"))
totCycles<-length(list.files(pattern="total"))
sampID<-1:ens
mcmcCycles<-vector("numeric")
for(j in 1:totCycles){
  mcmcCycles[j]<-length(list.files(pattern=paste("MCMC_",j,"_",sep="")))/ens
  
}

MasterList<-list()
for(kk in 1:length(mcmcCycles)){
  
# Load Initial Parameters for Each Cycle
if(kk==1){
  load("mhParameters_0.RData")
}else{
    load(paste("mhParameters_",kk-1,".RData",sep=""))
  }
InitialPar<-parMat[sampID,parInd]
# Load Parameters from the Importance Sampling Step
vectPF1<-vector("numeric")
for(i in 1:length(sampID)){
  h<-sampID[i]
  load(paste("PF_",kk,"_",h,".RData",sep=""))
  vectPF1[i]<-(jobPar[parInd])
}
#Match Resampled Particle to Original Sample
load(paste("rsParameters_",kk,".RData",sep=""))
MasterVect<-cbind(InitialPar,vectPF1,parWeightMat[sampID,parInd],parMat[sampID,parInd])
#Match 
vectMCMC1<-matrix(NA,nrow=length(sampID),ncol=mcmcCycles[kk]*2)

for(i in 1:length(sampID)){
  for(j in 1:mcmcCycles[kk]){
  h<-sampID[i]
  load(paste("MCMC_",kk,"_",j,"_",h,".RData",sep=""))
  if(j==1){
    IntVect<-amcmc.out[[1]][c(1,nrow(amcmc.out[[1]])),parInd]
  }else{
    IntVect<-c(IntVect,amcmc.out[[1]][c(1,nrow(amcmc.out[[1]])),parInd])
  }
  }
  vectMCMC1[i,]<-IntVect
}

MasterVect<-cbind(MasterVect,vectMCMC1)
if(mcmcCycles[kk]==3){
  colnames(MasterVect)<-c("Initial","PF","RS-A","RS-B","MCMC1-A","MCMC1-Z","MCMC2-A","MCMC2-Z","MCMC3-A","MCMC3-Z")  
}else{
  colnames(MasterVect)<-c("Initial","PF","RS-A","RS-B","MCMC1-A","MCMC1-Z","MCMC2-A","MCMC2-Z")  
}

MasterList[[kk]]<-MasterVect
}

##########################################################################################################################################################
# Check If Certain Samples are Staying the same 
#(i.e. Resampled Particle is used as intial state for MCMC + Last MCMC sample is used as the Particle for IS)
check3<-rbind(c(1,2),c(1,3),c(2,3),c(4,5),c(6,7),c(8,9))
check2<-rbind(c(1,2),c(1,3),c(2,3),c(4,5),c(6,7))



checkEquality<-function(x,kk){MasterList[[kk]][,x[1]]==MasterList[[kk]][,x[2]]}
CheckResults<-list();CheckResults[[1]]<-vector("logical");CheckResults[[2]]<-vector("logical")
for(kk in 1:length(MasterList)){
  if(ncol(MasterList[[kk]])==13){
    CheckResults[[1]][kk]<-all(apply(apply(check3,1,checkEquality,kk=kk),2,sum)==length(sampID))
  }else{
    CheckResults[[1]][kk]<-all(apply(apply(check2,1,checkEquality,kk=kk),2,sum)==length(sampID))
  }
  
  if(kk!=1){
    checkInter<-c(ncol(MasterList[[kk-1]]),1)
    CheckResults[[2]][kk]<-all(MasterList[[kk-1]][,checkInter[1]]==MasterList[[kk]][,checkInter[2]])
  }
}
print("Check Flow of Particles")
print(CheckResults)
#SUMMARY: EVERYTHING CHECKS OUT

##########################################################################################################################################################
# Check if Total Particles Match 
totalList<-list.files(pattern = "total")
orderNum<-sapply(totalList,function(x) as.numeric(strsplit(x,split="_|[/.]")[[1]][2]))
totalList<-totalList[sort(orderNum)]
iter=45
TotalParticlesMaster<-list()
for(kk in 1:length(totalList)){
load(totalList[kk])  
TotalParticlesMaster[[kk]]<-matrix(NA,nrow=ens,ncol=2*nrow(TotalParMat)/(iter*ens))

for(j in 1:ens){
  if(ncol(TotalParticlesMaster[[kk]])==4){
    orderJ<-c(iter*(j-1)+1 , iter*(j-1)+iter ,
              (iter*ens)+iter*(j-1)+1 , (iter*ens)+iter*(j-1)+iter) 
  }else{
    orderJ<-c(iter*(j-1)+1 , iter*(j-1)+iter ,
              (iter*ens)+iter*(j-1)+1 , (iter*ens)+iter*(j-1)+iter,
              (2*iter*ens)+iter*(j-1)+1 , (2*iter*ens)+iter*(j-1)+iter )  
  }
  
  TotalParticlesMaster[[kk]][j,]<-TotalParMat[orderJ,parInd]
}
TotalParticlesMaster[[kk]]<-TotalParticlesMaster[[kk]][sampID,]
}

#### Check if this is the same as the MasterList
totalMasterCheck<-vector("logical")
for(kk in 1:length(TotalParticlesMaster)){
  totalMasterCheck[kk]<-all(TotalParticlesMaster[[kk]]==MasterList[[kk]][,5:ncol(MasterList[[kk]])])
}
print(totalMasterCheck)

#SUMMARY: EVERYTHING CHECKS OUT
################################################################################################################

##########################################################################################################################################################
# Check if Importance Sampling is properly being done (i.e. higher weights -> higher chance of being resampled)
rsList<-list.files(pattern = "rsParameters_")
orderNum<-sapply(rsList,function(x) as.numeric(strsplit(x,split="_|[/.]")[[1]][2]))
rsList<-rsList[sort(orderNum)]
isList<-list()
for(kk in 1:length(rsList)){
  load(rsList[kk])
  isList[[kk]]<-cbind(parMat[,parInd],parWeightMat[,parInd],reSampleInd,weights,weightVect)
}

checkIS<-sapply(isList,function(x) c(length(unique(x[,3])),sum(x[,4])))
rownames(checkIS)<-c("UniqueSamp","WeightSum")
print(checkIS)

return(list(MasterList,TotalParticlesMaster,isList))
}


################################################################################################################
#Checks across all parameters
CheckResultsList<-apply(matrix(1:11,nrow=1),2,CheckEntireData)
save(CheckResultsList,file="../processing/Checkresults.RData")
