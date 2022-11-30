rm(list=ls())
setwd("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/famos/Official_Fast/")
load("output/PF_1_1.RData")
llhd_t[[1]]


obs<-read.table("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/SBYP1_obs.txt" , stringsAsFactors = FALSE)
obs<-obs[-(1:2),]
obsInd<-which(obs$V2%in%c(2010) & obs$V3%in%12)
obs<-obs[obsInd,]
obs<-as.numeric(obs[,5])*0.0283 # Need to convert

output<-llhd_t[[2]]
plot(x=1:31,y=obs ,typ="l", ylim=range(obs,output[,2]))
lines(x=1:31,y=output[,2], col="blue")

n=length(obs)
alpha_p<-0.001
beta_p<-0.001
alpha<-alpha_p+n/2+1
beta<-0.5*(sum((obs-output[,2])^2)+2*beta_p)

invX<-exp(seq(-20,10, length.out = 100))
y<-dinvgamma(x=invX , shape = alpha , rate = beta)
plot(x=invX, y=y, typ="l")



llhd<-vector("numeric")
for(k in 1:1000){
  llhd[k]<-sum(dnorm(x=obs, mean = output[,2] , sd= sqrt(k), log = TRUE)) 
}

plot(x=1:1000,y=llhd,typ="l")

dnorm

x=10
mu=1
sigma2=1
1/sqrt(2*pi*sigma2)*exp(-0.5*(x-mu)^2/sigma2)

dnorm(x=x,mean = mu, sd = sqrt(sigma2))



invX<-exp(seq(-20,2, length.out = 100))
y<-dinvgamma(x=invX , shape = 0.2 , rate = 0.2)
plot(x=invX, y=y, typ="l")
dev.off()



###############################################################################################
# Run Hydromodel
setwd("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/famos/Official_Fast/")
source("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/famos/Official_Fast/run/rWrapper.R")
load("output/mhParameters_0.RData")
parVal<-c(parMat[1,]) ;  j=1 ; 
inputDir<-"/gpfs/group/kzk10/default/private/hydrocalib/SGrove/famos/Official_Fast/input"
outputDir<-'/gpfs/group/kzk10/default/private/hydrocalib/SGrove/famos/Official_Fast/output'
pt<-proc.time()
output<-modelEval(par=parVal, j=j , inputDir=inputDir , outputDir=outputDir)
ptFinal<-proc.time()-pt
print(ptFinal)

save(output,ptFinal,file="output/testOutput.RData")