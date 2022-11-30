# Generate Latin Hypercube sample for model parameters
setwd("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/LHSoutput/input/")
library(lhs)
set.seed(2020)
p<-12 # parameters
ens<-1000 # samples

A <- randomLHS(ens, p)
# A1 <- optimumLHS(ens, p, maxSweeps = 4, eps = 0.01)
A <- maximinLHS(ens, p)

# Sample frm the Priors
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

rangeMat<-rbind(t(boundMat),A)

sampRange<-function(x){
  qunif(x[-c(1,2)],min = x[1],max=x[2])
}

# Latin Hypercube Samples
parMat<-apply(rangeMat,2,sampRange)
colnames(parMat)<-c("PCTIM" , "ADIMP" , "UZTWM" ,"LZTWM" , "LZFSM" , "LZFPM" , "LZSK" , "UZFWM" , "REXP" , "UZK" , "Q0CHN" , "QMCHN")

save(parMat,file="/gpfs/group/kzk10/default/private/hydrocalib/SGrove/LHSoutput/input/LHS_Inputs.RData")
