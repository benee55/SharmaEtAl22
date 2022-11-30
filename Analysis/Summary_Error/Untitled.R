rm(list=ls())
# Target: Normal Distribution with mean 10 and variance 20^2
# Proposal: Normal distribution with variable mean and variance. 
# Objective: Find the mean and variance of the proposal distribution that minimizes the KL divergence

KL_fun<-function(par){
  trueMean<-3
  trueSD<-5
  mean_q<-par[1]
  sd_q<-par[2]
  f<-function(x){
    dnorm(x,mean=trueMean , sd=trueSD, log=FALSE)*(dnorm(x,mean=trueMean , sd=trueSD, log=TRUE)-dnorm(x,mean=mean_q , sd=sd_q, log=TRUE))
    }
  foo<-integrate(f, lower = -200, upper = 200)

  return(foo$value)
  }

# KL_fun(par=c(20,2))


res_Optim<-optim(par = c(-3,1), fn =KL_fun,method="L-BFGS-B" , lower = c(-Inf,0.1), upper = c(Inf,Inf))
res_Optim


parMat<-expand.grid(seq(-5,5,length.out=50),seq(0.01,10,length.out=50))
vals<-apply(parMat, 1, KL_fun)
parMat[which.min(vals),]


