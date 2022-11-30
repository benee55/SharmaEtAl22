source("~/Box Sync/LeeEtal-PSUICE-3D/calibration/Official/run/mcmc_source.R")

windows<-matrix(c(lower.window,upper.window),nrow=2,ncol=5,byrow = TRUE)

par(mfrow=c(2,3))
x<-seq(windows[1,1]-1,windows[2,1]+1,length.out = 10000)
plot(x=x,y=dtnorm(x,mean=windows[1,1]+abs(windows[2,1]-windows[1,1])/2,
       sd=30,lower=windows[1,1],upper=windows[2,1]),typ="l",main="Pliocene")

x<-seq(windows[1,2]-1,windows[2,2]+1,length.out = 10000)
plot(x=x,y=dtnorm(x,mean=windows[1,2]+abs(windows[2,2]-windows[1,2])/2,
       sd=10,lower=windows[1,2],upper=windows[2,2]),typ="l",main="LIG")

x<-seq(windows[1,3]-1,windows[2,3]+1,length.out = 10000)
plot(x=x,y=dtnorm(x,mean=windows[1,3]+abs(windows[2,3]-windows[1,3])/2,
       sd=20,lower=windows[1,3],upper=windows[2,3]),typ="l",main="LGM")

x<-seq(windows[1,4]-1,windows[2,4]+1,length.out = 10000)
plot(x=x,y=dtnorm(x,mean=windows[1,4]+abs(windows[2,4]-windows[1,4])/2,
       sd=abs((windows[2,4]-windows[1,4])/4),
       lower=windows[1,4],upper=windows[2,4]),typ="l",main="Modern Vol") # assume window contains 95%% so STD is window/4

x<-seq(windows[1,5]-1,windows[2,5]+1,length.out = 10000)
plot(x=x,y=dtnorm(x,mean=windows[1,5]+abs(windows[2,5]-windows[1,5])/2,
       sd=abs((windows[2,5]-windows[1,5])/4),
       lower=windows[1,5],upper=windows[2,5]),typ="l",main="Modern Area")

