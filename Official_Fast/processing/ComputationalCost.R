multiplierParticle<-2*15 # cycles * mutations
multiplierNaive<-10*100
start<-10/60/60
end<-15*16/60+1
X<-seq(start,end,length.out = 1000)
particleY<-multiplierParticle*X
iter<-100000
MCMCY<-X*iter
naivePF<-X*multiplierNaive

dev.off()
pdf(file="ComputationalCost.pdf", height=8.5, width=11)
par(mar=c(4.5,6,2,2))
plot(x=X,y=particleY,typ="l",ylim=range(MCMCY+1e+7,particleY),xlim=range(-0.1,X),
     xlab="Computer Model Run Time (hours)",ylab="",lwd=2,
     log=c("y"),yaxt='n', main="Calibration Walltime vs. Computer Model Run Time")

timeSeq<-c(1,6,24,168,30*24,365*24,365*24*5,365*24*10,365*24*50)
timeLab<-c("1 hour", "6 hours", "1 day", "1 week", "1 month",
           "1 year", "5 years", "10 years", "50 years")
text(y=timeSeq, x=par("usr")[1] - 0.3, labels = timeLab,
     srt = 0, xpd = TRUE,cex=1.3)
text(y=365*24*500, x=par("usr")[1] - 0.3, labels= "Calibration \n Walltime",
     srt = 0, xpd = TRUE,cex=1.5)

text(y=1,x=0.65,labels = "PSU3D-ICE \n 80km",col="black",cex=1.3)
text(y=1,x=1+0.4,labels = "PSU3D-ICE \n 40km",col="blue",cex=1.3)
text(y=1,x=4.4,labels = "PSU3D-ICE \n 20km",col="red",cex=1.3)
axis(2, at=timeSeq, labels=NA,tick = TRUE)
lines(x=X,y=naivePF,lty=2,lwd=2)
lines(x=X,y=MCMCY,lty=3,lwd=2)

abline(v=15/60,lwd=2)
abline(v=15*4/60,lwd=2,col="blue")
abline(v=15*16/60,lwd=2,col="red")
legend("topright",  legend=c("MCMC","Naive","Adaptive"),
       lty=c(3,2,1),cex=1.5,lwd=c(2,2,2))

dev.off()

plot(x=X,y=particleY,lty=2,ylim=range(particleY,naivePF),typ="l")
lines(x=X,y=naivePF,lty=2,col="red")


abline(h=foo)
# amount of computing time available
standardCheyenne<-1000000*60
nParticles<-2000

foo<-standardCheyenne/nParticles
abline(h=foo,col="red")

# 2km resolution
standardCheyenne<-400000
nParticles<-seq(100,100000,length.out = 100)
plot(x=nParticles,y=standardCheyenne/nParticles,typ="l",
     log=c("y","x"))
abline(h=15/60*16,lwd=2)
