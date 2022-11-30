rm(list=ls())
setwd("~/Dropbox/FamosHydroModel/Analysis/")
load("calibrationParameters.RData")
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

# Parameter Names
parNames<-c("PCTIM" , "ADIMP" , "UZTWM" ,"LZTWM" , 
            "LZFSM" , "LZFPM" , "LZSK" , "snow_SCF" ,
            "REXP" , "UZK" , "Q0CHN" , "QMCHN")

# png(file = "Parameters.png", height = 1000, width=1200)
par(mfrow=c(4,4), mar=c(2,2,2,2))
plot(x=1,y=1,typ="n")
legend("center" , legend=c("Famous - 12 par","Famous - 4 par" ,
                           "Precalibration - 12 par" , "Handtune - 12 par" , "Prior"),
       lty=c(1,1,1,1,1), col=c("black","red","gray","green","blue"),
       lwd=c(2,2,2,2,2) , cex=1.3 , bty = "n")
for(k in 1:ncol(famosParMat)){
  d1<-density(famosParMat[,k])
  d4<-density(precalibrationParMat[,k])
  if(k %in% c(1,2,11,12)){
    hk<-ifelse(k%in%c(1,2), k, k-8)
    d2<-density(famos4ParMat[,hk])
    plot(d1, main=parNames[k] ,  lwd=2,
         xlim=range(d1$x, d2$x, d4$x , handTuneParMat[k],boundMat[k,]), 
         ylim=range(d1$y, d2$y, d4$y))  
    lines(d2, col="red", lwd=2)
  }else{
    plot(d1, main=parNames[k] , lwd=2,
         xlim=range(d1$x, d4$x , handTuneParMat[k],boundMat[k,]), 
         ylim=range(d1$y, d4$y))
    
  }
  lines(d4, col="gray", lwd=2)
  abline(v=handTuneParMat[k] , col="green", lwd=2)
  abline(v=boundMat[k,] , col="blue", lwd=2)
  abline(v=handTuneParMat[k] , col="green", lwd=2)
}

# dev.off()
