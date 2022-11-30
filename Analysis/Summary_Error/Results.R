rm(list=ls())
##########################################################################################
##########################################################################################
##########################################################################################
setwd("~/Dropbox/FamosHydroModel/manuscript/revisionCode/")
##########################################################################################
##########################################################################################
##########################################################################################
load("resultsStremflow_validation_extreme.RData")
length_validation<-length(subsetFinalValidation)
load("../../Official_Fast/output_rep_var/mhParameters_6.RData")
variance_rep_var<-mean(parMat[,1])
load("../../Official_Fast/output_rep/rsParameters_4.RData")
variance_rep<-mean(parMat[,1])
par(mfrow=c(4,5), mar=c(2,2,2,2))
for(k in 1:length_validation){
  foo<-famosOutput_rep_var[,k]
  foo<-foo+rnorm(n=nrow(famosOutput_rep_var), mean=0, sd=sqrt(variance_rep_var))
  d1<-density(foo)
  
  bar<-famosOutput[,k]
  bar<-bar+rnorm(n=nrow(famosOutput), mean=0, sd=sqrt(variance_rep))
  d2<-density(bar)
  
  baz<-precalibrationOutput_Window[,k]
  baz<-baz+rnorm(n=nrow(precalibrationOutput_Window), mean=0, sd=sqrt(variance_rep))
  d3<-density(baz)
  
  plot(d1, main=dateVect[k] , ylim=range(d1$y,d2$y,d3$y) , xlim=range(d1$x,d2$x,d3$x))
  lines(d2, col="blue")
  lines(d3, col="red")
  abline(v=subsetFinalValidation[k])
}




rm(list=ls())
##########################################################################################
##########################################################################################
##########################################################################################
setwd("~/Dropbox/FamosHydroModel/manuscript/revisionCode/")
##########################################################################################
##########################################################################################
##########################################################################################
load("resultsStremflow_validation_extreme.RData")
length_validation<-length(subsetFinalValidation)
load("../../Official_Fast/output_rep_var/mhParameters_6.RData")
variance_rep_var<-0.00000001
load("../../Official_Fast/output_rep/rsParameters_4.RData")
variance_rep<-0.00000001
par(mfrow=c(4,5), mar=c(2,2,2,2))
for(k in 1:length_validation){
  foo<-famosOutput_rep_var[,k]
  foo<-foo+rnorm(n=nrow(famosOutput_rep_var), mean=0, sd=sqrt(variance_rep_var))
  d1<-density(foo, bw=350)
  
  bar<-famosOutput[,k]
  bar<-bar+rnorm(n=nrow(famosOutput), mean=0, sd=sqrt(variance_rep))
  d2<-density(bar, bw=350)
  
  baz<-precalibrationOutput_Window[,k]
  baz<-baz+rnorm(n=nrow(precalibrationOutput_Window), mean=0, sd=sqrt(variance_rep))
  d3<-density(baz, bw=350)
  
  plot(d1, main=dateVect[k] , ylim=range(d1$y,d2$y,d3$y) , xlim=range(d1$x,d2$x,d3$x))
  lines(d2, col="blue")
  lines(d3, col="red")
  abline(v=subsetFinalValidation[k])
}
