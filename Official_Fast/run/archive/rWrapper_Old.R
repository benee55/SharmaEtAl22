# RWrapper for Hydro Model

library(invgamma)

##################################################################
# Begin R Wrapper Function
##################################################################
##################################################################
# Function 1: Input File Text
writeInput<-function(par , j , dir){
writeLabel<-c("ADD_PCTIM" ,"ADD_ADIMP" , "ADD_UZTWM" , "ADD_LZTWM" , 
                "ADD_LZFSM" , "ADD_LZFPM" , "ADD_LZSK" , "ADD_UZFWM" , 
                "ADD_REXP" , "ADD_UZK" , "ADD_Q0CHN" , "ADD_QMCHN")
readInputText  <- readLines(paste(dir,"/calsnow_Ben_template.card",sep="")) # Read Template
readInputText  <- gsub(pattern = "ADD_output", replace = paste("output",j,sep=""), # Replace output directory
                       x = readInputText)
  for(kk in 1:length(par)){
    readInputText  <- gsub(pattern = writeLabel[kk], replace = par[kk], x = readInputText) # Replace parameters
  }
  writeLines(readInputText, con=paste(dir,"/input",j,".card",sep=""))
}

# Function 2: Delete Directory and Create
writeOutput<- function(j,dir){
newOutputDirectory<-paste(dir,"/output",j,sep="")
system(paste("rm -rf ",newOutputDirectory,";","mkdir ",newOutputDirectory , sep=""))
}

# Function 3: Run Model + Save FIle
runHydroModel<- function(j,dir){
newInputFile<-paste(dir,'/input',j,'.card',sep="")
system(paste("/gpfs/group/kzk10/default/private/hydrocalib/SGrove/bin/rdhm ",newInputFile))
}

# Function 4: Read OutputFile
readOutput<-function(j,dir){
  output<-system(paste("cat ",dir,"/output",j,"/SBYP1_discharge_outlet.ts",sep=""), intern=TRUE)
  cbind(as.numeric(unlist(lapply(output[135:165],substr,start = 11 , stop = 16))),
  as.numeric(unlist(lapply(output[135:165],substr,start = 22 , stop = 35))))
  # cbind(as.numeric(unlist(lapply(output[166:530],substr,start = 11 , stop = 16))),
                     # as.numeric(unlist(lapply(output[166:530],substr,start = 22 , stop = 35))))
  
}

# Combine Run File + Read
runReadModel<-function( j , inputDir , outputDir){
  runHydroModel( j = j , dir = inputDir)
  readOutput( j = j , dir = outputDir)
}
  


# Wrapper for R -> Runs Hydromodel and outputs log posterior + model output


