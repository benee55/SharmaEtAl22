# RWrapper for Hydro Model


##################################################################
# Begin R Wrapper Function
##################################################################
##################################################################
# Function 1: Input File Text

writeInput<-function(par , # Parameters
                     j ,  # JobNumber
                     dir, # Directory
                     outputDir #  OutputDirectory
                     ){
intervalMat<-c("20030601T00","20080331T23") # 2003/06/01-2008/03/31
# intervalMat<-c("20040901T00","20040930T23") # 2004/09/01-2004/09/30
  writeLabel<-c("ADD_PCTIM" ,"ADD_ADIMP" , "ADD_UZTWM" , "ADD_LZTWM" , 
                "ADD_LZFSM" , "ADD_LZFPM" , "ADD_LZSK" , "ADD_snow_SCF" , 
                "ADD_REXP" , "ADD_UZK" , "ADD_Q0CHN" , "ADD_QMCHN")

# Replace placeholders with values
readInputText  <- readLines("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/precalibration/input/calsnow_template.card") # Read Template
readInputText  <- gsub(pattern = "ADD_output", replace = paste(outputDir,"/output",j,sep=""), # Replace output directory
                       x = readInputText)

  for(kk in 1:length(par)){
    readInputText  <- gsub(pattern = writeLabel[kk], replace = par[kk], x = readInputText) # Replace parameters
  }

readInputText  <- gsub(pattern = "ADD_DSTART", replace = intervalMat[1], # Replace Time Start
                       x = readInputText)
readInputText  <- gsub(pattern = "ADD_DEND", replace =intervalMat[2], # Replace Time End
                       x = readInputText)

system(paste("rm -f ",paste(dir,"/input",j,".card",sep="") , sep=""))
writeLines(readInputText, con=paste(dir,"/input",j,".card",sep="")) # Write Calsnow file
}

# Function 2: Delete Directory and Create
writeOutput<- function( j ,  # Job Number
                        dir # Directory
                       ){
newOutputDirectory<-paste(dir,"/output",j,sep="")
system(paste("rm -rf ",newOutputDirectory,";","mkdir ",newOutputDirectory , sep=""))
}

# Function 3: Run Model + Save FIle
runHydroModel<- function(j,dir){
newInputFile<-paste(dir,"/input",j,".card",sep="")
system(paste("/glade/u/home/sanjib/SBYN6/bin3/rdhm ",newInputFile))
}


# Function 4: Read OutputFile
readOutput<-function(j,dir){
  # Dates evaluated within each interval
  output<-system(paste("cat ",dir,"/output",j,"/SBYP1_discharge_outlet.ts",sep=""), intern=TRUE)
  endInd<-length(output) # Row of final date
  startInd<-89 #Row of first date
  flow<-as.numeric(unlist(lapply(output[startInd:endInd],substr,start = 22 , stop = 35))) #Streamflow for specified dates
  dates<-unlist(lapply(output[startInd:endInd],substr,start = 11 , stop = 16)) #Streamflow for specified dates
  useDate<-c(c("190904","200904"),
             c("150105","160105", "300305" , "310305" , "030405" , "040405" ,"050405" , "060405"),
             c("011205"),
             c("280606" , "290606" , "300606"),
             c("181106"),
             c("160307" , "170307"),
             c("080208" , "060308" , "090308" , "100308"))
  extremeFlow<-flow[which(dates%in%useDate)]
  return(list(extremeFlow,flow))
}



# Combine Run File + Read
modelEval<-function( par, j , inputDir , outputDir){
  pt<-proc.time()
  writeInput( par = par[-1] , j = j , dir = inputDir , outputDir=outputDir)  # Write Input
  writeOutput( j = j ,dir = outputDir) # Write Output
  runHydroModel( j = j , dir = inputDir) # Run Model
  output<-readOutput( j = j , dir = outputDir)
  ptFinal<-proc.time()-pt
  output[[3]]<-ptFinal[3]
  return(output) # Return output
}
  


####################################################################################################
## FOr Debugging
####################################################################################################
# Combine Run File + Read
modelEval_only<-function( par, j , inputDir , outputDir){
  writeInput( par = par[-1] , j = j ,  dir = inputDir)  # Write Input
  writeOutput( j = j , dir = outputDir) # Write Output
  runHydroModel( j = j , dir = inputDir) # Run Model
}

modelEval_read<-function( par, j , inputDir , outputDir){
  output<-readOutput( j = j , dir = outputDir)
  return(output) # Return output
}

 
# Function 4: Read OutputFile for subset
readOutput_subset<-function(j,dir){
  # Dates evaluated within each interval
  output<-system(paste("cat ",dir,"/output",j,"/SBYP1_discharge_outlet.ts",sep=""), intern=TRUE)
  endInd<-length(output) # Row of final date
  startInd<-89 #Row of first date
  flow<-as.numeric(unlist(lapply(output[startInd:endInd],substr,start = 22 , stop = 35))) #Streamflow for specified dates
  dates<-unlist(lapply(output[startInd:endInd],substr,start = 11 , stop = 16)) #Streamflow for specified dates
  useDate<-c(c("190904","200904"),
             c("150105","160105", "300305" , "310305" , "030405" , "040405" ,"050405" , "060405"),
             c("011205"),
             c("280606" , "290606" , "300606"),
             c("181106"),
             c("160307" , "170307"),
             c("080208" , "060308" , "090308" , "100308"))
  flow[which(dates%in%useDate)]
}

