# # Copyright (C) 2022 Ben S. Lee
#email: slee287@gmu.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
####### # Copyright (C) 2022 Ben S. Lee
#email: slee287@gmu.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#######
# For study: ens = 2015; niter = 6 ; Proposal Matrix scaling 0.5
######## # Copyright (C) 2022 Ben S. Lee
#email: slee287@gmu.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################################################
# rWrapper_Continuous.R
# This helper file contains the R Wrapper to run the hydrological model. 
# This file consists of multiple nested functions. 
################################################################################################################

##################################################################
##################################################################
# Function 1: Input File Text
writeInput<-function(par , # Parameters
                     j ,  # JobNumber
                     dir, # Directory
                     outputDir #  OutputDirectory
                     ){
intervalMat<-c("20030601T00","20080331T23") # 2003/06/01-2008/03/31
  writeLabel<-c("ADD_PCTIM" ,"ADD_ADIMP" , "ADD_UZTWM" , "ADD_LZTWM" , 
                "ADD_LZFSM" , "ADD_LZFPM" , "ADD_LZSK" , "ADD_snow_SCF" , 
                "ADD_REXP" , "ADD_UZK" , "ADD_Q0CHN" , "ADD_QMCHN")

# Replace placeholders with values
readInputText  <- readLines("/glade/u/home/sanjib/FamosHydroModel/Official_Fast/precalibration/input/calsnow_template.card") # Read Template
readInputText  <- gsub(pattern = "ADD_output", replace = paste(outputDir,"/output",j,sep=""), # Replace output directory
                       x = readInputText)
  for(kk in 1:length(par)){ # Replace placeholders with parameter values
    readInputText  <- gsub(pattern = writeLabel[kk], replace = par[kk], x = readInputText) # Replace parameters
  }
# Replace Start and End Dates
readInputText  <- gsub(pattern = "ADD_DSTART", replace = intervalMat[1], # Replace Time Start
                       x = readInputText)
readInputText  <- gsub(pattern = "ADD_DEND", replace =intervalMat[2], # Replace Time End
                       x = readInputText)
# Write Files into the input directory
system(paste("rm -f ",paste(dir,"/input",j,".card",sep="") , sep="")) # Delete previous copies
writeLines(readInputText, con=paste(dir,"/input",j,".card",sep="")) # Write Calsnow file
}
##################################################################
##################################################################
# Function 2: Create New Output Directory
writeOutput<- function( j ,  # Job Number
                        dir # Directory
                       ){
newOutputDirectory<-paste(dir,"/output",j,sep="") # Output Directory
system(paste("rm -rf ",newOutputDirectory,";","mkdir ",newOutputDirectory , sep="")) # Remove Previous Versions
}

##################################################################
##################################################################
# Function 3: Run the Hydrological Model
runHydroModel<- function(j,dir){
newInputFile<-paste(dir,"/input",j,".card",sep="") # Read input directory
system(paste("/glade/u/home/sanjib/SBYN6/bin3/rdhm ",newInputFile)) # Run Model
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
  output<-readOutput( j = j , dir = outputDir) # Read the output and save
  ptFinal<-proc.time()-pt # Record the Wall Time
  output[[3]]<-ptFinal[3] # Save Walltime
  return(output) # Return output
}
  


####################################################################################################
## FOr Debugging use ONLY
####################################################################################################
# Function A: Run Model only
modelEval_only<-function( par, j , inputDir , outputDir){
  writeInput( par = par[-1] , j = j ,  dir = inputDir)  # Write Input
  writeOutput( j = j , dir = outputDir) # Write Output
  runHydroModel( j = j , dir = inputDir) # Run Model
}
# Function B: Read the output only
modelEval_read<-function( par, j , inputDir , outputDir){
  output<-readOutput( j = j , dir = outputDir)
  return(output) # Return output
}

# Function 4: Read the output file for a subset of dates
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

