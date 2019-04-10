#This is a R script on how to extract data from the 3_Runs Folder
#The 3_Runs folder holds numerous subfolders - one for each site
#In each site subfolder there should be two files - sw_input.RData and sw_output.RData
#The sw_input file is all the input files that were given to the model
#The sw_output file is all the daily raw data. Includes water balance components (i.e. transpiration, evaporation) and water output (volumetric water content and soil water potential)
#### BEWARE - there are empty or defenct 'slots'. Do not be alarmed
library(DBI)
library(RSQLite)
library(Rsoilwat31) # Install from GitHub : 
#                             system2(command = "git", args = "clone -b master --single-branch --recursive https://github.com/Burke-Lauenroth-Lab/Rsoilwat.git Rsoilwat")
#                             tools::Rcmd(args = paste("INSTALL Rsoilwat"))


# Steps for installing the old version of Rsoilwat needed to get this to work 
#  In a Mac open a command line shell 
#
# 1. clone rSOILWAT2 git repo into a directory 
#   >> git clone https://github.com/Burke-Lauenroth-Lab/Rsoilwat.git Rsoilwat_v31
# 2. Then checkout the the last commit of the repo before it changed to the 
#  the new package name.  
#   >> cd Rsoilwat_v31 
#   >> git checkout v1.2.4 
# 3. Then check the README.md file for the old instructions 
# 4. Then you need to initialize and update the submodules
#   >> git submodule init
#   >> git submodule update 
# 5. Then you need to run the installer shell script 
#   >> bash tools/TarAndInstall_31n.sh 
# 6. Then back up outside of the directory 
#   >> cd .. 
# 7. Then install the package from the command line 
#   >> R CMD INSTALL Rsoilwat_v31 
# 
# And that should work to get the package installed 


#Step 1 - Set working directory
dir.prj <- "data/SW_files/"

#Step 2- load and save files
#Input
#load(file.path('~/Desktop/Results2/1_SheepStation_SheepStation1', "sw_input.RData"),verbose=TRUE)
#swInput<-swRunScenariosData[[1]]
#swInput@soils

#Structure of Data
#Input
# slotNames(swInput)
# swInput@cloud

#Output
load(file.path(dir.prj,"sw_output.RData"),verbose=T)

#NOTE - the number (X) in the runData object (so runData[[X]]) refers to different scenario ids
#So for example extracting runData[[1]] grabs the current scenario data
#For the 'historical' runs and the 'current' runs there will only be the current scenario data
#For the Future climate runs there will be numerous scenarios  (37 to be precise))
#To know what scenario the number represents you need to refer back to the scenario_labels table in either the Weather sqlite database OR the Aggregated ouptut database
#Also please note that the 'years' in the future scenarios are broken. They reflect 'the current time period' but rest assured the data was generated with future climate data.

swOutputCurrent <- runData[[1]]

#swOutputFuture<-runData[[37]]#will only work if you have future data

#Further definitions of outputs are in the 'outsetup_v31.in' file
#slotNames(swOutputFuture)
# str(swOutputFuture@VWCMATRIC)
# str(swOutputFuture@VWCMATRIC@Day)
# 

VWC <-data.frame(swOutputCurrent@VWCMATRIC@Day)

VWC <- subset(VWC, !(Year == 2016 & DOY >= 270))

write.csv(VWC, 'data/daily_VWC.csv', row.names = FALSE)
