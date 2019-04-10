# Precip 
Repository for field experiment seeking to predict plant community response to experimental precipitation regimes. 

# Setup 

## open .Rproj file 
This is an Rproject for use in Rstudio. Open the .Rproj file found in the home directory with Rstudio. 
This will automatically set the appropriate working directory. 

## update submodules 
This project requires the use of data and functions from two submodules 
- Rsoilwat31 
- driversdata

These submodules can be initialized and updated by running 
git submodule init 
and then 
git submodule update 

Once the Rsoilwat31 submodule is updated follow instructions to install the Rsoilwat31
package on your system. 

The driversdata submodule loads all the relevant historical and experimental demographic 
data for use in the analyses.  This is currently a private repository. 

## Run the analysis 
Run the code/analysis_pipeline.R script to process data and run all analyses

This takes a long time to loop through all the species, cross validation folds and model variants. 

