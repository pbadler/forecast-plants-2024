# Forecast-plants 
Data and analyses for an Adler lab project to predict plant population response to experimental precipitation regimes. 

# Setup 

## open forecast-plants.Rproj file 
This is an Rproject for use in Rstudio. Open the .Rproj file found in the home directory with Rstudio. 
This will automatically set the appropriate working directory. 

## update submodules 
This project requires the use of data from the following submodule 
- driversdata

This submodules can be initialized and updated by running 
git submodule init 
and then 
git submodule update 

The driversdata submodule loads all the relevant historical and experimental demographic 
data for use in the analyses.  This is currently a private repository. 

## Run the analysis 
Run the code/analysis_pipeline.R script to process data and run all analyses

This takes a long time to loop through all the species, cross validation folds and model variants. 

