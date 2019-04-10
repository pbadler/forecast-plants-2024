######################################################################################
#
# Make STAN datalist  
#
#####################################################################################

# These retrieve and aggregate all the climate data 

# This depends on soilwat data sent by Caitlin Andrews -------------- # 

# Depends Rsoilwat31 
#source('code/prepare_data/climate/ExtractData_3Runs.R') # don't run

# don't run this unless you can install Rsoilwat31, an old version of
# the soilwat package. 
# Instructions for installing the version are given in the script.  
# ------------------------------------------------------------------ #

source('code/prepare_data/climate/aggregate_spot_VWC.R')
source('code/prepare_data/climate/soilMoistureTreatmentEffects.R')
source('code/prepare_data/climate/aggregate_VWC_data.R')
source('code/prepare_data/climate/make_climate_variables.R') 
source('code/prepare_data/climate/prepare_climate_covariates.R')

# ----------------------------------------------------------------- 
source('code/prepare_data/get_all_demographic_data.R') # depends on access to driversdata 
source('code/prepare_data/calculate_cover_per_plot.R') # depends on access to driversdata 
source('code/prepare_data/prep_vital_rate_df.R')
