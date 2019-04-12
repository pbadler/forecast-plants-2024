######################################################################################
#
# Make STAN datalist  
#
#####################################################################################

# These retrieve and aggregate all the climate data 

source('code/prepare_data/make_rainfall.R')
source('code/prepare_data/aggregate_spot_VWC.R')
source('code/prepare_data/merge_decagon_with_climate_station_data.R')
source('code/prepare_data/soilMoistureTreatmentEffects.R')

source('code/prepare_data/aggregate_VWC_data.R')
source('code/prepare_data/make_climate_variables.R') 
source('code/prepare_data/prepare_climate_covariates.R')

# ----------------------------------------------------------------- 
source('code/prepare_data/get_all_demographic_data.R') # depends on access to driversdata 
source('code/prepare_data/calculate_cover_per_plot.R') # depends on access to driversdata 
source('code/prepare_data/prep_vital_rate_df.R')
