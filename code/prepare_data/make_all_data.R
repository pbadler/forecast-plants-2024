######################################################################################
#
# Make all data for analyses 
#
#####################################################################################
# devtools::install_gitub('akleinhesselink/sheepweather')
library(sheepweather)  # this package makes weather and soil moisture data available

# 1. Process the climate data  ------------------------- 
source('code/prepare_data/make_rainfall.R')
#source('code/prepare_data/aggregate_spot_VWC.R')
#source('code/prepare_data/merge_decagon_with_climate_station_data.R')
source('code/prepare_data/soilMoistureTreatmentEffects.R')
source('code/prepare_data/calc_treatment_effects_on_SOILWAT.R')

# source('code/prepare_data/aggregate_VWC_data.R')
# source('code/prepare_data/make_climate_variables.R') 
# source('code/prepare_data/prepare_climate_covariates.R')

# 2. Import and process the demographic data ------------
source('code/prepare_data/get_all_demographic_data.R') # depends on access to driversdata 
source('code/prepare_data/calculate_cover_per_plot.R') # depends on access to driversdata 
# source('code/prepare_data/prep_vital_rate_df.R')

# 3. clean up temp files

source('code/prepare_data/delete_temp_vital_rate_files.R')
