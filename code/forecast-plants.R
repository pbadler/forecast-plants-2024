############################################################# 
#
#
#
#############################################################

rm(list = ls() ) 
# library(sheepweather) # in-house package for weather data!!! 
# see https://github.com/akleinhesselink/sheepweather

# library(tidyverse)
# library(sheepweather)
# library(lubridate)
# library(zoo)
# library(texreg)
# library(xtable)
# library(gridExtra)
# library(MASS)
# library(emmeans)
# library(lme4)
# library(egg)
# library( optimx )

# run data preparation files first --------------------------- # 
source('code/figure_scripts/save_plot_theme.R')

source('code/prepare_data/make_all_data.R') # takes a minute or two 

# Some treatment figures  ------------------------------------------ # 

# 1a. Treatment effects plot --spot measurements 
source('code/figure_scripts/plot_spring_soil_moisture_spot_measures.R')

# 1b. Treatment effects plot -- daily soil raw average vs modeled average 
source('code/figure_scripts/plot_soil_moisture_effects.R')

# 2. Daily SoilWat vs. average ambient soil moisture 
source('code/figure_scripts/plot_soilwat_to_ambient_comparison.R')

# 3. Shelter temperature plot 
source('code/figure_scripts/test_for_shelter_temperature_effects.R')

# 4. Soil moisture aggregated by season 
source('code/figure_scripts/plot_seasonal_soil_moisture.R')

# 5. Temperature aggregated by season  
source('code/figure_scripts/plot_longterm_TAVG.R')

# 2. Cover trends

source('code/figure_scripts/plot_cover.R')

# ----- Model Fitting/Selection ---------------------------------- # 
# Use sliding windows from ClimWin package find best climate windows 
# for growth and survival of each species. 
source('code/analysis/prep_daily_weather.R')

## Warning ## takes a while to run for all 4 species !!! 
source('code/analysis/growth_daily_climWin.R')

## WARNING ## takes hours to run for all 4 species (especially PSSP) !!!! 
source('code/analysis/survival_daily_climWin.R') # !!!!!!!

# plot  windows and save best windows to tables 
source('code/analysis/plot_window_comparisons.R')

# Refit models with chosen climate windows and save to output folder 
source('code/analysis/train_growth_models.R')
source('code/analysis/train_survival_models.R')

# Validate model performance 
source('code/analysis/growth_validation.R')
source('code/analysis/survival_validation.R')


# ------ Simulation ---------------------------------------------  #

# 6. Generate IBM predictions based on top demographic models 

source('code/analysis/IBM_simulation.R')
source('code/figure_scripts/plot_IBM_predictions.R') # a work in progress

# ----- Validation ----------------------------------------------- # 

# 7. Evaluate vital rate predictions on the held-out validation data (2012 to 2016)

source('code/analysis/growth_validation.R')
source('code/analysis/survival_validation.R')

# 8. Generate cover predictions from IBMs 


# ----- Generate Figures ----------------------------------------- # 



# ----- Generate Tables  ----------------------------------------- # 


# ----- Knit manuscript  ------------------------------------------#

