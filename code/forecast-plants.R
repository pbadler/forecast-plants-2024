############################################################# 
#
# NOTE: Running all the K-fold cross validation takes days 
#
#############################################################

rm(list = ls() ) 

# library(sheepweather) # in-house package for weather data!!! 
# see https://github.com/akleinhesselink/sheepweather

# library(rstan)
# library(tidyverse)
# library(zoo)
# library(texreg)
# library(xtable)
# library(gridExtra)
# library(MASS)
# library(lsmeans)
# library(lme4)

# run data preparation files first --------------------------- # 

source('code/figure_scripts/save_plot_theme.R')

source('code/prepare_data/make_all_data.R') # takes a minute or two 

# Some treatment figures  ------------------------------------------ # 

# 1. Treatment effects plot 
source('code/figure_scripts/plot_spring_soil_moisture_spot_measures.R')

# 2. Shelter temperature plot 
source('code/figure_scripts/test_for_shelter_temperature_effects.R')

# 3. Soil moisture aggregated by season 
source('code/figure_scripts/plot_seasonal_soil_moisture.R')

# 4. Temperature aggregated by season  
source('code/figure_scripts/plot_longterm_TAVG.R')

# 2. Cover trends

source('code/figure_scripts/plot_cover.R')

# ----- Model Fitting/Selection ---------------------------------- # 

# 3. Fit candidate models and evaluate using K-Fold cross validation

source('code/analysis/fit_growth_models.R') # takes a while 

source('code/analysis/fit_survival_models.R') # takes a while 

source('code/analysis/fit_recruitment_models.R') # takes a while 

# 4. Rank models based on LPPD 

source('code/analysis/rank_models.R')

# 5. Re-run top models for each species and vital rate 

source('code/analysis/fit_top_survival_models.R')
source('code/analysis/fit_top_growth_models.R')
source('code/analysis/fit_top_recruitment_models.R')

# ------ Simulation ---------------------------------------------  #

# 6. Generate IBM predictions based on top demographic models 

source('code/analysis/IBM_simulation.R')
source('code/figure_scripts/plot_IBM_predictions.R') # a work in progress

# ----- Validation ----------------------------------------------- # 

# 7. Evaluate vital rate predictions on the held-out validation data (2012 to 2016)

source('code/analysis/model_validation.R')

# 8. Evalute cover predictions from IBMs 


# ----- Generate Figures ----------------------------------------- # 



# ----- Generate Tables  ----------------------------------------- # 


# ----- Knit manuscript  ------------------------------------------#

