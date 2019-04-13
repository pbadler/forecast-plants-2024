############################################################# 
#
# NOTE: Running all the K-fold cross validation takes days 
#
#############################################################

rm(list = ls() ) 

# library(rstan)
# library(tidyverse)
# library(zoo)
# library(texreg)
# library(xtable)
# library(gridExtra)
# library(MASS)
# library(lsmeans)
# library(lme4)
# library(sheepweather)

# run data preparation files first --------------------------- # 

source('code/figure_scripts/save_plot_theme.R')

source('code/prepare_data/make_all_data.R')

# analysis pipeline ------------------------------------------ # 

# 1. Soil moisture analysis 
source('code/figure_scripts/plot_spring_soil_moisture_spot_measures.R')

# 2. Cover trends

source('code/figure_scripts/plot_cover.R')

# ----- Model Fitting/Selection ---------------------------------- # 

# 3. Fit candidate models and evaluate using K-Fold cross validation

source('code/analysis/fit_growth_models.R') # takes ~ a while 

source('code/analysis/fit_survival_models.R') # takes ~ a while 

source('code/analysis/fit_recruitment_models.R') # takes ~ a while 

# 4. Rank models based on LPPD 

source('code/analysis/rank_models.R')

# 5. Re-run top models for each species and vital rate 

source('code/analysis/fit_top_survival_models.R')
source('code/analysis/fit_top_growth_models.R')
source('code/analysis/fit_top_recruitment_models.R')

# ------ Simulation ---------------------------------------------  #

# 6. Generate IBM predictions based on top demographic models 

source('code/analysis/IBM_simulation.R')

# ----- Validation ----------------------------------------------- # 

# 7. Evaluate vital rate predictions on the held-out validation data (2012 to 2016)


# 8. Evalute cover predictions from IBMs 


# ----- Generate Figures ----------------------------------------- # 



# ----- Generate Tables  ----------------------------------------- # 


# ----- Knit manuscript  ------------------------------------------#

