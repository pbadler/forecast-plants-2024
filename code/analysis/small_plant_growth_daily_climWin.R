rm(list = ls())
library(tidyverse)
library(lubridate)
library(climwin)
library(lme4)
library(optimx)
library(dfoptim)

# LMER optimization options
control_lmer = lmerControl(
  optimizer = "optimx",
  calc.derivs = FALSE,
  optCtrl = list(
    method = "nlminb",
    starttests = FALSE,
    kkt = FALSE
  )
)
control_lmer$optCtrl$eval.max <- 1e8
control_lmer$optCtrl$iter.max <- 1e8


source('code/analysis/functions.R')

# Variables -------------------------------- : 
last_year <- 2010 # last year of training data, everything earlier is used 
sp_list <- c('ARTR') #, 'HECO', 'POSE', 'PSSP')

# ClimWin Window Settings Monthly
window_open_max <- 24
window_open_min <- 1
window_exclude_dur <- 1
window_exclude_max <- 5

size_cutoff <- -1 # do small plants
species <- 'ARTR' # for testing 

# Climate and VWC data  ------------------- # 
quad_info <- read_csv( file = 'data/quad_info.csv')
daily_weather <- 
  read_csv('data/temp/daily_weather_for_models.csv') %>% 
  filter( Treatment == 'Control') %>% 
  filter( year <= last_year)

## ------------------------------------------------- 

for(species in sp_list){ 
  # loop Species 
  # 1. Find best univariate climate model 
  # 2. Redo ClimWin model selection with additional variable
  #     a. When best climate variable is temperature add VWC window 
  #     b. When best climate variable is VWC add temperature window 
  # 3. Save top model and data 
  
  growth <- prep_growth_for_climWin(species, 
                                    last_year = last_year, 
                                    quad_info = quad_info, 
                                    size_cutoff = -Inf)
  
  growth <- growth %>% filter( area0 <= size_cutoff) # Small plants 
  
  ##baseline model for growth/size
  m_baseline <- lmer( area ~ W.intra + (1|year),
                      data = growth,
                      REML = F,
                      control = control_lmer)
  model_type <- "mer"
  
  write_rds( m_baseline, paste0( 'output/growth_models/', species, '_growth_small_', model_type, '_baseline.rds'))
  
  growthWin <- slidingwin(xvar = list(TMAX_scaled = daily_weather$TMAX_scaled, 
                                      TAVG_scaled = daily_weather$TAVG_scaled,
                                      TMIN_scaled = daily_weather$TMIN_scaled, 
                                      VWC_scaled = daily_weather$VWC_scaled),
                          cdate = daily_weather$date_reformat,
                          bdate = growth$date_reformat,
                          baseline = m_baseline,
                          cinterval = "month",
                          range = c(window_open_max, window_open_min),
                          #exclude = c(window_exclude_dur, window_exclude_max),                              
                          type = "absolute", refday = c(15, 06), 
                          stat = 'mean', 
                          func = c('lin'))
  
  addVars_list <- addVars(growthWin, data1 = growth, responseVar = 'area')

  m_baseline <- update( m_baseline, paste0(  ". ~ . + ", addVars_list$bestVar), 
                        data = addVars_list$data2)
  
  newVars <- as.list( daily_weather[, addVars_list$addVars])
  
  growthWin2 <- slidingwin(xvar = newVars,
                           cdate = daily_weather$date_reformat,
                           bdate = addVars_list$data2$date_reformat,
                           baseline = m_baseline,
                           cinterval = 'month',
                           range = c(window_open_max, window_open_min),
                           #exclude = c(window_exclude_dur, window_exclude_max),
                           type = "absolute", refday = c(15, 06),
                           stat = 'mean', 
                           func = c('lin'))
  
  out_obj_name <- paste(species, 'growth_small', 'monthly_ClimWin', sep = '_')
  
  out <- list( growthWin, growthWin2 )
  names(out ) <- c('ClimWinFit1', 'ClimWinFit2')
  
  assign( 
    out_obj_name, 
    out
  )
  
  save(list = out_obj_name, 
       file = paste0( "output/growth_models/", species, "_growth_small_", model_type, "_monthly_ClimWin.rda"))
  
}
