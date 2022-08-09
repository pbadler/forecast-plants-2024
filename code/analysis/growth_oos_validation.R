rm(list = ls() )
library(lme4)
library(tidyverse)
library(caret)
source('code/analysis/functions.R')

species_list <- c('ARTR')#, 'HECO', 'POSE', 'PSSP')
out <- list( list()  , list(), list() , list() )
names( out ) <- species_list
sp <- "ARTR"

# Out of sample validation 
out <- list( list()  , list(), list() , list() )
names( out ) <- species_list
sp <- 'ARTR'

for( sp in species_list) { 
  
  g_model <- read_rds(paste0( 'output/growth_models/', sp, '_growth.rds'))
  g_model_null <- read_rds(paste0( 'output/growth_models/', sp, '_growth_null.rds'))  
  g_model_no_intxn <- read_rds(paste0('output/growth_models/', sp, '_growth_no_intxn.rds'))
  
  # Predict out of sample experimental data 
  # 2010 - 2016 
  testing <- read_csv( paste0( 'data/temp/', sp, '_growth_testing_data.csv'))
  testing_no_intxn <- read_csv(paste0( 'data/temp/', sp, '_growth_no_intxn_testing_data.csv'))
  
  # Check that testing data include all the same obs 
  stopifnot( nrow( testing_no_intxn ) == nrow(testing))
  stopifnot( all( testing_no_intxn$pid == testing$pid) & all( testing_no_intxn$year == testing$year))
  stopifnot( min( testing$year) > 2009 )  
  
  # Climate model ------------------------  # 
  testing$yhat_clim <- predict( g_model, newdata = testing, allow.new.levels = T)
  
  # Climate no intxn -----------------------# 
  testing_no_intxn$yhat_clim_no_intxn <- predict( g_model_no_intxn, newdata = testing_no_intxn, allow.new.levels = T)
  
  # Null Model no climate vars ------------ # 
  testing$yhat_null <- predict( g_model_null, newdata = testing, allow.new.levels = T)
  
  predictions <- 
    testing %>% 
    select( size_class, pid, year, Treatment, Group, quad, area , area0 , yhat_clim, yhat_null ) %>% 
    left_join( testing_no_intxn %>% select(size_class, pid, year, Treatment, area, yhat_clim_no_intxn)) %>% 
    pivot_longer( yhat_clim:yhat_clim_no_intxn )
  
  oos_stats <- 
    predictions %>% 
    group_by( size_class, Treatment, name ) %>%
    summarise( 
      R2 = caret::R2( pred  = value, obs = area ), 
      RMSE = caret::RMSE( pred = value, obs = area), 
      MAE = caret::MAE( pred = value, obs = area)) %>% 
    bind_rows(
      predictions %>% 
        group_by(size_class , name ) %>% 
        summarise( 
          R2 = caret::R2( pred  = value , obs = area ), 
          RMSE = caret::RMSE( pred = value, obs = area), 
          MAE = caret::MAE( pred = value, obs = area)) %>% 
        mutate( Treatment = 'Overall')
    )
  
  oos_stats$species <- sp 
  oos_stats$unit <- 'log_size'
  
  # do small plants ---------------------------------
  g_model_small <- read_rds(paste0( 'output/growth_models/', sp, '_growth_small.rds'))
  g_model_small_null <- read_rds(paste0( 'output/growth_models/', sp, '_growth_small_null.rds')) 
  
  # testing data 
  testing_small <- read_csv( paste0( 'data/temp/', sp, '_growth_small_testing_data.csv'))
  
  # Climate model ------------------------  # 
  testing_small$yhat_clim <- predict( g_model_small, newdata = testing_small, allow.new.levels = T)
  # Null Model no climate vars ------------ # 
  testing_small$yhat_null <- predict( g_model_small_null, newdata = testing_small, allow.new.levels = T)
  
  predictions_small <- 
    testing_small %>% 
    select( size_class, pid, year, Treatment, Group, quad, area , area0 , yhat_clim, yhat_null ) %>% 
    pivot_longer( yhat_clim:yhat_null )
  
  oos_stats_small <- 
    predictions_small %>% 
    filter( size_class == 'small') %>% 
    group_by( size_class, Treatment , name ) %>% 
    summarise( 
      R2= caret::R2( pred  = value, obs = area ), 
      RMSE = caret::RMSE( pred = value, obs = area), 
      MAE = caret::MAE( pred = value , obs = area)) %>% 
    bind_rows(
      predictions_small %>% 
        group_by(size_class, name ) %>% 
        summarise( 
          R2= caret::R2( pred  = value, obs = area ), 
          RMSE = caret::RMSE( pred = value, obs = area), 
          MAE = caret::MAE( pred = value, obs = area)) %>% 
        mutate( Treatment = 'Overall')
    )
  
  oos_stats_small$species <- sp 
  oos_stats_small$unit <- 'log_size'
  
  out[[sp]] <- bind_rows( oos_stats, oos_stats_small)
  
}

oos_performance <- do.call(bind_rows, out)

oos_performance %>% 
  mutate( Period = 'Experimental') %>% 
  write_csv('output/growth_models/oos_growth_validation.csv')

