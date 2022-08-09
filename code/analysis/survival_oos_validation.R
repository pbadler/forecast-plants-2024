rm(list = ls() )
library(lme4)
library(tidyverse)
library(caret)
library(pROC)
source('code/analysis/functions.R')

species_list <- c('ARTR')#, 'HECO', 'POSE', 'PSSP')
out <- list( list()  , list(), list() , list() )
names( out ) <- species_list
sp <- 'ARTR'

# Out of sample validation 
out <- list( list()  , list(), list() , list() )
names( out ) <- species_list
sp <- "ARTR"

for( sp in species_list) { 
  
  s_model <- read_rds(paste0( 'output/survival_models/', sp, '_survival.rds'))
  s_model_null <- read_rds(paste0( 'output/survival_models/', sp, '_survival_null.rds'))  
  s_model_no_intxn <- read_rds(paste0('output/survival_models/', sp, '_survival_no_intxn.rds'))
  
  # Predict out of sample experimental data 
  # 2010 - 2016 
  testing <- read_csv( paste0( 'data/temp/', sp, '_survival_testing_data.csv'))
  testing_no_intxn <- read_csv(paste0( 'data/temp/', sp, '_survival_no_intxn_testing_data.csv'))
  
  # Check that testing data include all the same obs 
  stopifnot( nrow( testing_no_intxn ) == nrow(testing))
  stopifnot( all( testing_no_intxn$pid == testing$pid) & all( testing_no_intxn$year == testing$year))
  stopifnot( min( testing$year) > 2009 )  
  
  # Climate model ------------------------  # 
  testing$yhat_clim <- predict( s_model, newdata = testing, allow.new.levels = T, type = 'response')
  
  # Climate no intxn -----------------------# 
  testing_no_intxn$yhat_clim_no_intxn <- predict( s_model_no_intxn, newdata = testing_no_intxn, allow.new.levels = T, type = 'response')
  
  # Null Model no climate vars ------------ # 
  testing$yhat_null <- predict( s_model_null, newdata = testing, allow.new.levels = T, type = 'response')
  
  predictions <- 
    testing %>% 
    select( size_class, pid, year, Treatment, Group, quad, survives, area0 , yhat_clim, yhat_null ) %>% 
    left_join( testing_no_intxn %>% select(size_class, pid, year, Treatment, survives, area0, yhat_clim_no_intxn)) %>%   
    pivot_longer( yhat_clim:yhat_clim_no_intxn )
  
  
  predictions %>% 
    filter( name == 'yhat_clim') %>% 
    summarise( AUC = as.numeric( pROC::auc( pROC::roc(survives, value))))
  
  
  oos_stats <- 
    predictions %>% 
    group_by( Treatment, name ) %>%
    summarise( 
      AUC = as.numeric( pROC::auc( pROC::roc(survives, value))), 
      RMSE = caret::RMSE( pred = value, obs = survives) ) %>% 
    bind_rows(
      predictions %>% 
        group_by( name ) %>% 
        summarise( 
          AUC = as.numeric( pROC::auc( pROC::roc(survives, value))), 
          RMSE = caret::RMSE( pred = value, obs = survives) ) %>% 
        mutate( Treatment = 'Overall')
    )
  
  oos_stats$species <- sp 
  oos_stats$unit <- 'survival'
  
  out[[sp]] <- oos_stats
  
}

oos_performance <- do.call(bind_rows, out)

oos_performance %>% 
  mutate( Period = 'Experimental') %>% 
  write_csv('output/survival_models/oos_survival_validation.csv')
