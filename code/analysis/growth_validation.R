rm(list = ls() )
library(lme4)
library(tidyverse)
library(caret)
source('code/analysis/functions.R')

species_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')
out <- list( list()  , list(), list() , list() )
names( out ) <- species_list
sp <- "ARTR"

for( sp in species_list) { 
  
  g_model <- read_rds(paste0( 'output/growth_models/', sp, '_growth.rds'))
  g_model_null <- read_rds(paste0( 'output/growth_models/', sp, '_growth_null.rds'))  
  
  all_dat <- g_model@frame
  
  years <- unique( all_dat$year )
  folds <- lapply( years , function( y ) { return( which(all_dat$year == y ))})

  # Leave one year out cross validation of historical data 
  null_cv <- cross_validate_growth(model = g_model_null, data = all_dat, folds = folds)
  null_cv$mtype <- "null"
  null_cv$species <- sp 
  null_cv$unit <- 'log size'
  null_cv$size_class <- 'large'
  
  # climate model 
  clim_cv <- cross_validate_growth(model = g_model, data = all_dat , folds = folds)
  clim_cv$mtype <- "clim"
  clim_cv$species <- sp 
  clim_cv$unit <- 'log size'
  clim_cv$size_class <- "large"
  
  # Small plants ------------------------------------------ 
  g_model_small <- read_rds(paste0( 'output/growth_models/', sp, '_small_plant_growth.rds'))
  g_model_small_null <- read_rds(paste0( 'output/growth_models/', sp, '_small_plant_growth_null.rds'))  
  
  all_dat_small <- g_model_small@frame
  
  # K-fold CV for small plants ------------------------------ # 
  years <- unique( all_dat_small$year )
  folds <- lapply( years , function( y ) { return( which(all_dat_small$year == y ))})
  
  small_plant_null_cv <- cross_validate_growth(model = g_model_small_null, data = all_dat_small, folds = folds)
  small_plant_null_cv$mtype <- "null"
  small_plant_null_cv$species <- sp 
  small_plant_null_cv$unit <- 'log size'
  small_plant_null_cv$size_class <- 'small'
  
  # climate model 
  small_plant_clim_cv <- cross_validate_growth(model = g_model_small, data = all_dat_small , folds = folds)
  small_plant_clim_cv$mtype <- "clim"
  small_plant_clim_cv$species <- sp 
  small_plant_clim_cv$unit <- 'log size'
  small_plant_clim_cv$size_class <- "small"
  
  # Average scores 
  out[[sp]] <- bind_rows(
    null_cv, 
    clim_cv, 
    small_plant_null_cv, 
    small_plant_clim_cv) %>% 
    group_by( mtype, species, unit , size_class) %>% 
    summarise( R2 = mean(R2, na.rm = T), RMSE = mean(RMSE), MAE = mean(MAE))  
  
  out[[sp]]$AICc <- NA
  out[[sp]]$BIC <- NA
  
  out[[sp]]$AICc[ out[[sp]]$mtype == 'clim' & out[[sp]]$size_class == 'large' ] <- MuMIn::AICc(update( g_model, data = all_dat, REML = F))
  out[[sp]]$AICc[ out[[sp]]$mtype == 'null' & out[[sp]]$size_class == 'large' ] <- MuMIn::AICc(update( g_model_null, data =all_dat, REML = F))
  out[[sp]]$BIC[ out[[sp]]$mtype == 'clim' & out[[sp]]$size_class == 'large' ] <- BIC(update( g_model, data = all_dat, REML = F))
  out[[sp]]$BIC[ out[[sp]]$mtype == 'null' & out[[sp]]$size_class == 'large' ] <- BIC(update( g_model_null, data =all_dat, REML = F))
  
  
  out[[sp]]$AICc[ out[[sp]]$mtype == 'clim' & out[[sp]]$size_class == 'small' ] <- MuMIn::AICc(update( g_model_small, data = all_dat_small, REML = F))
  out[[sp]]$AICc[ out[[sp]]$mtype == 'null' & out[[sp]]$size_class == 'small' ] <- MuMIn::AICc(update( g_model_small_null, data = all_dat_small, REML = F))
  out[[sp]]$BIC[ out[[sp]]$mtype == 'clim' & out[[sp]]$size_class == 'small' ] <- BIC(update( g_model_small, data = all_dat_small, REML = F))
  out[[sp]]$BIC[ out[[sp]]$mtype == 'null' & out[[sp]]$size_class == 'small' ] <- BIC(update( g_model_small_null, data =all_dat_small, REML = F))
  
}

cv_out <- do.call( rbind, out ) 

cv_out %>% 
  mutate( Treatment = 'Control', Period = 'Historical') %>% 
  pivot_longer(cols = R2:BIC) %>% unite( col = 'stat', c(name, mtype)) %>% 
  pivot_wider(names_from = stat, values_from = value) %>% 
  write_csv('output/growth_models/cross_validation_growth_models.csv')

# Out of sample validation 
out <- list( list()  , list(), list() , list() )
names( out ) <- species_list
sp <- 'ARTR'

for( sp in species_list) { 
  
  g_model <- read_rds(paste0( 'output/growth_models/', sp, '_growth.rds'))
  g_model_null <- read_rds(paste0( 'output/growth_models/', sp, '_growth_null.rds'))  
  
  # Predict out of sample experimental data 
  # 2010 - 2016 
  testing <- read_csv( paste0( 'data/temp/', sp, '_growth_testing_data.csv'))
  
  # Climate model ------------------------  # 
  testing$yhat_clim <- predict( g_model, newdata = testing, allow.new.levels = T)
  # Null Model no climate vars ------------ # 
  testing$yhat_null <- predict( g_model_null, newdata = testing, allow.new.levels = T)
  
  oos_stats <- testing %>% 
    filter( size_class == 'large') %>% 
    group_by( size_class, Treatment ) %>% 
    summarise( 
      R2_clim = caret::R2( pred  = yhat_clim, obs = area ), 
      RMSE_clim = caret::RMSE( pred = yhat_clim, obs = area), 
      MAE_clim = caret::MAE( pred = yhat_clim, obs = area), 
      R2_null = caret::R2( pred = yhat_null, obs = area ), 
      RMSE_null = caret::RMSE( pred = yhat_null, obs = area), 
      MAE_null = caret::MAE( pred = yhat_null, obs = area )) %>% 
    bind_rows(
      testing %>% 
        group_by(size_class) %>% 
        summarise( 
          R2_clim = caret::R2( pred  = yhat_clim, obs = area ), 
          RMSE_clim = caret::RMSE( pred = yhat_clim, obs = area), 
          MAE_clim = caret::MAE( pred = yhat_clim, obs = area), 
          R2_null = caret::R2( pred = yhat_null, obs = area ), 
          RMSE_null = caret::RMSE( pred = yhat_null, obs = area), 
          MAE_null = caret::MAE( pred = yhat_null, obs = area )) %>% 
        mutate( Treatment = 'Overall')
    )
  oos_stats$species <- sp 
  oos_stats$unit <- 'log_size'
  
  # do small plants ---------------------------------
  g_model <- read_rds(paste0( 'output/growth_models/', sp, '_small_plant_growth.rds'))
  g_model_null <- read_rds(paste0( 'output/growth_models/', sp, '_small_plant_growth_null.rds')) 
  
  # testing data 
  testing <- read_csv( paste0( 'data/temp/', sp, '_small_plant_growth_testing_data.csv'))
  
  # Climate model ------------------------  # 
  testing$yhat_clim <- predict( g_model, newdata = testing, allow.new.levels = T)
  # Null Model no climate vars ------------ # 
  testing$yhat_null <- predict( g_model_null, newdata = testing, allow.new.levels = T)
  
  oos_stats_small <- testing %>% 
    filter( size_class == 'small') %>% 
    group_by( size_class, Treatment ) %>% 
    summarise( 
      R2_clim = caret::R2( pred  = yhat_clim, obs = area ), 
      RMSE_clim = caret::RMSE( pred = yhat_clim, obs = area), 
      MAE_clim = caret::MAE( pred = yhat_clim, obs = area), 
      R2_null = caret::R2( pred = yhat_null, obs = area ), 
      RMSE_null = caret::RMSE( pred = yhat_null, obs = area), 
      MAE_null = caret::MAE( pred = yhat_null, obs = area )) %>% 
    bind_rows(
      testing %>% 
        group_by(size_class) %>% 
        summarise( 
          R2_clim = caret::R2( pred  = yhat_clim, obs = area ), 
          RMSE_clim = caret::RMSE( pred = yhat_clim, obs = area), 
          MAE_clim = caret::MAE( pred = yhat_clim, obs = area), 
          R2_null = caret::R2( pred = yhat_null, obs = area ), 
          RMSE_null = caret::RMSE( pred = yhat_null, obs = area), 
          MAE_null = caret::MAE( pred = yhat_null, obs = area )) %>% 
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







