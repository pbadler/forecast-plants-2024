rm(list = ls() )
library(lme4)
library(tidyverse)
library(caret)
source('code/analysis/functions.R')

g_model_files <- dir('output/growth_models/', 'growth.rds', full.names = T)
g_null_model_files <- dir('output/growth_models/', 'null.rds', full.names = T)
growth_training_data_fls <- dir('data/temp/', 'growth_training_data.csv', full.names = T)

species_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')
out <- list( list()  , list(), list() , list() )
names( out ) <- species_list
sp <- "ARTR"

for( sp in species_list) { 

  g_model <- read_rds( g_model_files[str_detect(g_model_files, sp)])
  g_model_null <- read_rds( g_null_model_files[str_detect(g_null_model_files, sp)])
  
  all_dat <- g_model@frame
  
  years <- unique( all_dat$year )
  folds <- lapply( years , function( y ) { return( which(all_dat$year == y ))})

  # Leave one year out cross validation of historical data 
  null_cv <- cross_validate_growth(model = g_model_null, data = all_dat, folds = folds)
  null_cv$mtype <- "null"
  null_cv$species <- sp 
  null_cv$unit <- 'log size'
  
  # climate model 
  clim_cv <- cross_validate_growth(model = g_model, data = all_dat , folds = folds)
  clim_cv$mtype <- "clim"
  clim_cv$species <- sp 
  clim_cv$unit <- 'log size'

  # Average scores 
  out[[sp]] <- bind_rows(
    null_cv, 
    clim_cv) %>% 
    group_by( mtype, species, unit ) %>% 
    summarise( R2 = mean(R2), RMSE = mean(RMSE), MAE = mean(MAE))  
  
  out[[sp]]$AIC <- NA
  out[[sp]]$BIC <- NA
  
  out[[sp]]$AIC[ out[[sp]]$mtype == 'clim' ] <- AIC( g_model )
  out[[sp]]$AIC[ out[[sp]]$mtype == 'null' ] <- AIC( g_model_null )
  
  out[[sp]]$BIC[ out[[sp]]$mtype == 'clim' ] <- BIC( g_model )
  out[[sp]]$BIC[ out[[sp]]$mtype == 'null' ] <- BIC( g_model_null )
  
}

cv_out <- do.call( rbind, out ) 

cv_out %>% 
  mutate( Treatment = 'Control', Period = 'Historical') %>% 
  pivot_longer(cols = R2:BIC) %>% unite( col = 'stat', c(name, mtype)) %>% 
  pivot_wider(names_from = stat, values_from = value) %>% 
  write_csv('output/growth_models/in_sample_cv_growth_models.csv')

# Out of sample validation 
out <- list( list()  , list(), list() , list() )
names( out ) <- species_list

for( sp in species_list) { 
  
  g_model <- read_rds( g_model_files[str_detect(g_model_files, sp)])
  g_model_null <- read_rds( g_null_model_files[str_detect(g_null_model_files, sp)])
  
  # Predict out of sample experimental data 
  # 2010 - 2016 
  testing <- read_csv( paste0( 'data/temp/', sp, '_growth_testing_data.csv'))
  
  # Climate model ------------------------  # 
  testing$yhat_clim <- predict( g_model, newdata = testing, allow.new.levels = T)
  # Null Model no climate vars ------------ # 
  testing$yhat_null <- predict( g_model_null, newdata = testing, allow.new.levels = T)
  
  oos_stats <- testing %>% 
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
  
  out[[sp]] <- oos_stats 
  
}

oos_performance <- do.call(bind_rows, out)

oos_performance %>% 
  mutate( Period = 'Experimental') %>% 
  write_csv('output/oos_growth_validation.csv')







