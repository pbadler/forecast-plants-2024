rm(list = ls() )
library(lme4)
library(tidyverse)
library(caret)
source('code/analysis/functions.R')

species_list <- c('ARTR')#, 'HECO', 'POSE', 'PSSP')
out <- list( list()  , list(), list() , list() )
names( out ) <- species_list
sp <- "ARTR"

for( sp in species_list) { 
  
  all_dat  <- read_csv( paste0( 'data/temp/', sp, '_growth_training_data.csv'))
  all_dat_no_intxn <- read_csv(paste0( 'data/temp/', sp, '_growth_no_intxn_training_data.csv'))
  all_dat_small <- read_csv(paste0( 'data/temp/', sp, '_growth_small_training_data.csv') ) 
  
  stopifnot( all( all_dat$pid == all_dat_no_intxn$pid) & all( all_dat$year == all_dat_no_intxn$year) ) # check that training data are for the same plants  
  
  g_model <- read_rds(paste0( 'output/growth_models/', sp, '_growth.rds'))
  g_model_null <- read_rds(paste0( 'output/growth_models/', sp, '_growth_null.rds'))  
  
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
  
  # No intxn model 
  g_model_no_intxn <- read_rds(paste0( 'output/growth_models/', sp, '_growth_no_intxn.rds'))
  
  # climate model 
  clim_no_intxn_cv <- cross_validate_growth(model = g_model_no_intxn, data = all_dat_no_intxn, folds = folds)
  clim_no_intxn_cv$mtype <- "clim no intxn"
  clim_no_intxn_cv$species <- sp 
  clim_no_intxn_cv$unit <- 'log size'
  clim_no_intxn_cv$size_class <- "large"
  
  # Small plants ------------------------------------------ 
  g_model_small <- read_rds(paste0( 'output/growth_models/', sp, '_growth_small.rds'))
  g_model_small_null <- read_rds(paste0( 'output/growth_models/', sp, '_growth_small_null.rds'))  
  
  # K-fold CV for small plants ------------------------------ # 
  years <- unique( all_dat_small$year )
  folds <- lapply( years , function( y ) { return( which(all_dat_small$year == y ))})
  
  # null model
  null_small_cv <- cross_validate_growth(model = g_model_small_null, data = all_dat_small, folds = folds)
  null_small_cv$mtype <- "null"
  null_small_cv$species <- sp 
  null_small_cv$unit <- 'log size'
  null_small_cv$size_class <- 'small'
  
  # climate model 
  clim_small_cv <- cross_validate_growth(model = g_model_small, data = all_dat_small, folds = folds)
  clim_small_cv$mtype <- "clim"
  clim_small_cv$species <- sp 
  clim_small_cv$unit <- 'log size'
  clim_small_cv$size_class <- "small"
  
  # Average scores 
  out[[sp]] <- bind_rows(
    null_cv, 
    clim_cv, 
    clim_no_intxn_cv,
    null_small_cv, 
    clim_small_cv) %>% 
    group_by( mtype, species, unit , size_class) %>% 
    summarise( R2 = mean(R2, na.rm = T), RMSE = mean(RMSE), MAE = mean(MAE))  
  
  out[[sp]]$AICc <- NA
  out[[sp]]$BIC <- NA
  
  out[[sp]]$AICc[ out[[sp]]$mtype == 'clim' & out[[sp]]$size_class == 'large' ] <- MuMIn::AICc(update( g_model, data = all_dat, REML = F))
  out[[sp]]$AICc[ out[[sp]]$mtype == 'null' & out[[sp]]$size_class == 'large' ] <- MuMIn::AICc(update( g_model_null, data = all_dat, REML = F))
  out[[sp]]$AICc[ out[[sp]]$mtype == 'clim no intxn' & out[[sp]]$size_class == 'large' ] <- MuMIn::AICc(update( g_model_no_intxn, data =all_dat_no_intxn, REML = F))
  
  out[[sp]]$BIC[ out[[sp]]$mtype == 'clim' & out[[sp]]$size_class == 'large' ] <- BIC(update( g_model, data = all_dat, REML = F))
  out[[sp]]$BIC[ out[[sp]]$mtype == 'null' & out[[sp]]$size_class == 'large' ] <- BIC(update( g_model_null, data = all_dat, REML = F))
  out[[sp]]$BIC[ out[[sp]]$mtype == 'clim no intxn' & out[[sp]]$size_class == 'large' ] <- BIC(update( g_model_no_intxn, data = all_dat_no_intxn, REML = F))
  
  
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







