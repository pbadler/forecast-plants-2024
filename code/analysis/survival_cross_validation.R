rm(list = ls() )
library(lme4)
library(tidyverse)
library(caret)
library(pROC)
source('code/analysis/functions.R')

species_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')
out <- list( list()  , list(), list() , list() )
names( out ) <- species_list
sp <- 'ARTR'

for( sp in species_list) { 
  
  all_dat  <- read_csv( paste0( 'data/temp/', sp, '_survival_training_data.csv'))
  all_dat_no_intxn <- read_csv(paste0( 'data/temp/', sp, '_survival_no_intxn_training_data.csv'))

  stopifnot( all( all_dat$pid == all_dat_no_intxn$pid) & all( all_dat$year == all_dat_no_intxn$year) ) # check that training data are for the same plants  
  
  s_model <- read_rds(paste0( 'output/survival_models/', sp, '_survival.rds'))
  s_model_null <- read_rds(paste0( 'output/survival_models/', sp, '_survival_null.rds'))  
  s_model_no_intxn <- read_rds(paste0( 'output/survival_models/', sp, '_survival_no_intxn.rds'))  
  
  years <- unique( all_dat$year )
  folds <- lapply( years , function( y ) { return( which(all_dat$year == y ))})
  # Leave one year out cross validation of historical data 
  null_cv <- cross_validate_survival(model = s_model_null, data = all_dat, folds = folds)
  null_cv$mtype <- "null"
  null_cv$species <- sp 
  null_cv$unit <- 'survival'
  
  # climate model 
  clim_cv <- cross_validate_survival(model = s_model, data = all_dat , folds = folds)
  clim_cv$mtype <- "clim"
  clim_cv$species <- sp 
  clim_cv$unit <- 'survival'

  # climate no_intxn model 
  clim_no_intxn_cv <- cross_validate_survival(model = s_model_no_intxn, 
                                              data = all_dat_no_intxn, folds = folds)
  clim_no_intxn_cv$mtype <- "clim no intxn"
  clim_no_intxn_cv$species <- sp 
  clim_no_intxn_cv$unit <- 'survival'
  
  
  # Average AUC scores 
  out[[sp]] <- bind_rows(
    null_cv, 
    clim_cv, 
    clim_no_intxn_cv) %>% 
    group_by( mtype, species, unit ) %>% 
    summarise( AUC = mean(AUC, na.rm = T), RMSE = mean(RMSE) )  
  
  out[[sp]]$AICc <- NA
  out[[sp]]$BIC <- NA
  out[[sp]]$DIC <- NA
  
  out[[sp]]$AICc[ out[[sp]]$mtype == 'clim' ] <- MuMIn::AICc(  s_model )
  out[[sp]]$AICc[ out[[sp]]$mtype == 'null' ] <- MuMIn::AICc( s_model_null )
  out[[sp]]$AICc[ out[[sp]]$mtype == 'clim no intxn' ] <- MuMIn::AICc( s_model_no_intxn )
  
  out[[sp]]$BIC[ out[[sp]]$mtype == 'clim' ] <- BIC( s_model )
  out[[sp]]$BIC[ out[[sp]]$mtype == 'null' ] <- BIC( s_model_null )
  out[[sp]]$BIC[ out[[sp]]$mtype == 'clim no intxn' ] <- BIC( s_model_no_intxn )
  
  out[[sp]]$DIC[ out[[sp]]$mtype == 'clim' ] <- MuMIn::DIC( s_model )
  out[[sp]]$DIC[ out[[sp]]$mtype == 'null' ] <- MuMIn::DIC( s_model_null )
  out[[sp]]$DIC[ out[[sp]]$mtype == 'clim no intxn' ] <- MuMIn::DIC( s_model_no_intxn )
  
}
cv_out <- do.call( rbind, out ) 

cv_out %>% 
  mutate( Treatment = 'Control', Period = 'Historical') %>% 
  pivot_longer(cols = AUC:DIC) %>% unite( col = 'stat', c(name, mtype)) %>% 
  pivot_wider(names_from = stat, values_from = value) %>% 
  write_csv('output/survival_models/cross_validation_survival_models.csv')
