rm(list = ls() )
library(lme4)
library(tidyverse)
library(caret)
library(pROC)
source('code/analysis/functions.R')

s_model_files <- dir('output/survival_models/', 'survival.rds', full.names = T)
s_null_model_files <- dir('output/survival_models/', '_survival_null.rds', full.names = T)
survival_training_data_fls <- dir('data/temp/', 'survival_training_data.csv', full.names = T)

species_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')
out <- list( list()  , list(), list() , list() )
names( out ) <- species_list
sp <- 'ARTR'

for( sp in species_list) { 
  s_model <- read_rds( s_model_files[str_detect(s_model_files, sp)])
  s_model_null <- read_rds( s_null_model_files[str_detect(s_null_model_files, sp)])
  
  years <- unique( s_model@frame$year)
  folds <- lapply( years , function( y ) { return( which(s_model@frame$year == y ))})
  
  all_dat  <- s_model@frame
  
  # Leave one year out cross validation of historical data 
  null_cv <- cross_validate_survival(model = s_model_null, data = all_dat, folds = folds )
  
  null_cv$mtype <- "null"
  null_cv$species <- sp 
  null_cv$unit <- 'survival'
  
  clim_cv <- cross_validate_survival(model = s_model, data = all_dat, folds = folds )
  clim_cv$mtype <- "clim"
  clim_cv$species <- sp 
  clim_cv$unit <- 'survival'
  
  # Average AUC scores 
  out[[sp]] <- bind_rows(
    null_cv, 
    clim_cv) %>% 
    group_by( mtype, species, unit ) %>% 
    summarise( AUC = mean(AUC, na.rm = T), RMSE = mean(RMSE) )  
  
  out[[sp]]$AICc <- NA
  out[[sp]]$BIC <- NA
  out[[sp]]$DIC <- NA
  
  out[[sp]]$AICc[ out[[sp]]$mtype == 'clim' ] <- MuMIn::AICc(  s_model )
  out[[sp]]$AICc[ out[[sp]]$mtype == 'null' ] <- MuMIn::AICc( s_model_null )
  
  out[[sp]]$BIC[ out[[sp]]$mtype == 'clim' ] <- BIC( s_model )
  out[[sp]]$BIC[ out[[sp]]$mtype == 'null' ] <- BIC( s_model_null )

  out[[sp]]$DIC[ out[[sp]]$mtype == 'clim' ] <- MuMIn::DIC( s_model )
  out[[sp]]$DIC[ out[[sp]]$mtype == 'null' ] <- MuMIn::DIC( s_model_null )
  
}
cv_out <- do.call( rbind, out ) 

cv_out %>% 
  mutate( Treatment = 'Control', Period = 'Historical') %>% 
  pivot_longer(cols = AUC:BIC) %>% unite( col = 'stat', c(name, mtype)) %>% 
  pivot_wider(names_from = stat, values_from = value) %>% 
  write_csv('output/survival_models/cross_validation_survival_models.csv')

# Out of sample validation 
out <- list( list()  , list(), list() , list() )
names( out ) <- species_list
sp <- "POSE"
for( sp in species_list) { 
  
  s_model <- read_rds( s_model_files[str_detect(s_model_files, sp)])
  s_model_null <- read_rds( s_null_model_files[str_detect(s_null_model_files, sp)])
  
  # Predict out of sample experimental data 
  # 2010 - 2016 
  testing <- read_csv( paste0( 'data/temp/', sp, '_survival_testing_data.csv'))
  
  # Climate model ------------------------  # 
  testing$yhat_clim <- predict( s_model, newdata = testing, allow.new.levels = T, type = 'response')
  # Null Model no climate vars ------------ # 
  testing$yhat_null <- predict( s_model_null, newdata = testing, allow.new.levels = T, type = 'response')
  
  
  oos_stats <- 
    testing  %>%  
    group_by( Treatment ) %>% 
    summarise( 
      RMSE_clim = caret::RMSE( pred = yhat_clim, obs = survives), 
      RMSE_null = caret::RMSE( pred = yhat_null, obs = survives),
      AUC_clim =  as.numeric(auc(roc(survives, yhat_clim))),  # pROC package functions  
      AUC_null = as.numeric( auc(roc(survives, yhat_null)))) %>% 
    bind_rows(
      testing %>% 
        ungroup() %>% 
        summarise( 
          RMSE_clim = caret::RMSE( pred = yhat_clim, obs = survives), 
          RMSE_null = caret::RMSE( pred = yhat_null, obs = survives),
          AUC_clim =  as.numeric(auc(roc(survives, yhat_clim))),  # pROC package functions  
          AUC_null = as.numeric( auc(roc(survives, yhat_null)))) %>% 
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


read_csv('output/survival_models/oos_survival_validation.csv')





