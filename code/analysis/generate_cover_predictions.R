rm(list = ls() )

library(tidyverse)
library(lme4)
source('code/analysis/functions.R')


spList <- c('ARTR', 'HECO', 'POSE', 'PSSP')


for( sp in spList ){ 
  dat <- get_cover_prediction_data(sp, type = NULL )
  models <- get_model_list(sp, type = NULL )
  cover_predictions <- make_cover_predictions(sp, models, dat , 'clim', re.form = NA)
  
  null_models <- get_model_list(sp, 'null')
  cover_predictions_null <- make_cover_predictions(sp, null_models, dat, 'null', re.form = NA)

  no_intxn_models <- get_model_list(sp, 'no_intxn')
  dat <- get_cover_prediction_data( sp, type = 'no_intxn')
  cover_predictions_no_intxn <- make_cover_predictions(sp, no_intxn_models, dat, 'no_intxn', re.form = NA)

  # Join Predictions 
  all_cover_predictions <- 
    cover_predictions %>% 
    left_join(cover_predictions_no_intxn) %>% 
    left_join(cover_predictions_null) %>% 
    left_join(read_csv( 'data/quad_info.csv'))


  all_cover_predictions %>% 
    write_csv( paste0( 'data/temp/', sp, '_cover_predictions.csv'))

} 
