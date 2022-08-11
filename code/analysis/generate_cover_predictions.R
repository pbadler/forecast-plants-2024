rm(list = ls() )

library(tidyverse)
library(lme4)
source('code/analysis/functions.R')

sp <- 'ARTR'
dat <- get_cover_prediction_data(sp, type = NULL )
models <- get_model_list(sp, type = NULL )

cover_obs <- 
  read_csv('data/temp/all_cover.csv') %>% 
  filter( spp == sp )

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

cover <- 
  cover_obs %>% 
  rename( 'observed' = cover ) %>% 
  left_join(all_cover_predictions)


cover %>% 
  pivot_longer( c(predicted_clim, predicted_no_intxn, predicted_null, observed)) %>% 
  filter( name %in% c( 'predicted_no_intxn', 'observed') ) %>%
  filter( year < 1960 ) %>%  
  group_by( QuadName ) %>% 
  filter( !all(is.na(value))) %>% 
  ungroup() %>% 
  ggplot( aes( x = year, y = value, color = name )) + 
  geom_line() + 
  facet_wrap( ~ QuadName ) + 
  scale_color_viridis_d() + 
  ggtitle(paste0( sp, " cover (%)"))


cover %>% 
  pivot_longer( c(predicted_clim, predicted_no_intxn, predicted_null, observed)) %>% 
  filter( year > 2009 ) %>%  
  group_by( QuadName ) %>% 
  filter( !all(is.na(value))) %>% 
  ungroup() %>% 
  ggplot( aes( x = year, y = value, color = name )) + 
  geom_line() + 
  facet_wrap( ~ QuadName ) + 
  scale_color_viridis_d() + 
  ggtitle(paste0( sp, " cover (%)"))


cover %>% 
  ggplot( aes( x = observed, y = predicted_clim) ) + 
  geom_point() + 
  geom_point( aes( y = predicted_no_intxn ), color = 'gray') + 
  geom_point( aes( y = predicted_null), color = 'lightblue') + 
  facet_wrap( ~ Treatment)


