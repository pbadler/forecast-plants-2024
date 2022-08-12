rm(list = ls())
library(tidyverse)
source('code/analysis/functions.R')

sp <- 'PSSP'
split_year <- 2010
spList <- c('ARTR', 'HECO', 'POSE', 'PSSP')

for ( sp in spList ) { 
  cover_predictions <- read_csv(paste0( 'data/temp/', sp, '_cover_predictions.csv'))

  cover <- 
    read_csv('data/temp/all_cover.csv') %>% 
    filter( spp == sp ) %>% 
    rename( 'observed' = cover ) %>% 
    left_join(cover_predictions)

  cover %>% 
    mutate( Period = factor(ifelse(year > split_year, 'Testing', 'Training' ))) %>%  
    pivot_longer( c(predicted_clim, predicted_no_intxn, predicted_null) ) %>% 
    filter( !( Treatment != 'Control' & Period == 'Training')) %>% 
    group_by( Period, Treatment, name ) %>%
    summarise( RMSE = caret::RMSE(value , observed, na.rm = T), 
               MAE = caret::MAE(value, observed, na.rm = T), 
               R2 = caret::R2( value, observed, na.rm = T))  %>% 
    mutate( species = sp ) %>% 
    write_csv( paste0 ( 'output/', sp, '_cover_validation.csv'))
}

#     
# cover %>% 
#   filter( name %in% c( 'predicted_no_intxn', 'observed') ) %>%
#   filter( year < 1960 ) %>%  
#   group_by( QuadName ) %>% 
#   filter( !all(is.na(value))) %>% 
#   ungroup() %>% 
#   ggplot( aes( x = year, y = value, color = name )) + 
#   geom_line() + 
#   facet_wrap( ~ QuadName ) + 
#   scale_color_viridis_d() + 
#   ggtitle(paste0( sp, " cover (%)"))
# 
# 
# cover %>% 
#   pivot_longer( c(predicted_clim, predicted_no_intxn, predicted_null, observed)) %>% 
#   filter( year > 2009 ) %>%  
#   group_by( QuadName ) %>% 
#   filter( !all(is.na(value))) %>% 
#   ungroup() %>% 
#   ggplot( aes( x = year, y = value, color = name )) + 
#   geom_line() + 
#   facet_wrap( ~ QuadName ) + 
#   scale_color_viridis_d() + 
#   ggtitle(paste0( sp, " cover (%)"))
# 
# 
# cover %>% 
#   ggplot( aes( x = observed, y = predicted_clim) ) + 
#   geom_point() + 
#   geom_point( aes( y = predicted_no_intxn ), color = 'gray') + 
#   geom_point( aes( y = predicted_null), color = 'lightblue') + 
#   facet_wrap( ~ Treatment)
# 
# 
