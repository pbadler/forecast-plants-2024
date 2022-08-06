rm(list = ls())
library(tidyverse)
library(lubridate)
library(climwin)
library(imputeTS)
library(sheepweather)

source('code/analysis/functions.R')

# Climate and VWC data  ------------------- # 
quad_info <- read_csv( file = 'data/quad_info.csv')
daily_vwc <- read_csv('data/temp/daily_swVWC_treatments.csv') 

daily_vwc <- 
  daily_vwc %>% 
  mutate( VWC_log = log( VWC)) %>% 
  mutate( date_reformat = strftime( date, "%d/%m/%Y"))

daily_weather <- sheepweather::usses_weather %>% 
  pivot_wider(values_from = value, names_from = ELEMENT) %>% 
  mutate_at(.vars = c('PRCP', 'TMAX', 'TMIN'), function( x ) na_interpolation(ts(x)) ) %>% 
  rowwise() %>%
  mutate( TAVG = (TMAX + TMIN)/2) %>%
  mutate( date_reformat = strftime(date, '%d/%m/%Y')) %>% 
  ungroup() %>% 
  mutate( year = year(date)) 

daily_weather <- 
  daily_vwc %>% 
  left_join(daily_weather) %>% 
  mutate_at( .vars  = c( 'PRCP', 'TMAX', 'TMIN', 'TAVG',  'VWC' ), 
             list( scaled = function(x) as.numeric(scale(x)))) 

write_csv(daily_weather, file = 'data/temp/daily_weather_for_models.csv')
