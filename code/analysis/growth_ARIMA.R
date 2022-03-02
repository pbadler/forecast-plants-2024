rm(list = ls())
par(mfrow= c(1,1))
library(tidyverse)
library(imputeTS)
library(astsa)
library(forecast)
library(lme4)
library(zoo)
library(lubridate)


clim <- read_csv('lib/driversdata/data/idaho/climateData/Climate.csv')
monthly_clim <- readRDS('data/temp_data/monthly_climate.RDS')
monthly_vwc <- readRDS('data/temp_data/daily_swVWC_treatments.RDS') 

monthly_vwc <- 
  monthly_vwc %>% 
  mutate( month = lubridate::month(date), 
          year  = lubridate::year( date)) %>% 
  filter( Treatment == 'Control') %>% 
  group_by( year, month ) %>% 
  summarise( VWC = mean(VWC )) %>% 
  ungroup() 

monthly_clim <- 
  monthly_clim %>% 
  ungroup() %>% 
  rename( 'year' = YEAR, 
          'month' = MONTH ) %>% 
  select( year, month, PRCP, TAVG) %>% 
  arrange( year, month) %>% 
  mutate_at( c('PRCP', 'TAVG'), ts ) %>% 
  mutate_at( c('PRCP', 'TAVG'), na_interpolation)

# Early Adler climate definitions: 
# Precipitation is October through June  
# Temperature is April through June 

windows <- 
  expand.grid( month = 1:12, year = unique( monthly_clim$year), end = 0:24, len = 1:24 )

st_anom <- 
  monthly_clim %>% 
  left_join( monthly_vwc ) %>% 
  filter( !is.na(month))  %>% 
  filter( year <= 1960) %>% 
  select( year, month, PRCP:VWC) %>% 
  arrange( year, month ) %>% 
  left_join(windows, by = c('year', 'month')) %>% 
  gather( var, value, PRCP:VWC ) %>% 
  group_by( var, end, len) %>% 
  arrange( var, end, len, year ) %>% 
  mutate( rllm = rollmean( value, unique(end + len), fill = NA, align = 'right')) %>% 
  mutate( rllm = lag(rllm, unique( end ))) %>% 
  filter( month == 6 ) %>%  
  mutate( monitor_date = ymd( paste( year, month, 1 )), 
          start_date = monitor_date %m-% months( end + len), 
          end_date = monitor_date %m-% months( end ))

save('st_anom', 'windows', file = 'data/temp_data/st_anom.rda')

# Size data 
