rm(list = ls())
library(lmer)
library(tidyverse)
source('code/analysis/functions.R')

max_year <- 2000
daily_vwc <- read_csv('data/temp/daily_swVWC_treatments.csv') %>% 
  mutate( date_reformat = strftime( date, "%d/%m/%Y")) %>% 
  filter( year( date) < max_year) %>% 
  filter( Treatment == 'Control')

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
  left_join(daily_weather) 

# Construct Singlewin dataframe to test climate aggregation 
size <- read_csv( 'data/temp/ARTR_size.csv') %>% filter( year < 2000 , year > 1927)
size <- size %>%
  mutate( ref_day = lubridate::mdy( paste( '06-15', year, sep = '-'))) %>% 
  mutate( date_reformat = strftime( ref_day, "%d/%m/%Y"))

size <- size %>% 
  group_by( pid ) %>% 
  arrange( year ) %>% 
  mutate( area = log(area)) %>% 
  mutate( area0 = lag(area)) %>% 
  filter( !is.na( area), !is.na( area0))

model <- lmer( area ~ area0 + (1|year) , data = size )

open <- 10
close <- 1

test <- singlewin(
  list( TMAX = daily_weather$TMAX) , 
  daily_weather$date_reformat, 
  bdate = size$date_reformat, 
  baseline = model, 
  range = c(open, close), 
  cinterval = 'month',
  stat = 'mean', 
  func = 'lin', 
  refday = c(15, 06),
  type = 'absolute')

test_case <- 
  test$BestModelData %>% 
  mutate( climate = as.numeric(climate)) %>% 
  distinct( year, climate )

results <- getWindowAvg(daily_weather, "TMAX", open, close)

compare <- test_case %>% 
  left_join(
    results %>% distinct( open_date, close_date, year, TMAX) )

all.equal(compare$climate, compare$TMAX)

compare %>%
  ggplot( aes( x = climate, y = TMAX))  + 
  geom_point() + 
  xlab( 'ClimWin calculation') + 
  ylab( 'My calculation') 

