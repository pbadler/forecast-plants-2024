rm(list = ls())
library(tidyverse)
library(lubridate)
library(lme4)
library(sheepweather)

# output --------------------------------------------------------------- # 
outfile <- 'data/temp_data/daily_swVWC_treatments.RDS'

# input ---------------------------------------------------------------- # 
rain_data <- readRDS('data/temp_data/daily_station_dat_rainfall.RDS')
sw <- sheepweather::usses_soilwat
sm_model <- readRDS('output/treatment_effects_on_soil_moisture.RDS')
seasons <- read_csv('data/season_table.csv')

# ----------------------------------------------------------------------- # 

sw <- 
  sw %>% 
  mutate( date = date(parse_date_time(paste( Year, Day, sep = '-') , orders = '%Y-%j'))) %>% 
  rename( 'year' = Year) %>% 
  mutate( month = month( date )) %>%
  gather( layer, VWC, starts_with('Lyr')) 

sw_avg <- 
  sw %>% 
  filter( layer %in% c('Lyr_1')) %>% # only use top layer ( upper 20 cm ) 
  group_by(date) %>%
  summarise(swVWC = mean(VWC) * 100) %>% 
  ungroup() 

rain_data <- 
  rain_data %>%
  ungroup() %>% 
  select( date, inv_days_since_rain) 

pred_df <- 
  rain_data %>% 
  mutate(Treatment = 'Control') %>% 
  bind_rows(rain_data %>% 
              mutate( Treatment = 'Irrigation')) %>% 
  bind_rows(rain_data %>% 
              mutate( Treatment = 'Drought')) %>% 
  mutate(month = month(date)) %>% 
  left_join(seasons, by = 'month') %>% 
  data.frame()

pred_df$pred <- predict( sm_model, newdata = pred_df, re.form = NA)

sw_treatment <- 
  sw_avg %>%
  left_join(pred_df, by ='date') %>% 
  mutate( VWC = swVWC*pred)

saveRDS(sw_treatment, outfile)


