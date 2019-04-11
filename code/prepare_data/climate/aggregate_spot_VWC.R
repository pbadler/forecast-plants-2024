# calculate treatment effects compared to control soil moisture 

rm(list = ls()) 

library( tidyverse) 

daily_clim <- readRDS('data/daily_station_dat_rainfall.RDS')
seasons <- read.csv('data/season_table.csv')
spotVWC <- readRDS('data/spring_spot_measurements.RDS')

spotVWC <- 
  spotVWC %>% 
  mutate( month = as.numeric( strftime( date, '%m'))) %>% 
  left_join(seasons, by = 'month')

spot_weights <- 
  spotVWC %>% 
  group_by( date, PrecipGroup ) %>% 
  summarise( weight = n())

spotVWC <- merge( spotVWC, daily_clim [ , c('date', 'rainfall')])

spotVWC <- 
  spotVWC %>% 
  group_by( season, date, PrecipGroup,rainfall, Treatment ) %>% 
  summarise( avg_VWC = mean(VWC, na.rm = TRUE)) %>% 
  group_by(PrecipGroup) %>% 
  mutate( avg_VWC = scale(avg_VWC, mean(avg_VWC[Treatment == 'Control'], na.rm = T), sd(avg_VWC[Treatment == 'Control'], na.rm = T))) %>%  # scale within Precip Group and Depth 
  spread( Treatment, avg_VWC) %>%
  mutate( Drought = Drought - Control, Irrigation = Irrigation - Control ) %>% 
  arrange( PrecipGroup, date) 

spotVWC <- merge( spotVWC, spot_weights)

spotVWC <- spotVWC %>% mutate( simple_date = date ) 

saveRDS(spotVWC, 'data/temp_data/spotVWC.RDS')

