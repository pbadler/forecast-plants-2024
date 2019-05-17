#######################################################################################
#
# Setup seasonal SoilWAT variables for demographic rate models 
#
#######################################################################################

rm(list = ls()) 

library(tidyverse)
library(zoo)

# input ---------------------------------------------------- # 

load('code/figure_scripts/my_plotting_theme.Rdata')

df <- readRDS('data/temp_data/daily_swVWC_treatments.RDS')
  # comes from the soilMoistureTreatmentEffects script 

# output ---------------------------------------------------- # 

seasonal_outfile <- 'data/temp_data/seasonal_VWC.RDS'
annual_outfile <- 'data/temp_data/annual_VWC.RDS'
monthly_outfile <- 'data/temp_data/monthly_avg.RDS'

# make time periods --------------------------------------------------------------------

p1 <- data.frame( Period = 'Modern', year = 2007:2016)
p2 <- data.frame( Period = 'not monitored', year = 1958:2006)
p3 <- data.frame( Period = 'Historical', year = 1925:1957)
periods <- data.frame( rbind( p1, p2, p3 )) 

# set-up aggregate seasonal variables for model ----------------------------------------#

df <- 
  df %>% 
  ungroup() %>% 
  mutate( year = year(date)) %>% 
  mutate( water_year = year + lag_year ) %>% 
  dplyr::select(year, month, season, season_label, precip_seasons, water_year, Treatment, date, VWC)

# ---------- annual soil moisture -------------------------------------------------#
annual_VWC <- 
  df %>% 
  group_by( year, Treatment ) %>%
  summarise (avg_VWC = mean(VWC, na.rm = TRUE)) %>% 
  left_join(periods, by = 'year')

# ---------- seasonal soil moisture  -----------------------------------------------#
seasonal_VWC <- 
  df %>% 
  mutate(year = ifelse(month == 12 , year + 1, year  )) %>% # December of year x is in winter of year x + 1 
  group_by(year, season, Treatment) %>% 
  summarise( avg = mean(VWC, na.rm = TRUE) ) %>% 
  left_join(periods, by = 'year')

# ---------- seasonal soil moisture  -----------------------------------------------#
monthly_avg <- 
    df %>%
    filter( Treatment == 'Control') %>%
    group_by( month, year ) %>%
    summarise (avg_VWC = mean(VWC, na.rm = TRUE) ) %>% 
    left_join(periods, by = 'year')

# -------- output -----------------------------------------------------------------------------#

saveRDS( seasonal_VWC, seasonal_outfile) 
saveRDS( annual_VWC, annual_outfile) 
saveRDS( monthly_avg, monthly_outfile)
