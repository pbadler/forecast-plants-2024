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

# ---------- annual soil moisture -------------------------------------------------#
annual_VWC <- 
  df %>% 
  group_by( Period, year, Treatment ) %>%
  summarise (avg_VWC = mean(VWC_raw, na.rm = TRUE))

# ---------- seasonal soil moisture  -----------------------------------------------#
seasonal_VWC <- 
  df %>% 
  mutate(year = ifelse(month == 12 , year + 1, year  )) %>% # account for December 
  group_by(Period, year, season, Treatment) %>% 
  summarise( avg = mean(VWC_raw, na.rm = TRUE) )

# ---------- seasonal soil moisture  -----------------------------------------------#
monthly_avg <- 
    df %>%
    filter( Treatment == 'Control') %>%
    group_by( month, year ) %>%
    summarise (avg_VWC = mean(VWC_raw, na.rm = TRUE) )

# -------- output -----------------------------------------------------------------------------#

saveRDS( seasonal_VWC, seasonal_outfile) 
saveRDS( annual_VWC, annual_outfile) 
saveRDS( monthly_avg, monthly_outfile)