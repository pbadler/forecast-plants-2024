#######################################################################################
#
# Setup seasonal climate variables for demographic rate models 
# 
#  This script aggregates the daily weather data to monthly averages for TAVG and 
#  monthly total PRCP.  These are then aggregated to seasonal average TAVG and total 
#  PRCP.  
#
#######################################################################################

rm(list = ls()) 

library(tidyverse)
library(lubridate)
library(sheepweather)

# input --------------------------------------------------------------------#

season <- read_csv('data/season_table.csv')

weather <- usses_weather

# output --------------------------------------------------------------------#

seasonal_output <- 'data/temp/seasonal_climate.csv' # need this one for prepare climate coveriates script 
monthly_output <- 'data/temp/monthly_climate.csv' # used by plot modern climate 
annual_output <- 'data/temp/annual_climate.csv' # used by plot modern climate 

# --------------------------------------------------------------------#

dm <- c(4, 11 ) # Drought months 
im <- c(5, 10 ) # Irrigation months

p.treatments <- c(0.5, 1.5) # Drought and Irrigation adjustments to precip 
t.treatments <- c(1, 1)     # Drought and Irrigation adjustments to temperature 

# make time periods --------------------------------------------------------------------

p1 <- data.frame( Period = 'Modern', YEAR = 2007:2016)
p2 <- data.frame( Period = 'not monitored', YEAR = 1958:2006)
p3 <- data.frame( Period = 'Historical', YEAR = 1925:1957)
periods <- data.frame( rbind( p1, p2, p3 )) 

# -------------------------------------------------------------------- # 
monthly <- 
  weather %>% 
  spread( ELEMENT, value ) %>% 
  mutate( tmean = (TMAX + TMIN) / 2 ) %>% 
  mutate( MONTH = month(date)) %>%                                      ### I think this is where I made the timezone mistake
  mutate( YEAR = year(date)) %>% 
  group_by(YEAR, MONTH) %>% 
  summarise( PRCP = sum(PRCP, na.rm = T), TAVG = mean(tmean, na.rm = T))

monthly <- 
  monthly %>% 
  left_join(season, 
            by = c('MONTH' = 'month')) %>% 
  mutate( water_year = YEAR + lag_year ) %>% 
  select(YEAR, MONTH, season, precip_seasons, water_year, PRCP, TAVG)

# annual climate  -------------------------------------------------#

annual_clim <-
  monthly %>% 
  group_by( YEAR ) %>%
  summarise (TAVG = mean(TAVG, na.rm = T), 
             PRCP = sum( PRCP, na.rm = T), n = n())

# ---------- monthly climate  -----------------------------------------------------------#
# 
# incorporate Drought and Irrigation effects only in specificied months 
# 

monthly_clim <- 
  monthly %>% 
  ungroup() %>%  
  mutate(Control = PRCP) %>% 
  dplyr::select(-PRCP) %>% 
  mutate( Drought    = ifelse( YEAR > 2011 & MONTH %in% c(dm[1]:dm[2]), Control*p.treatments[1], Control) ) %>% 
  mutate( Irrigation = ifelse( YEAR > 2011 & MONTH %in% c(im[1]:im[2]), Control*p.treatments[2], Control) ) %>% 
  gather(Treatment, PRCP, Control, Drought, Irrigation)

monthly_clim %>% filter( is.na(TAVG))
monthly_clim %>% filter( is.na(PRCP))
# -------------------------------------------------------------------------------------------#
# -------------- aggregate monthly climate with Treatment effects by season -----------------#

seasonal_clim <-
  monthly_clim %>% 
  mutate(YEAR = ifelse(MONTH == 12 , YEAR + 1, YEAR  )) %>% # account for December 
  gather( var, val , TAVG, PRCP ) %>% 
  group_by(Treatment, var, MONTH)  %>% 
  mutate( val = ifelse(is.na(val), mean(val, na.rm = TRUE), val)) %>% # !!!!!!! fill in missing monthly averages after 1925 with monthly average !!!!!!!! 
  group_by( Treatment, var, YEAR, season ) %>%                        # !!!!!!! note missing values TAVG  !!!!!!!!
  summarise( avg = mean(val), ttl = sum(val) ) %>% 
  group_by(var) %>% 
  gather( stat, val, avg, ttl ) %>% 
  filter( (var == 'PRCP' & stat == 'ttl')| (var == 'TAVG' & stat == 'avg')) %>% 
  group_by( Treatment, var, stat) %>% 
  arrange(YEAR, season) %>%
  ungroup() %>% 
  unite(var, c(var, stat) , sep = '_')

# -------join periods -------------------------------------------------------------------------#

seasonal_clim <- left_join( seasonal_clim, periods , by = "YEAR")
monthly_clim <- left_join(monthly, periods, by = "YEAR" ) 
annual_clim <- left_join( annual_clim, periods, by = "YEAR" ) 

# -------- output -----------------------------------------------------------------------------#

write_csv( seasonal_clim, file = seasonal_output)
write_csv( monthly_clim, file = monthly_output)
write_csv( annual_clim, file = annual_output)
