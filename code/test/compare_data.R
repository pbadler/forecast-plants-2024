rm(list = ls())
library(tidyverse)

old_dir <- '~/Dropbox/projects/old_USSES_projects/'

df1 <- readRDS('data/temp_data/ARTR_growth_survival_dataframe.RDS')
df2 <- readRDS(paste0(old_dir, 'ExperimentTests/precip/data/temp_data/ARTR_growth_survival_dataframe.RDS'))

dim(df1)
dim(df2)
all.equal(df1, df2)

plot( df1$C.P.f.w.sp.1, df2$C.P.f.w.sp.1 )
plot( df1$C.T.sp.1, df2$C.T.sp.1)

ctdiff <- abs(df1$C.T.1 - df2$C.T.1)
hist(ctdiff)
df1 %>% 
  filter( ctdiff > 0.2) %>% 
  group_by( year ) %>% 
  distinct(year)

# Note: there is some difference in the climate data being used !!!!! 

acc1 <- readRDS('data/temp_data/all_clim_covs.RDS')
acc2 <- readRDS(paste0(old_dir, 'ExperimentTests/precip/data/temp_data/all_clim_covs.RDS'))

plot( acc1$P.a.1, acc2$P.a.1)
plot( acc1$T.sp.1, acc2$T.sp.1)

diffT <- abs( acc1$T.sp.1 - acc2$T.sp.1 )
hist(diffT)

acc1 %>% 
  filter( diffT > 1.0) %>% 
  group_by( year ) %>% 
  distinct(year)

# Recent years show a strong difference in temperature!!!! 

seasonal1 <- readRDS('data/temp_data/seasonal_climate.RDS')
seasonal2 <- readRDS(paste0(old_dir, 'ExperimentTests/precip/data/temp_data/seasonal_climate.RDS'))

dim(seasonal1)
dim(seasonal2)

# new data is slightly more complete 
seas1 <- 
  seasonal1 %>% 
  filter( Treatment == 'Control', var == 'TAVG_avg') %>% 
  distinct() 

seas2 <- 
  seasonal2 %>% 
  rename('YEAR' = year) %>% 
  filter( Treatment == 'Control', var == 'TAVG_avg') %>%  
  distinct()

seas_joined <- 
  seas1 %>% 
  left_join(seas2, by = c('YEAR', 'season', 'Treatment', 'var')) %>% 
  gather( type , val, val.x, val.y) %>% 
  mutate( type = factor( type, labels = c('new', 'old')))

seas_joined %>% 
  ggplot( aes( x = YEAR, y = val, color = type)) + 
  geom_line() + 
  facet_wrap(~season)

# Very weird, old data mostly new data until the 2000's 
# Also summer matches throughout! 

# Solved the recent spring and fall discrepancy.  

# What is the cause of the the seasonal difference in the fall in the 1950's? 

seas_joined %>% 
  filter( YEAR > 1940, YEAR < 1965) %>% 
  filter( season == 'fall' ) %>% 
  spread( type, val ) %>% 
  mutate( diff = new - old) %>% 
  filter( abs(diff) > 0.05)

# Check the monthly data 
monthly1 <- readRDS('data/temp_data/monthly_climate.RDS')
monthly2 <- readRDS(paste0(old_dir, 'ExperimentTests/precip/data/temp_data/monthly_climate.RDS'))

monthly2 %>% nrow()
monthly2 %>% distinct() %>% nrow() # duplicates in the data  could be throwing off the averages 

##
##  In the old data september is doubled in the 1950's and October is doubled since the 1960's 
## 


monthly2 %>% filter( year == 1952, season == 'fall')  ### duplicate rows for september !!!! This is the error 
monthly1 %>% filter( YEAR == 1952, season == 'fall')

monthly2 %>% filter( year == 1952, season == 'fall') %>% group_by( year, season) %>% summarise( TAVG = mean(TAVG, na.rm = T))
monthly1 %>% filter( YEAR == 1952, season == 'fall') %>% group_by( YEAR, season) %>% summarise( TAVG = mean(TAVG, na.rm = T))

# recreate seasonal from monthly to check 

seasonal_from_monthly_joined <- 
  monthly1 %>% 
  group_by(YEAR, season) %>% 
  summarise( new = mean(TAVG, na.rm = T)) %>% 
  left_join( 
    monthly2 %>% 
      group_by( year, season) %>% 
      summarise( old = mean(TAVG, na.rm = T)) %>% 
      rename('YEAR' = year), 
    by = c('YEAR', 'season')) %>% 
  gather( type, TAVG, new, old )


seasonal_from_monthly_joined %>% 
  ggplot( aes(x = YEAR, y = TAVG, color = type )) + geom_line() + facet_wrap(~ season )
# see the same issue 
# now do distinct on old data before aggregating 

seasonal_from_monthly_joined2 <- 
  monthly1 %>% 
  distinct() %>% 
  group_by(YEAR, season) %>% 
  summarise( new = mean(TAVG, na.rm = T)) %>% 
  left_join( 
    monthly2 %>% 
      distinct() %>% 
      group_by( year, season) %>%
      summarise( old = mean(TAVG, na.rm = T)) %>% 
      rename('YEAR' = year), 
    by = c('YEAR', 'season')) %>% 
  gather( type, TAVG, new, old )


seasonal_from_monthly_joined2 %>% 
  ggplot( aes(x = YEAR, y = TAVG, color = type )) + geom_line() + facet_wrap(~ season )
##
## FIXES the problem!!! 
##





### --------------------------------- # 
###

oldC1 <- read_csv(paste0(old_dir, 'ExperimentTests/precip/data/USSES_climate_monthly_new.csv'))
oldCZ1 <- read_csv(paste0(old_dir,'ExperimentTests/precip/data/Zachman_monthly_mean_temp.csv'))
oldCZ2 <- read_csv(paste0(old_dir, 'ExperimentTests/precip/data/Zachman_ppt.csv'))

oldC1 <- read_csv(paste0(old_dir,'ExperimentTests/precip/data/USSES_climate_monthly_new.csv'))
old_seasons <- read_csv(paste0(old_dir, 'ExperimentTests/precip/data/season_table.csv'))
new_seasons <- read_csv('data/season_table.csv')
all.equal(old_seasons, new_seasons) # seasons match ! 

# merge Zachman data 
oldCZ1 <- 
  oldCZ1 %>% 
  dplyr::select( - ANNUAL) %>% 
  gather(Month_name, TAVG, JAN:DEC) %>% 
  mutate( TAVG = (TAVG - 32)*5/9) # convert to celsius 

oldCZ2 <- 
  oldCZ2 %>% 
  dplyr::select( - ANNUAL) %>% 
  gather(Month_name, PRCP, JAN:DEC) %>% 
  mutate( PRCP = PRCP*25.4)       # convert to mm 

months <- data.frame( month = 1:12, Month_name = toupper( month.abb))

oldCZ <- merge( oldCZ1, oldCZ2, by = c('YEAR', 'Month_name'))
oldCZ <- merge( oldCZ, months)

# ---process dates----------------------------------------------------------------------#


oldC1$date <-  as.POSIXct( strptime( paste0(as.character(oldC1$DATE), '-01'), '%Y-%m-%d', tz = 'MST')  ) 
oldC1$month <- as.numeric(strftime( oldC1$date, '%m'))
oldC1$year <- as.numeric(strftime( oldC1$date, '%Y'))
oldCZ$year <- as.numeric(oldCZ$YEAR)

oldC1 %>% 
  select( DATE, date, year, month)

oldCZ %>% 
  select(YEAR, month, Month_name, year)

# set-up aggregate seasonal variables for model ----------------------------------------#

month_data <- merge( 
  oldCZ[ , c('year', 'month', 'PRCP', 'TAVG')], 
  oldC1[, c('year', 'month', 'PRCP', 'TAVG')], 
  by = c('year', 'month'), all = TRUE) 

month_data <- 
  month_data %>% 
  mutate( TAVG = ifelse(is.na(TAVG.x), TAVG.y, TAVG.x)) %>% 
  mutate( PRCP = ifelse(is.na(PRCP.x), PRCP.y, PRCP.x))

old_monthly <- merge( month_data, old_seasons, by = 'month')

old_monthly <- 
  old_monthly %>% 
  mutate( water_year = year + lag_year ) %>% 
  mutate( quarter = cut(month, 4, labels = paste0('Q', 1:4))) %>%
  dplyr::select(year, quarter, month, year, season, season_label, precip_seasons, water_year, PRCP, TAVG)

# get new monthly data to compare 
new_monthly <- readRDS('data/temp_data/monthly_climate.RDS')

new_monthly <- new_monthly %>% 
  select( YEAR, MONTH, PRCP, TAVG) %>% 
  gather( stat, value, PRCP, TAVG )

old_monthly <- old_monthly %>%
  rename('YEAR' = year, 
         'MONTH' = month) %>% 
  select( YEAR, MONTH, PRCP, TAVG ) %>% 
  gather(stat, value, PRCP, TAVG)

monthly_compare <- 
  new_monthly %>% 
  left_join(old_monthly, by = c('YEAR', 'MONTH', 'stat')) %>% 
  gather( type, value,  value.x, value.y) %>% 
  mutate( type = factor(type, labels = c('new', 'old')))

monthly_compare %>% 
  filter( stat == 'TAVG') %>% 
  ggplot( aes( x = YEAR, y = value, color = type )) +
  geom_line() + 
  facet_wrap(~ MONTH)

monthly_compare %>% 
  filter( stat == 'PRCP') %>% 
  ggplot( aes( x = YEAR, y = value, color = type )) +
  geom_line() + 
  facet_wrap(~ MONTH)

# 
rm(list = ls())
library(tidyverse)
library(lubridate)

old <- readRDS(
  '~/Dropbox/projects/old_USSES_projects/driversdata/data/idaho_modern/soil_moisture_data/data/processed_data/decagon_data_with_station_data.RDS'
)

new <- readRDS('data/temp_data/decagon_data_with_station_data.RDS')

old_daily  <- 
  old %>% 
  ungroup() %>% 
  select(datetime, plot, id, depth, measure, v) %>% 
  filter( measure != 'C') %>% 
  filter( !is.na(v)) %>% 
  mutate( date = date(datetime)) %>% 
  group_by(date, plot, id, depth, measure) %>% 
  summarise( v = mean(v, na.rm = T))

new_daily <- 
  new %>% 
  ungroup() %>% 
  select( datetime, plot, id, depth, measure, v) %>% 
  filter( measure != 'C') %>% 
  filter( !is.na(v)) %>% 
  mutate( date = date(datetime)) %>% 
  group_by( date, plot, id, depth, measure) %>% 
  summarise( v = mean(v, na.rm = T)) 

dim(old_daily)
dim( new_daily)

old_daily$plot <- paste0( 'X', old_daily$plot )
all.equal(data.frame( old_daily), data.frame( new_daily) )


