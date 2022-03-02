rm(list = ls())

library(tidyverse)
library(lme4)
library(zoo)
library(lubridate)
library(climwin)
library(nlme)

library(pollen)
library(glmmML)

# Climate and VWC data  ------------------- # 
quad_info <- read_csv( file = 'data/quad_info.csv')

load( file = 'data/temp_data/daily_weather.rda')
daily_vwc <- readRDS('data/temp_data/daily_swVWC_treatments.RDS') 

last_year <- 2010 # last year of training data, everything earlier is used 

gdd_low_temp <- 5
gdd_hi_temp <- 30

daily_vwc <- 
  daily_vwc %>% 
  filter( Treatment == 'Control', year(date) < last_year) %>% 
  mutate( date_reformat = paste( str_pad( day( date) , 2, pad = '0'), 
                                 str_pad( month( date ), 2, pad = '0'), 
                                 year(date), sep = '/' ))

daily_weather <- 
  daily_weather %>% 
  filter( year(date) < last_year) %>% 
  rowwise() %>% 
  mutate( date_reformat = paste( str_pad( day( date) , 2, pad = '0'), 
                                 str_pad( month( date ), 2, pad = '0'), 
                                 year(date), sep = '/' )) %>%
  mutate( TAVG = (tmax_ts + tmin_ts)/2 ) %>% 
  mutate( TAVG_0 = ifelse( TAVG <= 0, 0, TAVG)) %>%   # Temperatures below zero are set to zero
  mutate( TMAX_K = tmax_ts + 273.15) %>% 
  mutate( gdd = gdd(tmax_ts, tmin_ts, tbase = gdd_low_temp, tbase_max = gdd_hi_temp)) %>% # Growing degree days from pollen
  mutate( PRCP_pos = PRCP + 0.000001)  # ensure PRCP > 0 when log trans.

## ------------------------------------------------- 

sp_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')
iter <- 10
species <- 'PSSP'

for(species in sp_list){ 
  
  intra <- paste0( 'W.' , species )
  dat <- readRDS( paste0 ( 'data/temp_data/', paste0( species, '_survival.RDS')))
  
  dat <- 
    dat %>% 
    mutate( pid = paste( quad, trackID, sep = "_")) %>%
    mutate( year = year + 1) %>% # Add one year so that the survival "response" coincides with the current year, not the future
    filter( year < last_year )  
    
  mnyear <- min(dat$year)
  mxyear <- max(dat$year)
  pid <- unique( dat$pid )

  my_dat <- 
    dat %>% 
    dplyr::select(quad, pid, age, year, logarea, survives, starts_with('W.'))  %>% 
    mutate( W.total = W.ARTR + W.HECO + W.POSE + W.PSSP - eval(parse(text = intra)), 
            W.intra = eval(parse(text = intra)))
  
  survival <- 
    my_dat %>% 
    group_by( pid ) %>% 
    arrange(pid, year) %>% 
    filter( !is.na(logarea), 
            !is.na(survives),
            is.finite(logarea), 
            year > 1925) %>%
    mutate( date = paste0( '15/06/', year)) # date for climwin 
  
  survival <- 
    survival %>% 
    left_join(quad_info, by = 'quad') %>%
    mutate( group_year = paste( Group, year , sep = '_'))
  
  survivesWin_VWC <- slidingwin(xvar = list(VWC = daily_vwc$VWC),
                              cdate = daily_vwc$date_reformat,
                              bdate = survival$date,
                              baseline = glm(survives ~ 1 + logarea, 
                                                  data = survival, 
                                                  family = 'binomial'), 
                              cinterval = 'week',
                              range = c(66, 0),
                              type = "absolute", refday = c(15, 06),
                              stat = 'mean', 
                              func = c('lin', 'log'))
  
  
  survivesWin_TEMP <- slidingwin(xvar = list(TMAX = daily_weather$TMAX, 
                                             TAVG = daily_weather$TAVG, 
                                             TAVG_0 = daily_weather$TAVG_0),
                               cdate = daily_weather$date_reformat,
                               bdate = survives$date,
                               baseline = glm(survives ~ 1 + logarea, 
                                              data = survival, 
                                              family = 'binomial'),
                               cinterval = "week",
                               range = c(66, 0),
                               type = "absolute", refday = c(15, 06), 
                               stat = 'mean', 
                               func = c('lin'))
  
  survivesWin_GDD <- slidingwin(xvar = list(GDD = daily_weather$gdd),
                              cdate = daily_weather$date_reformat,
                              bdate = survives$date,
                              baseline = lme(survives ~ 1, data = survives, random = ~ 1|group_year),
                              cinterval = "week",
                              range = c(66, 0),
                              type = "absolute", refday = c(15, 06), 
                              stat = 'sum', 
                              func = c('lin', 'log', 'quad'))
  
  
  survivesWin_PRCP <- slidingwin(xvar = list(PRCP = daily_weather$PRCP_pos),
                               cdate = daily_weather$date,
                               bdate = survives$date,
                               baseline = lme(survives ~ 1, data = survives, random = ~ 1|group_year),
                               cinterval = "week",
                               range = c(66, 0),
                               type = "absolute", refday = c(15, 06),
                               stat = "sum",
                               func = c("lin", "log"))
  
  
  save(survivesWin_VWC, 
       survivesWin_TEMP, 
       survivesWin_GDD, 
       survivesWin_PRCP, 
       file = paste0( "data/temp_data/", species, "_weekly_ClimWin.rda"))
}
