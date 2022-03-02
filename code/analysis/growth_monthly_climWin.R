rm(list = ls())

library(tidyverse)
library(lme4)
library(zoo)
library(lubridate)
library(climwin)
library(imputeTS)

load( file = 'data/temp_data/st_anom.rda')

# Climate and VWC data  ------------------- # 
monthly_clim <- read_rds('data/temp_data/monthly_climate.RDS')
daily_vwc <- readRDS('data/temp_data/daily_swVWC_treatments.RDS') 

monthly_vwc <- 
  daily_vwc %>% 
  mutate( month = lubridate::month(date), 
          year  = lubridate::year( date)) %>% 
  filter( Treatment == 'Control') %>% 
  group_by( year, month ) %>% 
  summarise( VWC = mean(VWC )) %>% 
  ungroup() 

rm(daily_vwc)

monthly_clim <- 
  monthly_clim %>% 
  ungroup() %>% 
  rename( 'year' = YEAR, 
          'month' = MONTH ) %>% 
  select( year, month, PRCP, TAVG) %>% 
  arrange( year, month) %>% 
  mutate_at( c('PRCP', 'TAVG'), ts ) %>% 
  mutate_at( c('PRCP', 'TAVG'), na_interpolation)

monthly_clim <- 
  monthly_clim %>% 
  left_join( monthly_vwc ) %>% 
  filter( !is.na(month))  %>% 
  filter( year <= 1960) %>% 
  select( year, month, PRCP:VWC) %>% 
  arrange( year, month ) 

## ------------------------------------------------- 
# ClimWin Analysis to find best window for growth 

monthly_clim <- 
  monthly_clim %>% 
  group_by( month) %>% 
  mutate( TAVG_K = TAVG + 273.15) %>% 
  mutate( PRCP_pos = PRCP + 0.00001) %>%  # ensure PRCP > 0 for log trans.
  mutate( date = paste0( '15/', str_pad( month, 2, pad = 0), '/', year))


PRCP_TAVG_cor <- crosswin(xvar = list( PRCP = monthly_clim$PRCP), 
                          xvar2 = list( TAVG = monthly_clim$TAVG), 
                          cdate = monthly_clim$date, 
                          bdate = growth$date, 
                          stat = 'sum', 
                          stat2 = 'mean', 
                          type = 'absolute', refday = c(15, 06), range = c(15, 0), cinterval = 'month')


VWC_TAVG_cor <- crosswin(xvar = list( VWC = monthly_clim$VWC), 
                         xvar2 = list( TAVG = monthly_clim$TAVG), 
                         cdate = monthly_clim$date, 
                         bdate = growth$date, 
                         stat = 'sum', 
                         stat2 = 'mean', 
                         type = 'absolute', refday = c(15, 06), range = c(15, 0), cinterval = 'month')

save(PRCP_TAVG_cor, VWC_TAVG_cor, file = file.path('data', 'temp_data', 'climwin_correlations.rda'))

sp_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')
iter <- 5
species <- 'PSSP'

for(species in sp_list){ 

  intra <- paste0( 'W.' , species )
  dat <- readRDS( paste0 ( 'data/temp_data/', paste0( species, '_size.RDS')))
  dat <- dat[ dat$year < 1960 , ]
  
  mnyear <- min(dat$year)
  mxyear <- max(dat$year)
  pid <- unique( dat$pid )
  
  dat <- 
    expand.grid( year = mnyear:mxyear, pid = unique(pid)) %>% 
    left_join(dat, by = c('pid', 'year'))  %>% 
    arrange(pid, year) 
  
  my_dat <- 
    dat %>% 
    group_by( pid ) %>% 
    mutate( age = ifelse( year == min(year[!is.na(area)]), 1, age )) %>% 
    dplyr::select(quad, pid, age, year, area, starts_with('W.'))  %>% 
    ungroup() %>% 
    mutate( W.total = W.ARTR + W.HECO + W.POSE + W.PSSP - eval(parse(text = intra)), 
            W.intra = eval(parse(text = intra)))
   
  growth <- 
    my_dat %>% 
    group_by( pid ) %>% 
    arrange(pid, year) %>% 
    mutate( area = log(area)) %>% 
    mutate( area0 = lag(area)) %>% 
    mutate( growth = area - area0 ) %>% 
    filter( !is.na(growth), 
            is.finite(growth), 
            area > -1.3, 
            area0 > -1.3, year > 1927) %>%
    mutate( date = paste0( '15/06/', year)) # date for climwin 
  
  
  growthWin_TAVG <- slidingwin(xvar = list(TAVG = monthly_clim$TAVG_K),
                        cdate = monthly_clim$date,
                        bdate = growth$date,
                        baseline = lm(growth ~ 1 , data = growth),
                        cinterval = "month",
                        range = c(15, 0),
                        type = "absolute", refday = c(15, 06),
                        stat = c("mean"),
                        func = c("lin", "log"))
  
  growthWin_PRCP <- slidingwin(xvar = list(PRCP = monthly_clim$PRCP_pos),
                            cdate = monthly_clim$date,
                            bdate = growth$date,
                            baseline = lm(growth ~ 1 , data = growth),
                            cinterval = "month",
                            range = c(15, 0),
                            type = "absolute", refday = c(15, 06),
                            stat = "sum",
                            func = c("lin", "log"))
  
  growthWin_VWC <- slidingwin(xvar = list(VWC = monthly_clim$VWC),
                                 cdate = monthly_clim$date,
                                 bdate = growth$date,
                                 baseline = lm(growth ~ 1 , data = growth),
                                 cinterval = "month",
                                 range = c(15, 0),
                                 type = "absolute", refday = c(15, 06),
                                 stat = "mean",
                                 func = c("lin", "log"))
  
  
  TAVG_rand <- randwin(repeats = iter, 
                      xvar = list(TAVG = monthly_clim$TAVG_K),
                      cdate = monthly_clim$date,
                      bdate = growth$date,
                      baseline = lm(growth ~ 1 , data = growth),
                      cinterval = "month",
                      range = c(15, 0),
                      type = "absolute", refday = c(15, 06),
                      stat = c("mean"),
                      func = c("lin", "log"))
  
  PRCP_rand <- randwin(repeats = iter, 
                       xvar = list(PRCP = monthly_clim$PRCP_pos),
                       cdate = monthly_clim$date,
                       bdate = growth$date,
                       baseline = lm(growth ~ 1 , data = growth),
                       cinterval = "month",
                       range = c(15, 0),
                       type = "absolute", refday = c(15, 06),
                       stat = c("sum"),
                       func = c("lin", "log"))
  
  VWC_rand <- randwin(repeats = iter, 
                       xvar = list(VWC = monthly_clim$VWC),
                       cdate = monthly_clim$date,
                       bdate = growth$date,
                       baseline = lm(growth ~ 1 , data = growth),
                       cinterval = "month",
                       range = c(15, 0),
                       type = "absolute", refday = c(15, 06),
                       stat = c("mean"),
                       func = c("lin", "log"))
  
  save(growthWin_TAVG, 
       growthWin_PRCP, 
       growthWin_VWC, 
       TAVG_rand, 
       PRCP_rand, 
       VWC_rand, 
       file = paste0( "data/temp_data/", species, "_no_competition.rda"))
}



