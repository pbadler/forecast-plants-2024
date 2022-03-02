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
# ClimWin Analysis to find best window for survival 

monthly_clim <- 
  monthly_clim %>% 
  group_by( month) %>% 
  mutate( TAVG_K = TAVG + 273.15) %>% 
  mutate( PRCP_pos = PRCP + 0.00001) %>%  # ensure PRCP > 0 for log trans.
  mutate( date = paste0( '15/', str_pad( month, 2, pad = 0), '/', year))

sp_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')
iter <- 1

species <- 'ARTR'

for(species in sp_list){ 

  intra <- paste0( 'W.' , species )
  dat <- readRDS( paste0 ( 'data/temp_data/', paste0( species, '_survival.RDS')))
  dat <- dat[ dat$year < 1960 , ]
  
  dat <- 
    dat %>% 
    mutate( pid = paste( quad, trackID, sep = "_")) %>%
    group_by( pid ) %>% 
    arrange( year) %>% 
    mutate( year = year + 1)
  
  mnyear <- min(dat$year)
  mxyear <- max(dat$year)
  pid <- unique( dat$pid )
  
  dat <- 
    expand.grid( year = mnyear:mxyear, pid = unique(pid)) %>% 
    left_join(dat, by = c('pid', 'year'))  %>% 
    arrange(pid, year) 
  
  my_dat <- 
    dat %>% 
    ungroup() %>% 
    mutate( W.total = W.ARTR + W.HECO + W.POSE + W.PSSP - eval(parse(text = intra)), 
            W.intra = eval(parse(text = intra)))
   
  survives <- 
    my_dat %>% 
    group_by( pid ) %>% 
    arrange(pid, year) %>% 
    filter( !is.na(logarea), 
            !is.na(survives),
            is.finite(logarea), 
            year > 1927) %>%
    mutate( date = paste0( '15/06/', year)) # date for climwin 
  
  survWin_TAVG <- slidingwin(xvar = list(TAVG = monthly_clim$TAVG_K),
                        cdate = monthly_clim$date,
                        bdate = survives$date,
                        baseline = glm(survives ~ 1 + logarea , data = survives, family = 'binomial'),
                        cinterval = "month",
                        range = c(15, 0),
                        type = "absolute", refday = c(15, 06),
                        stat = c("mean"),
                        func = c("lin", "log"))
  
  survWin_PRCP <- slidingwin(xvar = list(PRCP = monthly_clim$PRCP_pos),
                            cdate = monthly_clim$date,
                            bdate = survives$date,
                            baseline = glm(survives ~ 1 + logarea , data = survives, family = 'binomial'),
                            cinterval = "month",
                            range = c(15, 0),
                            type = "absolute", refday = c(15, 06),
                            stat = "sum",
                            func = c("lin", "log"))
  
  survWin_VWC <- slidingwin(xvar = list(VWC = monthly_clim$VWC),
                              cdate = monthly_clim$date,
                              bdate = survives$date,
                              baseline = glm(survives ~ 1 + logarea , data = survives, family = 'binomial'),
                              cinterval = "month",
                              range = c(15, 0),
                              type = "absolute", refday = c(15, 06),
                              stat = "mean",
                              func = c("lin", "log"))
  
  
  survTAVG_rand <- randwin(repeats = iter, 
                      xvar = list(TAVG = monthly_clim$TAVG_K),
                      cdate = monthly_clim$date,
                      bdate = survives$date,
                      baseline = glm(survives ~ 1 + logarea , data = survives, family = 'binomial'),
                      cinterval = "month",
                      range = c(15, 0),
                      type = "absolute", refday = c(15, 06),
                      stat = c("mean"),
                      func = c("lin", "log"))
  
  survPRCP_rand <- randwin(repeats = iter, 
                       xvar = list(PRCP = monthly_clim$PRCP_pos),
                       cdate = monthly_clim$date,
                       bdate = survives$date,
                       baseline = glm(survives ~ 1 + logarea , data = survives, family = 'binomial'),
                       cinterval = "month",
                       range = c(15, 0),
                       type = "absolute", refday = c(15, 06),
                       stat = c("sum"),
                       func = c("lin", "log"))
  
  survVWC_rand <- randwin(repeats = iter, 
                       xvar = list(VWC = monthly_clim$VWC),
                       cdate = monthly_clim$date,
                       bdate = survives$date,
                      baseline = glm(survives ~ 1 + logarea , data = survives, family = 'binomial'),
                      cinterval = "month",
                       range = c(15, 0),
                       type = "absolute", refday = c(15, 06),
                       stat = c("mean"),
                       func = c("lin", "log"))
  
  save(survWin_TAVG, 
       survWin_PRCP, 
       survWin_VWC, 
       survTAVG_rand, 
       survPRCP_rand, 
       survVWC_rand, 
       file = paste0( "data/temp_data/", species, "survival_no_competition.rda"))
}



