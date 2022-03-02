rm(list = ls())
library(tidyverse)

species <- 'PSSP'

intra <- paste0( 'W.' , species )
dat <- readRDS( paste0 ( 'data/temp_data/', paste0( species, '_size.RDS')))
dat <- dat[ dat$year < 1960 , ]

sdat <- readRDS( paste0 ( 'data/temp_data/', paste0( species, '_survival.RDS')))
sdat <- sdat[ sdat$year < 1960 , ]

mnyear <- min(c(sdat$year, dat$year))
mxyear <- max(c(sdat$year, dat$year))

dat <- 
  dat %>% 
  mutate( pid = paste(quad, trackID, sep = '-'))

sdat <- 
  sdat %>% 
  mutate( pid = paste( quad, trackID, sep = '-'))
          
pid <- unique(c( dat$pid, sdat$pid ))

dat <- 
  expand.grid( year = mnyear:mxyear, pid = pid) %>% 
  left_join(dat, by = c('pid', 'year'))  %>% 
  arrange(pid, year) 

sdat <- 
  expand.grid( year = mnyear:mxyear, pid = pid) %>% 
  left_join(sdat, by = c('pid', 'year')) %>% 
  arrange( pid , year )

dat <- 
  dat %>% 
  select( year, pid, area, starts_with('W')) %>% 
  left_join(sdat %>% select( year,  pid, survives), by = c('pid', 'year')) %>% 
  arrange( pid, year ) 

dat <- 
  dat %>% 
  group_by( pid ) %>% 
  arrange( pid, year ) %>%  
  mutate( max_year = max(year[!is.na(area)]), min_year = min(year[!is.na(area)])) %>% 
  mutate( fill_in = ifelse( is.na(survives) & year < max_year & year > min_year, 1, survives))

dat %>%
  mutate( survives = fill_in ) %>% 
  write_csv(path = paste0('~/Desktop/', species, '_growth_surv.csv'))




