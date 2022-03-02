rm(list = ls())

dat <- read_csv('~/Desktop/PSSP_growth_surv.csv')


dat %>% 
  group_by( pid ) %>% 
  filter( year >= min_year, year <= max_year ) %>% 
  mutate( area_censored = as.numeric(area <= 0.3)) %>% 
  mutate( survival_missing = as.numeric(is.na( survives )))


for( i in 1:nrow(data)){ } 
