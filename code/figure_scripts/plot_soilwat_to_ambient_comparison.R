rm(list = ls())
library(tidyverse)
library(lubridate)
library(lme4)
library(sheepweather)


# output --------------------------------------------------------------- # 
fig_file <- 'figures/soilwat_v_ambient_obs.png'

# input ---------------------------------------------------------------- # 
load('code/figure_scripts/my_plotting_theme.Rdata')

sw <- sheepweather::usses_soilwat
daily_sm <- read_csv('data/temp/daily_sm.csv')
seasons <- read_csv('data/season_table.csv')

# ----------------------------------------------------------------------- # 

sw_avg <- 
  sw %>% 
  mutate( date = date(parse_date_time(paste( Year, Day, sep = '-') , orders = '%Y-%j'))) %>% 
  rename( 'year' = Year) %>% 
  mutate( month = month( date )) %>% 
  gather(layer, VWC, Lyr_1:Lyr_4) %>%
  filter(layer %in% c('Lyr_1')) %>%
  group_by(date) %>%
  summarise(swVWC = mean(VWC) * 100) %>% 
  ungroup() 


compare_sw_df <- 
  sw_avg %>%   
  filter( date > '2012-01-01', date < '2017-01-01') %>% 
  left_join(
    daily_sm %>%
      filter( Treatment == 'Control') %>% 
      group_by( date) %>%
      summarise( obsVWC = mean(raw))
  ) %>% 
  rename('Observed' = obsVWC) %>% 
  rename('SOILWAT' = swVWC) %>% 
  gather( type, value, SOILWAT, Observed) 

gg_compare <- compare_sw_df %>% 
  mutate( year = year(date), DOY = yday(date)) %>% 
  ggplot(aes( x = DOY, y = value, color = type )) + 
  geom_line() + 
  facet_wrap(~year, nrow = 4) + 
  facet_wrap(~ year, nrow = 5) + 
  scale_linetype_manual(values = c(1, 2)) +
  scale_alpha_manual(values = c(1, 0.7)) +
  my_theme + 
  ylab('Volumetric soil moisture (ml/ml)') + 
  xlab('Day of year')


ggsave(gg_compare, filename = fig_file, width = 7, height = 7)

