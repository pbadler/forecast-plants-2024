rm(list = ls()) 

library(tidyverse)

# input ---------------------------------------------------- # 

load('code/figure_scripts/my_plotting_theme.Rdata')

seasonal <- read_csv('data/temp/seasonal_climate.csv')

seasonal <- 
  seasonal %>% 
  filter( YEAR > 1926, YEAR < 2017) 

seasonal$Period[ seasonal$YEAR > 2010 ] <- 'Modern'
seasonal$Period[ seasonal$YEAR < 2011 ] <- 'Historical'

Tavgs <- 
  seasonal %>% 
  filter( var == 'TAVG_avg', Treatment == 'Control') %>% 
  mutate( season = factor( season, levels = c('winter', 'spring', 'summer', 'fall'), ordered = T)) # order for plotting 

historical_range <- 
  Tavgs %>% 
  filter( Period == 'Historical') %>% 
  group_by( season ) %>% 
  summarise(lower = quantile(val, 0.05), upper = quantile(val, 0.95), lt_mean = mean(val, na.rm = T))

gg_lt_Tavg <- 
  Tavgs %>% 
  left_join(historical_range, by= 'season') %>% 
  ggplot( aes( x = YEAR, y = val )) + 
  geom_line() + 
  geom_hline( aes( x = YEAR, yintercept = lower), linetype = 2, alpha = 0.5) +
  geom_hline( aes( x = YEAR, yintercept = upper ), linetype = 2, alpha = 0.5) +
  geom_vline(aes( xintercept = YEAR[ which.max(YEAR[Period == 'Historical'])] ), color = 'blue', linetype = 2) + 
  facet_wrap(~season, scale = 'free_y') + 
  ylab( 'Average daily average temperature') +
  theme_bw() +
  theme( panel.grid =  element_blank()) 

gg_lt_Tavg

Tavgs %>% 
  left_join(historical_range, by= 'season') %>% 
  filter(YEAR > 2010) %>% 
  ggplot( aes( x = YEAR, y = val )) + 
  geom_line() + 
  geom_hline( aes( x = YEAR, yintercept = lower), linetype = 2, alpha = 0.5) +
  geom_hline( aes( x = YEAR, yintercept = upper ), linetype = 2, alpha = 0.5) +
  #geom_vline(aes( xintercept = YEAR[ which.max(YEAR[Period == 'Historical'])] ), color = 'blue', linetype = 2) + 
  facet_wrap(~season, scale = 'free_y') + 
  ylab( 'Average daily average temperature') +
  theme_bw() +
  theme( panel.grid =  element_blank()) 

ggsave(gg_lt_Tavg, filename = 'figures/lt_TAVG.png', height = 6, width = 6)
