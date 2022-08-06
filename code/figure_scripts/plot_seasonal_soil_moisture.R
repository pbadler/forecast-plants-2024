rm(list = ls()) 

library(tidyverse)
library(zoo)

# input ---------------------------------------------------- # 

load('code/figure_scripts/my_plotting_theme.Rdata')

seasonal_VWC <- read.csv('data/temp/seasonal_VWC.csv') 
annual_VWC <- read.csv('data/temp/annual_VWC.csv')
monthly_avg <- read.csv('data/temp/monthly_avg.csv')

seasonal_VWC <- 
  seasonal_VWC %>% 
  filter( year > 1926, year < 2017)

annual_VWC <- 
  annual_VWC %>% 
  filter( year > 1926, year < 2017)

monthly_avg <- 
  monthly_avg %>% 
  filter( year > 1926, year < 2017)

# ---------- monthly soil moisture --------------------------------------------------# 
#

gg_monthly <- 
  ggplot( monthly_avg, aes( x = factor( month) , y = avg_VWC) ) +
  geom_boxplot()  +
  ylim( 0, 35) +
  ylab( 'Soil water content (ml/ml)') +
  xlab( 'Month') +
  my_theme

ggsave(gg_monthly, filename = 'figures/longterm_monthly_soilwat.png', height = 5, width = 5)

# Seasonal VWC ---------------------------------------------------------- # 

seasonal_VWC$Period[ seasonal_VWC$year > 2010 ] <- 'Modern'
seasonal_VWC$Period[ seasonal_VWC$year < 2011 ] <- 'Historical'

historical_avgs <-
  seasonal_VWC %>%
  filter( Period == 'Historical')%>%
  group_by( Period, season ) %>%
  summarise(lower = quantile(avg, 0.05), upper = quantile(avg, 0.95), avg = mean(avg) ) %>%
  spread( Period, avg )

modern <-
  seasonal_VWC %>%
  filter( Period == 'Modern') %>%
  group_by( Period, season, year )

modern <- merge( modern, historical_avgs )

modern$season <- factor(modern$season, levels = c('winter', 'spring', 'summer', 'fall'), ordered = T)

seasonal_VWC$season <- factor(seasonal_VWC$season, levels = c('winter', 'spring', 'summer', 'fall'), ordered = T)

gg_lt_seasonal_historic <- 
  seasonal_VWC %>% 
  filter( season != 'winter') %>% 
  filter( Treatment == 'Control' & year < 1960) %>% 
  ggplot( aes(x = year, y = avg, color = Treatment)) + 
  geom_line() + 
  facet_wrap(~season, nrow = 4) + 
  my_theme + 
  scale_color_manual(values = my_colors[2]) + 
  ylab( 'Avg. seasonal soil moisture (%)') + 
  ylim(0, 28)  + 
  guides(color = 'none')

gg_lt_seasonal_modern <- 
  seasonal_VWC %>% 
  filter( season != 'winter') %>% 
  filter( year > 2006) %>% 
  filter( Treatment == 'Control' | year > 2011 ) %>% 
  ggplot( aes(x = year, y = avg, color = Treatment)) + 
  geom_line() + 
  facet_wrap(~season, nrow = 4) +
  scale_color_manual(values = my_colors[2:4]) + 
  my_theme + 
  theme(axis.title.y = element_blank()) + 
  ylim(0,28) 

library(gridExtra)
gg_lt_sm <- grid.arrange(gg_lt_seasonal_historic, gg_lt_seasonal_modern, nrow = 1 )

ggsave(gg_lt_sm, filename = 'figures/lt_seasonal_sm.png', width = 7, height =10 )

# make plot of seasonal moisture during the experiment  ------------------------------------------- #

gg_soil_moist <- 
  modern %>% 
  filter( year > 2011) %>% 
  ggplot( aes( x = year, y = avg, color = Treatment ) ) +
  geom_line() +
  geom_point() +
  facet_grid(season ~ . ) +
  geom_hline( aes( x = year, yintercept  = lower), linetype = 2, alpha = 0.5) +
  geom_hline( aes( x = year, yintercept = upper ), linetype = 2, alpha = 0.5) +
  scale_color_manual(values = my_colors[2:4]) +
  ylab( 'Average seasonal soil moisture content (ml/ml)') +
  theme_bw() +
  theme( panel.grid =  element_blank())


ggsave(gg_soil_moist, 
       filename = 'figures/modern_soil_moisture_comparison.png', 
       width = 6, height = 6)

# Just focus on spring 
gg_spring <- 
  modern %>% 
  filter(season == 'spring' &  year > 2011) %>% 
  ggplot(aes( x = year, y = avg, color = Treatment ) ) +
    geom_line() +
    geom_point() +
    geom_hline( aes( x = year, yintercept  = lower), linetype = 2, alpha = 0.5)  +
    geom_hline( aes( x = year, yintercept = upper ), linetype = 2, alpha = 0.5)  +
    scale_color_manual(values = my_colors[2:4])  +
    ylab( 'Average seasonal soil moisture content (ml/ml)')  +
    theme_bw() +
    theme( panel.grid =  element_blank())

ggsave(gg_spring, 
      filename = 'figures/modern_soil_moisture_comparison_spring.png', 
       width = 5, height = 3)




