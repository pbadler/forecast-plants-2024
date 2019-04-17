rm(list = ls()) 

library(tidyverse)
library(zoo)

# input ---------------------------------------------------- # 

load('code/figure_scripts/my_plotting_theme.Rdata')

seasonal_VWC <- readRDS('data/temp_data/seasonal_VWC.RDS') 
annual_VWC <- readRDS('data/temp_data/annual_VWC.RDS')
monthly_avg <- readRDS('data/temp_data/monthly_avg.RDS')


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
       filename = 'figures/modern_soil_moisture_comparison_spring.png', 
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
