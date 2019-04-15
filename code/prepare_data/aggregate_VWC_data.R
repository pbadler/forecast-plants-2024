#######################################################################################
#
# Setup seasonal SoilWAT variables for demographic rate models 
#
#######################################################################################

rm(list = ls()) 

library(tidyverse)
library(zoo)

# input ---------------------------------------------------- # 

load('code/figure_scripts/my_plotting_theme.Rdata')

df <- readRDS('data/temp_data/daily_swVWC_treatments.RDS')
  # comes from the soilMoistureTreatmentEffects script 

# output ---------------------------------------------------- # 

seasonal_outfile <- 'data/temp_data/seasonal_VWC.RDS'
#quarterly_outfile <- 'data/temp_data/quarterly_VWC.RDS'
annual_outfile <- 'data/temp_data/annual_VWC.RDS'

# ---------- annual soil moisture -------------------------------------------------#
annual_VWC <- 
  df %>% 
  group_by( Period, year, Treatment ) %>%
  summarise (avg_VWC = mean(VWC_raw, na.rm = TRUE))

# ---------- seasonal soil moisture  -----------------------------------------------#
seasonal_VWC <- 
  df %>% 
  mutate(year = ifelse(month == 12 , year + 1, year  )) %>% # account for December 
  group_by(Period, year, season, Treatment) %>% 
  summarise( avg = mean(VWC_raw, na.rm = TRUE) )

# ------------ aggregate VWC with Treatment effects by quarter ---------------#

quarterly_VWC <-
  df %>% 
  group_by( Period, Treatment, year, quarter ) %>%                       
  summarise( avg = mean(VWC_raw, na.rm = TRUE)) 

# -------- output -----------------------------------------------------------------------------#

saveRDS( seasonal_VWC, seasonal_outfile) 
saveRDS( annual_VWC, annual_outfile) 

# ---------- monthly soil moisture --------------------------------------------------# 
# 
# monthly_avgs <- 
#   df %>% 
#   filter( Treatment == 'Control') %>% 
#   group_by( month, year ) %>%
#   summarise (avg_VWC = mean(VWC_raw, na.rm = TRUE) )  
# 
# #
# 
# png( 'figures/longterm_monthly_soilwat.png', height = 5, width = 5, res = 300, units = 'in')
# print( 
#   ggplot( monthly_avgs, aes( x = factor( month) , y = avg_VWC) ) + 
#     geom_boxplot()  + 
#     ylim( 0, 35) + 
#     ylab( 'Soil water content (ml/ml)') + 
#     xlab( 'Month') + 
#     my_theme 
# )
# dev.off()
# 
# 
# seasonal_VWC$Period[ seasonal_VWC$year > 2010 ] <- 'Modern'
# seasonal_VWC$Period[ seasonal_VWC$year < 2011 ] <- 'Historical'
# 
# historical_avgs <- 
#   seasonal_VWC %>% 
#   filter( Period == 'Historical')%>%
#   group_by( Period, season ) %>% 
#   summarise(lower = quantile(avg, 0.05), upper = quantile(avg, 0.95), avg = mean(avg) ) %>% 
#   spread( Period, avg )
# 
# modern <- 
#   seasonal_VWC %>% 
#   filter( Period == 'Modern') %>% 
#   group_by( Period, season, year ) 
# 
# modern <- merge( modern, historical_avgs )
# 
# modern$season <- factor(modern$season, levels = c('winter', 'spring', 'summer', 'fall'), ordered = T)
# 
# # make plot of seasonal moisture during the experiment  ------------------------------------------- # 
# 
# png('figures/modern_soil_moisture_comparison.png', width = 6, height = 6, res = 300, units  = 'in')
# 
# print( 
# ggplot( modern %>% filter( year > 2011), aes( x = year, y = avg, color = Treatment ) ) + 
#   geom_line() + 
#   geom_point() + 
#   facet_grid(season ~ . ) + 
#   geom_hline( aes( x = year, yintercept  = lower), linetype = 2, alpha = 0.5) + 
#   geom_hline( aes( x = year, yintercept = upper ), linetype = 2, alpha = 0.5) + 
#   scale_color_manual(values = my_colors[2:4]) + 
#   ylab( 'Average seasonal soil moisture content (ml/ml)') + 
#   theme_bw() + 
#   theme( panel.grid =  element_blank()) 
# )
# 
# dev.off()
# 
# png('figures/modern_soil_moisture_comparison_spring.png', width = 5, height = 3, res = 300, units  = 'in')
# 
# print( 
#   ggplot( data = subset(modern, modern$season == 'spring' & modern$year > 2011), aes( x = year, y = avg, color = Treatment ) ) + 
#     geom_line() + 
#     geom_point() + 
#     geom_hline( aes( x = year, yintercept  = lower), linetype = 2, alpha = 0.5)  + 
#     geom_hline( aes( x = year, yintercept = upper ), linetype = 2, alpha = 0.5)  + 
#     scale_color_manual(values = my_colors[2:4])  + 
#     ylab( 'Average seasonal soil moisture content (ml/ml)')  + 
#     theme_bw() + 
#     theme( panel.grid =  element_blank()) 
# )
# 
# dev.off()
# 
#   