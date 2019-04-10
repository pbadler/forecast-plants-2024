##### Spring 2015 soil moisture
##### 

rm(list = ls())
library(tidyverse)

load('code/figure_scripts/my_plotting_theme.Rdata')

soil_moist_dir <- 'data'
q_info <-read.csv(paste0(soil_moist_dir, '/quad_info.csv') )

calibration <- read.csv(paste0(soil_moist_dir, '/2015-05-07_soil_probe_calibration.csv'))
calibration2 <- read.csv(paste0(soil_moist_dir, '/2015-04-29_soil_probe_calibration.csv'), skip = 8)

p1 <- read.csv(paste0( soil_moist_dir, '/2012-06-06_spot_measurements.csv'), skip = 3)

p1$date <- '2012-06-06'
p1$Plot <- gsub( p1$Plot, pattern = '-', replacement = '_')
p1$rep <- c(1:2)
p1 <- p1 %>% rename( plot = Plot )

p2 <- read.csv(paste0(soil_moist_dir, '/2015-04-29_spot_measurements.csv'), skip = 2)
p3 <- read.csv(paste0(soil_moist_dir, '/2015-05-07_spot_measurements.csv'))
p4 <- read.csv(paste0(soil_moist_dir, '/2016-05-10_spot_measurements.csv'))
p5 <- read.csv(paste0(soil_moist_dir, '/2015-06-09_spot_measurements.csv'))

p2$date <- '2015-04-29'

df <- rbind( p2, p3, p4, p5)

df <- df %>% gather( key = rep, PCT, E1:W3 )

df <- rbind( p1, df )

q_info$plot <- gsub( q_info$QuadName, pattern = 'X', replacement = '')

df <- merge (df, q_info , by = 'plot')

df$date <- as.POSIXct(df$date, tz = 'MST')
df <- df %>% rename(VWC = PCT)

soil_density <- read.csv(paste0(soil_moist_dir, '/exclosure_soil_samples.csv'))

soil_density <- soil_density %>% 
  mutate( layer = ifelse( depth > 15, 'deep', 'shallow') ) %>% 
  mutate( wet_soil = wet_weight - bag_weight, soil = dry_weight - bag_weight ) %>% 
  mutate( water_weight = wet_soil - soil, soil_density = soil/soil_volume_cm3, GWC = water_weight/soil) %>% 
  group_by(layer) %>% 
  summarise( avg_soil_density = mean(soil_density))

soil_d_value <- soil_density$avg_soil_density[ soil_density$layer == 'shallow']

# run calibration ---------------------------------------------------------------------------------------------------- 
calib <- 
  calibration %>% 
  mutate( wet_weight = (wet_weight- bag_weight) , 
          dry_weight = (dry_weight - bag_weight), 
          water_weight = wet_weight - dry_weight, 
          GWC = 100*water_weight/dry_weight, 
          VWC_true = GWC*soil_d_value ) %>% 
  gather( rep, VWC, T1:T4 ) %>% 
    group_by( date, plot ) %>% 
    summarise( VWC_true = mean(VWC_true), VWC = mean(VWC))

m1 <- lm( data = calib, VWC_true ~ poly(VWC, 1))
summary(m1)


ggplot( data = calib, aes( x = VWC, y = VWC_true)) + geom_point() + geom_smooth(method = 'lm', formula = 'y ~ poly(x, 1)') + xlim( 0, 35) + ylim(0,35)

predict(m1, newdata = data.frame(VWC = 0))

df$VWC_corrected <- predict(m1, newdata = data.frame( VWC = df$VWC) )

# spot measurements --------------------------------------------------------------------------------------------------

plot_by_date <- ggplot( df , aes ( x = factor(date), y = VWC_corrected , color = Treatment )) + 
  geom_boxplot() + 
  ylab( 'Volumetric water content (ml/ml)') + 
  geom_point(position = position_dodge(width = 0.8), alpha = 0.5) + 
  scale_color_manual(values = my_colors[2:4]) + 
  my_theme + labs( color= '')

plot_by_treat <- ggplot( df, aes( x = Treatment, y = VWC_corrected, color = Treatment)) + 
  geom_boxplot() + 
  ylab( 'Soil volumetric water content (ml/ml)') + 
  geom_point( position = position_dodge(width = 0.8), alpha = 0.5) + 
  scale_color_manual(values = my_colors[2:4]) + 
  my_theme + labs( color= '')

plot_by_group <- 
  ggplot( df, aes( x = Treatment , y  = VWC, color = Treatment )) + 
  geom_boxplot( ) + 
  facet_wrap(~PrecipGroup)

plot_by_group

plot_by_date + ggtitle( 'Spot measurements')

# analysis -------------------------------------------------------------------------------------------------------------

spot_form <- formula(VWC ~ (1|date) + (1|PrecipGroup) + (1|plot) + Treatment)

spot_m <- lmer( data = df , spot_form)
spot_m_2016 <- lmer(data = subset(df, date > '2016-01-01'), VWC ~ (1|plot) + Treatment)

summary(spot_m)
summary(spot_m_2016)

drop1(spot_m)
anova(spot_m)

aggregate(data = df, VWC ~ Treatment, FUN = 'mean')

4.19/8.16
11.42/8.16

# make table for stats output 
library(xtable)
library(lsmeans)
statsOutput <- 'manuscript/spot_measurements.tex'

test <- lsmeans(spot_m , "Treatment" )

#xtlsm_spot <- xtable(test, caption = '')
#print(xtlsm_spot, )

texreg(spot_m,caption="Estimated parameters from a mixed effects model fit to the spring soil moisture data. Parameter with s.e. are shown. Intercept gives soil moisture in ambient controls.",
        caption.above=TRUE,file=statsOutput, label = 'table:spotVWC')


aggregate( data = df, VWC ~ Treatment, FUN = 'mean')


# print figures --------------------------------------------------------------------------------------------------------

png( 'figures/VWC_spot_measurements.png', height = 5 , width = 7, res = 300, units = 'in' ) 

print( plot_by_treat ) 

dev.off()

# Save spot measurements -----------------------------------------------------------------------------------------------

#saveRDS( df, 'data/processed_data/spring_spot_measurements.RDS')
