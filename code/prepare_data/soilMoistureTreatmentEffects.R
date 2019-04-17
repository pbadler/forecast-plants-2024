# calculate treatment effects compared to control soil moisture

rm(list = ls())

library(tidyverse)
library(lubridate)
library(lme4)
library(lsmeans)
library(texreg)
library(sheepweather)

# input ---------------------------------------------------- #

load('code/figure_scripts/my_plotting_theme.Rdata')

swVWC <- usses_soilwat  # soilwat output
spotVWC <- readRDS('data/temp_data/spotVWC.RDS') # from aggregate spot script 

myVWC <- readRDS('data/temp_data/decagon_data_with_station_data.RDS')   # cleaned decagon data
daily_clim <- readRDS('data/temp_data/daily_station_dat_rainfall.RDS')  # climate station rainfall

seasons <- read.csv('data/season_table.csv')
quads <- read_csv('data/quad_info.csv')

# output ---------------------------------------------------- #

daily_swVWC_treatment_outfile <- 'data/temp_data/daily_swVWC_treatments.RDS'
fig1_file <- 'figures/avg_daily_soil_moisture.png'
fig2_file = 'figures/modeled_soilwat_soil_moisture_example.pdf'

# --------------------------------------------------------------------------------------#
myVWC <-
  myVWC %>%
  mutate(v = ifelse(measure == 'VWC', v * 100, v)) %>%  # convert to percent
  mutate( year = year(date))

# Steps for treatment effects standardization:
# 1. aggregate soil moisture by PrecipGroup, Treatment, depth and date
# 2. then standardize soil moisture measurements within PrecipGroup and depth
# 3. then subtract standardized Drought and Irrigation from Control to get
#       standardized treatment effects.
# 4. output

myVWC <-
  myVWC %>%
  filter(measure == 'VWC') %>%
  filter(!is.na(rainfall)) %>%
  filter(!(plot == 'X16')) # drop plot 16

myVWC <-
  myVWC %>%
  ungroup() %>%
  mutate( unique_position = paste0(plot, '.', position)) %>%
  mutate( month = month(date)) %>% 
  left_join(seasons, by = 'month')

VWC_weights <-
  myVWC %>%
  group_by(year, PrecipGroup, date) %>%
  summarise(weight = n_distinct(unique_position))

df_soil_moist <-
  myVWC %>%
  group_by(year, PrecipGroup, Treatment, season, date, rainfall) %>%
  summarise(avg_VWC = mean(v, na.rm = TRUE)) %>%
  group_by(PrecipGroup) %>%
  mutate(avg_VWC = scale(avg_VWC, mean(avg_VWC[Treatment == 'Control'], na.rm = T), sd(avg_VWC[Treatment == 'Control'], na.rm = T))) %>%  # scale within Precip Group and Depth
  spread(Treatment, avg_VWC) %>%
  mutate(Drought = Drought - Control, Irrigation = Irrigation - Control) %>%
  arrange(PrecipGroup, date)

df_soil_moist <- 
  df_soil_moist %>% 
  left_join(VWC_weights, by = c('year', 'PrecipGroup', 'date'))

# -------- add in the spot measurements ------------------------------------------------#
df_soil_moist$type <- 'logger'
VWC_df <- df_soil_moist

spotVWC$type <- 'spot'
spotVWC$year <- year(spotVWC$date)

VWC_df <- rbind(VWC_df, spotVWC[,-1])

#
#  Model the standardized treatment differences.
#  Treatment effects vary with season and whether it's a rainy period.
#  See definition of rainy periods in the datadrivers data preparation scripts.
#

# fit models

VWC_test <- VWC_df %>%
  gather(Treatment, VWC, Drought , Irrigation)

mTreatment <-
  lm(data = VWC_test,
     VWC ~ Treatment * rainfall * season,
     weights = VWC_test$weight)
# select models

mTreatment <-
  MASS::stepAIC(mTreatment,
          scope = list(upper = ~ . , lower = ~ 1),
          trace = T)

coef(mTreatment)
summary(mTreatment)

mTreatment.mer <-
  lmer(update(formula(mTreatment) , . ~ . + (1 | date) + (1 | PrecipGroup)),
       data = VWC_test,
       weights = VWC_test$weight)


summary(mTreatment.mer)

test <- lsmeans(mTreatment.mer,  ~ Treatment + season + rainfall)

tab <- summary(test)
tab <-  data.frame(tab)

tab <-
  tab %>%
  dplyr::select(season, rainfall,  Treatment, lsmean, SE, asymp.LCL, asymp.UCL) %>%
  arrange(season, rainfall, Treatment)

statsOutput <- 'tables/soil_moisture_model.tex'

texreg(
  mTreatment,
  caption = "Treatment effects on soil moisture. Intercept refers to drought effects in fall not rainy conditions.  Model fit to the continuously logged soil moisture data as well as the spot measurements collected from all plots in the spring.",
  caption.above = TRUE,
  file = statsOutput,
  label = 'table:soil_moisture_model'
)


# data frame to view predictions
pred_df <-
  expand.grid(
    rainfall = unique(VWC_test$rainfall),
    Treatment = unique(VWC_test$Treatment),
    season = unique(VWC_test$season)
  )

pred_df$mu <- predict(mTreatment, newdata = pred_df,  re.form = NA)

#

daily_VWC <-
  myVWC %>%
  ungroup() %>%
  filter(!is.na(v)) %>%
  dplyr::select(date, year, season, rainfall, Treatment, v) %>%
  group_by(Treatment, date, year, season, rainfall) %>%
  summarise(v = mean(v , na.rm = T)) %>%
  ungroup()

# table of seasonal differences on raw scale
daily_VWC %>%
  group_by(Treatment, season, rainfall) %>%
  summarise(avg = mean(v) , sd = sd(v)) %>%
  arrange(season, rainfall, Treatment) %>%
  gather(stat, v, avg:sd) %>%
  spread(Treatment , v) %>%
  mutate (dp = (Drought - Control) / Control,
          ip = (Irrigation - Control) / Control) %>% filter(stat == 'avg')

daily_control <-
  daily_VWC %>%
  filter(Treatment == 'Control') %>%
  mutate(v = as.numeric(scale(v))) %>%
  spread(Treatment, v)

daily_VWC2 <-
  daily_VWC %>%
  filter(Treatment == 'Control') %>%
  spread(Treatment, v)

pred_df <-
  daily_VWC %>%
  filter(Treatment != 'Control') %>%
  dplyr::select(-v)

pred_df$predicted <-
  predict(mTreatment.mer, pred_df, re.form = NA)

pred_df <-
  pred_df %>%
  spread(Treatment, predicted)

pred_df <-
  merge(pred_df, daily_control)  %>%
  mutate(Drought = Drought + Control, Irrigation = Irrigation + Control)  %>%
  gather(Treatment, predicted, Drought:Control)

observed_df  <-
  daily_VWC %>%
  rename(observed = v)

plot_df <- 
  pred_df %>% 
  left_join(observed_df, by = c('date', 'year', 'season', 'rainfall', 'Treatment')) 

plot_df <-
  plot_df %>%
  mutate(back_scaled_pred = predicted * sd(observed[Treatment == 'Control']) + mean(observed[Treatment == 'Control']))

plot_df <-
  plot_df %>%
  dplyr::select(-predicted) %>%
  distinct() %>%
  rename(predicted = back_scaled_pred)  %>%
  gather(type, VWC, predicted , observed)

everyday <-
  expand.grid(
    date = seq(ymd('2012-01-01'), ymd('2016-12-30'), 1),
    Treatment = c('Control', 'Drought', 'Irrigation'),
    type = c('predicted', 'observed')
  )

plot_df <- 
  everyday %>% 
  left_join(plot_df, by = c('date', 'Treatment', 'type'))

plot_df <-
  plot_df %>% mutate(julian_date = yday(date),
                     year = year(date))


ggplot(plot_df,
       aes(
         x = julian_date,
         y = VWC,
         color = Treatment,
         linetype = type,
         alpha = type
       )) +
  geom_line() +
  facet_grid(year ~ .) +
  scale_color_manual(values = my_colors[2:4]) +
  scale_linetype_manual(values = c(2, 1)) +
  scale_alpha_manual(values = c(1, 0.7)) +
  my_theme

png(
  fig1_file,
  width = 5,
  height = 6,
  res = 300,
  units = 'in'
)

print(
  ggplot(
    subset(plot_df, type != 'predicted'),
    aes(x = julian_date, y = VWC, color = Treatment)
  ) +
    geom_line(alpha = 0.8) +
    facet_grid(year ~ .) +
    scale_color_manual(values = my_colors[2:4]) +
    scale_alpha_manual(values = c(0.7)) +
    labs(color = '') +
    ylab('Soil volumetric water content (ml/ml)') +
    xlab('Day of year') +
    my_theme
)
dev.off()

# predict the treatment effects from the scaled model VWC ----------------------------#
# ---process dates----------------------------------------------------------------------#

swVWC <- 
  swVWC %>% 
  mutate( date = date(parse_date_time(paste( Year, DOY, sep = '-') , orders = '%Y-%j'))) %>% 
  rename( 'year' = Year) %>% 
  mutate( month = month( date )) # use lubridate function 

lyr3 <-
  swVWC %>% 
  group_by(month, year) %>% 
  summarise(avg = mean(Lyr_3))

lyr3 <- lyr3 %>% spread(month, avg)
lyr3 <- lyr3[complete.cases(lyr3),]
mydata <- lyr3[, 2:13]
mydata <- scale(mydata)

# set-up aggregate seasonal variables for model ----------------------------------------#

swVWC <-
  swVWC %>%
  gather(layer, VWC, Lyr_1:Lyr_6) %>%
  filter(layer %in% c('Lyr_1', 'Lyr_2', 'Lyr_3', 'Lyr_4'))

swVWC <- 
  swVWC %>%
  group_by(date) %>%
  summarise(modelVWC = mean(VWC) * 100)

swVWC$month <- month(swVWC$date)
swVWC$year <- year(swVWC$date)

daily_clim <-
  daily_clim %>%
  ungroup() %>%
  dplyr::select(-year)

swVWC <- left_join(swVWC, seasons, by = 'month')
swVWC <- left_join(swVWC, daily_clim, by = 'date')

soilWAT <-
  swVWC %>% 
  dplyr::select(date, modelVWC, year, season, rainfall)

soilWAT$SW_predicted <- soilWAT$modelVWC
daily_VWC2$observed <- daily_VWC2$Control

obs_predicted <-
  merge(daily_VWC2,
        soilWAT,
        by = c('date', 'year', 'season', 'rainfall'))

ggplot(obs_predicted, aes(x = SW_predicted, y = observed)) + geom_point()

obs_predicted <-
  obs_predicted %>% 
  gather(type, val, observed, SW_predicted)

ggplot(obs_predicted, aes(x = date, y = val, color = type)) + geom_line()

swVWC$Control <- scale(swVWC$modelVWC) # standardize control SWC
Control_mean <- mean(swVWC$modelVWC)
Control_sd <- sd(swVWC$modelVWC)

swVWC2 <- swVWC
swVWC$Treatment <-  'Drought'
swVWC2$Treatment <- 'Irrigation'

swVWC <- rbind(swVWC, swVWC2)

swVWC$predicted <- predict(mTreatment,  swVWC, re.form = NA)

swVWC <- 
  swVWC %>% 
  spread(Treatment, predicted)

# unscale the Control Drought and Irrigation VWC --------------------------------------------------------  #

swVWC <-
  swVWC %>%
  mutate(Drought = Control + Drought,
         Irrigation = Control + Irrigation) %>%
  gather(Treatment, VWC, Control:Irrigation)

swVWC$VWC_raw <- swVWC$VWC * Control_sd + Control_mean

# make time periods --------------------------------------------------------------------

p1 <- data.frame( Period = 'Modern', year = 2007:2016)
p2 <- data.frame( Period = 'not monitored', year = 1958:2006)
p3 <- data.frame( Period = 'Historical', year = 1925:1957)
periods <- data.frame( rbind( p1, p2, p3 )) 

# set-up aggregate seasonal variables for model ----------------------------------------#

swVWC <- 
  swVWC %>% 
  ungroup() %>% 
  mutate( water_year = year + lag_year ) %>% 
  dplyr::select(year, month, season, season_label, precip_seasons, water_year, Treatment, date, VWC, VWC_raw)

swVWC <- 
  swVWC %>% 
  left_join(periods, by = 'year')


saveRDS(swVWC, daily_swVWC_treatment_outfile)

#  Additional figures ---------------------------------------- # 

pdf(
  fig2_file,
  width = 8,
  height = 6
)

print(
  ggplot(swVWC, aes(
    x = date, y = VWC_raw, color = Treatment
  )) +
    geom_line() +
    scale_color_manual(values = my_colors[2:4]) +
    xlim(ymd(c(
      '2016-01-01', '2016-10-01'
    )))
)
dev.off()
print(
  ggplot(swVWC, aes(
    x = date, y = VWC, color = Treatment
  )) +
    geom_line() +
    scale_color_manual(values = my_colors[2:4]) +
    xlim(ymd(c(
      '2014-01-01', '2015-01-01'
    )))
)

