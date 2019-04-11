# calculate treatment effects compared to control soil moisture

rm(list = ls())

library(tidyverse)
library(lme4)
library(zoo)
library(MASS)
library(emmeans)
library(texreg)
library(xtable)

swVWC <- read.csv('data/daily_VWC.csv')                       # soilwat output 
spotVWC <- readRDS('data/temp_data/spotVWC.RDS')              # spot check soil moisture

myVWC <- readRDS('data/decagon_data_with_station_data.RDS')   # cleaned decagon data 
daily_clim <- readRDS('data/daily_station_dat_rainfall.RDS')  # climate station rainfall

seasons <- read.csv('data/season_table.csv')

# --------------------------------------------------------------------------------------#
myVWC <- myVWC %>%
  mutate(v = ifelse(measure == 'VWC', v * 100, v)) # convert to percent

myVWC$year <-
  as.numeric(strftime(as.Date(myVWC$simple_date), '%Y'))

# Steps for treatment effects standardization:
# 1. aggregate soil moisture by PrecipGroup, Treatment, depth and date
# 2. then standardize soil moisture measurements within PrecipGroup and depth
# 3. then subtract standardized Drought and Irrigation from Control to get
#       standardized treatment effects.
# 4. output

VWC_weights <-
  myVWC %>%
  filter(bad_values == 0 , stat == 'raw', measure == 'VWC') %>%
  filter(!is.na(rainfall)) %>%
  filter(!(plot == 16)) %>%           ##### Drop plot 16
  group_by(year, PrecipGroup, simple_date) %>%
  summarise(weight = n_distinct(unique_position))

myVWC <-
  myVWC %>%
  filter(bad_values == 0 , stat == 'raw', measure == 'VWC') %>%
  filter(!is.na(rainfall)) %>%
  filter(!(plot == 16)) # drop plot 16

df_soil_moist <-
  myVWC %>%
  group_by(year, PrecipGroup, Treatment, season, simple_date, rainfall) %>%
  summarise(avg_VWC = mean(v, na.rm = TRUE)) %>%
  group_by(PrecipGroup) %>%
  mutate(avg_VWC = scale(avg_VWC, mean(avg_VWC[Treatment == 'Control'], na.rm = T), sd(avg_VWC[Treatment == 'Control'], na.rm = T))) %>%  # scale within Precip Group and Depth
  spread(Treatment, avg_VWC) %>%
  mutate(Drought = Drought - Control, Irrigation = Irrigation - Control) %>%
  arrange(PrecipGroup, simple_date)

df_soil_moist <- merge(df_soil_moist, VWC_weights)

# -------- add in the spot measurements ------------------------------------------------#
df_soil_moist$type <- 'logger'
VWC_df <- df_soil_moist

spotVWC$type <- 'spot'
spotVWC$year <- strftime(spotVWC$date, '%Y')

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
  stepAIC(mTreatment,
          scope = list(upper = ~ . , lower = ~ 1),
          trace = T)
summary(mTreatment)

mTreatment <-
  lmer(update(formula(mTreatment) , . ~ . + (1 | simple_date) + (1 | PrecipGroup)),
       data = VWC_test,
       weights = VWC_test$weight)

summary(mTreatment)


test <- lsmeans(mTreatment,  ~ Treatment + season + rainfall)
#test <- lsm(mTreatment, ~ Treatment + season + rainfall )


tab <- summary(test)

tab <-  data.frame(tab)

tab <-
  tab %>%
  dplyr::select(season, rainfall,  Treatment, lsmean, SE, asymp.LCL, asymp.UCL) %>%
  arrange(season, rainfall, Treatment)

statsOutput <- 'manuscript/soil_moisture_model.tex'

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
  dplyr::select(simple_date, year, season, rainfall, Treatment, v) %>%
  group_by(Treatment, simple_date, year, season, rainfall) %>%
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
  predict(mTreatment, pred_df, re.form = NA)

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

plot_df <- merge(pred_df, observed_df)

plot_df <-
  plot_df %>%
  mutate(back_scaled_pred = predicted * sd(observed[Treatment == 'Control']) + mean(observed[Treatment == 'Control']))

head(plot_df)

plot_df <-
  plot_df %>% 
  dplyr::select(-predicted) %>% 
  distinct() %>% 
  rename(predicted = back_scaled_pred)  %>% 
  gather(type, VWC, predicted , observed)

everyday <-
  expand.grid(
    simple_date = seq.Date(as.Date('2012-01-01'), as.Date('2016-12-30'), 1),
    Treatment = c('Control', 'Drought', 'Irrigation'),
    type = c('predicted', 'observed')
  )

plot_df <- merge(everyday, plot_df, all.x = T)

plot_df <-
  plot_df %>% mutate(julian_date = as.numeric(strftime(simple_date, '%j')),
                     year = as.numeric(strftime(simple_date, '%Y')))

load('code/figure_scripts/my_plotting_theme.Rdata')

subset(plot_df, type != 'predicted')

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
  'figures/avg_daily_soil_moisture.png',
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

swVWC$date <-
  as.POSIXct(strptime(paste(swVWC$Year, swVWC$DOY, sep = '-') , '%Y-%j'))
swVWC$year <- swVWC$Year
swVWC$month <- as.numeric(strftime(swVWC$date, '%m'))

lyr3 <-
  swVWC %>% group_by(month, year) %>% summarise(avg = mean(Lyr_3))
lyr3 <- lyr3 %>% spread(month, avg)
lyr3 <- lyr3[complete.cases(lyr3),]
mydata <- lyr3[, 2:13]
mydata <- scale(mydata)

pca <- princomp(mydata)
biplot(pca)
pca$loadings

# set-up aggregate seasonal variables for model ----------------------------------------#

swVWC <-
  swVWC %>%
  gather(layer, VWC, Lyr_1:Lyr_6) %>%
  filter(layer %in% c('Lyr_1', 'Lyr_2', 'Lyr_3', 'Lyr_4'))

swVWC <- swVWC %>%
  group_by(date) %>%
  summarise(modelVWC = mean(VWC) * 100)

daily_clim$date <- as.Date(daily_clim$date)
swVWC$date <- as.Date(swVWC$date)
swVWC$month <- as.numeric(strftime(swVWC$date, '%m'))
swVWC$year <- as.numeric(strftime(swVWC$date, '%Y'))

daily_clim <-
  daily_clim %>%
  ungroup() %>%
  dplyr::select(-year)

swVWC <- left_join(swVWC, seasons, by = 'month')
swVWC <- left_join(swVWC, daily_clim, by = c('date'))

head(swVWC)
head(daily_VWC2)

soilWAT <-
  swVWC %>% dplyr::select(simple_date, modelVWC, year, season, rainfall)
soilWAT$SW_predicted <- soilWAT$modelVWC
daily_VWC2$observed <- daily_VWC2$Control

obs_predicted <-
  merge(daily_VWC2,
        soilWAT,
        by = c('simple_date', 'year', 'season', 'rainfall'))
ggplot(obs_predicted, aes(x = SW_predicted, y = observed)) + geom_point()

obs_predicted <-
  obs_predicted %>% gather(type, val, observed, SW_predicted)
ggplot(obs_predicted, aes(x = simple_date, y = val, color = type)) + geom_line()


swVWC$Control <- scale(swVWC$modelVWC) # standardize control SWC
Control_mean <- mean(swVWC$modelVWC)
Control_sd <- sd(swVWC$modelVWC)

swVWC2 <- swVWC
swVWC$Treatment <-  'Drought'
swVWC2$Treatment <- 'Irrigation'

swVWC <- rbind(swVWC, swVWC2)

swVWC$predicted <- predict(mTreatment,  swVWC, re.form = NA)

swVWC <- swVWC %>% spread(Treatment, predicted)

# unscale the Control Drought and Irrigation VWC --------------------------------------------------------  #

swVWC <-
  swVWC %>%
  mutate(Drought = Control + Drought,
         Irrigation = Control + Irrigation) %>%
  gather(Treatment, VWC, Control:Irrigation)

swVWC$VWC_raw <- swVWC$VWC * Control_sd + Control_mean

pdf(
  'figures/modeled_soilwat_soil_moisture_example.pdf',
  width = 8,
  height = 6
)

print(
  ggplot(swVWC, aes(
    x = date, y = VWC_raw, color = Treatment
  )) +
    geom_line() +
    scale_color_manual(values = my_colors[2:4]) +
    xlim(as.Date(c(
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
    xlim(as.Date(c(
      '2014-01-01', '2015-01-01'
    )))
)


saveRDS(swVWC, 'data/temp_data/daily_swVWC_treatments.RDS')
