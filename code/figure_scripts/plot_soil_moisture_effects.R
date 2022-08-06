rm(list = ls())
library(tidyverse)
library(lubridate)
library(lme4)

load('code/figure_scripts/my_plotting_theme.Rdata')

# output --------------------------------------------------------------- # 
fig1_file <- 'figures/avg_daily_soil_moisture.png'
fig2_file <- 'figures/modeled_soilwat_soil_moisture_example.pdf'

# input ---------------------------------------------------------------- # 
sm_model <- readRDS('data/temp/treatment_sm_model.rds')
daily_sm <- read_csv('data/temp/daily_sm.csv')
weather <- read_csv('data/temp/daily_station_dat_rainfall.csv')
seasons <- read_csv('data/season_table.csv')
# ---------------------------------------------------------------------- # 

daily_avg_sm <- 
  daily_sm %>% 
  group_by(date, Treatment) %>% 
  summarise( obs =  mean(v), raw = mean(raw)) 

rain_data <- 
  weather %>%
  ungroup() %>% 
  filter(date > '2012-01-01', date < '2017-01-01') %>% 
  select( date, inv_days_since_rain) 
  
pred_df <- 
  rain_data %>% 
  mutate(Treatment = 'Control') %>% 
  bind_rows(rain_data %>% 
              mutate( Treatment = 'Irrigation')) %>% 
  bind_rows(rain_data %>% 
              mutate( Treatment = 'Drought')) %>% 
  left_join(daily_avg_sm, by = c('date', 'Treatment')) %>% 
  mutate(month = month(date)) %>% 
  left_join(seasons %>% select(month, season), by = 'month') %>% 
  data.frame()

pred_df$pred <- predict(sm_model, newdata = pred_df, re.form = NA )

ambient_sm <- 
  pred_df %>% 
  filter( Treatment == 'Control') %>% 
  select( date, raw)   %>% 
  distinct()

predicted_soil_moisture <- 
  pred_df %>% 
  select( date, Treatment, season, obs, pred) %>% 
  arrange( date) %>% 
  gather(type, value, c('obs', 'pred')) %>%  
  left_join(ambient_sm, by = 'date') %>% 
  mutate( value_bt = value*raw) %>% 
  mutate( DOY = yday(date), year = year(date))


gg_predicted <- 
  predicted_soil_moisture %>% 
  filter(Treatment != 'Control') %>% 
  ggplot( aes( x = DOY, y = value_bt, linetype = type, color = Treatment)) + 
  geom_line() + 
  facet_wrap(~ year, nrow = 5) + 
  scale_color_manual(values = my_colors[3:4]) +
  scale_linetype_manual(values = c(1, 2)) +
  scale_alpha_manual(values = c(1, 0.7)) +
  my_theme + 
  ylab('Volumetric soil moisture (ml/ml)') + 
  xlab('Day of year')

ggsave(plot = gg_predicted,filename =fig1_file, width = 7, height = 7)


