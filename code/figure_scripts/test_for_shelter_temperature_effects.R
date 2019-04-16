rm(list = ls())

library(lubridate)
library(tidyverse)
library(lme4)
library(emmeans)
# load data ------------------------------------------------ # 

library(sheepweather)  # get weather, decagon, ibutton

quads   <- read_csv('data/quad_info.csv')
seasons <- read_csv('data/season_table.csv')

# ---------------------------------------------------------- # 

air_temps <-
  usses_decagon %>%
  filter( measure == 'C' & position == 'air') %>%
  mutate( hour = hour(datetime), month = month(date) , year = year(date))

soil_temps <-
  usses_decagon %>%
  filter( measure == 'C' & position != 'air') %>%
  mutate( hour = hour(datetime))

daily_air <-
  air_temps %>%
  distinct(datetime, date, plot, v) %>%
  group_by(date, plot) %>%
  summarise( TMAX  = max(v), TMIN = min(v), TMEAN = (TMAX + TMIN)/2)  %>%
  gather( stat, decagon , TMAX, TMIN, TMEAN )

daily_soil <-
  soil_temps %>%
  filter( depth == '5 cm deep') %>%
  distinct( datetime, date, plot, v) %>%
  group_by(date, plot) %>%
  summarise( TMAX  = max(v), TMIN = min(v), TMEAN = (TMAX + TMIN)/2)  %>%
  gather( stat, decagon , TMAX, TMIN, TMEAN )

station <-
  usses_weather %>%
  filter(  ELEMENT != 'PRCP') %>%
  spread( ELEMENT, value ) %>%
  mutate( TMEAN = (TMAX + TMIN)/2 )  %>%
  gather( stat, station, TMAX, TMIN, TMEAN)


daily_anoms <-
  station %>%
  left_join(
    bind_rows(
      daily_air %>% rename('value' = decagon) %>% mutate(type = 'decagon'),
      usses_ibutton %>% rename( 'value' = ibutton)  %>% mutate( type = 'ibutton') %>% select(-n) ),
    by = c('date', 'stat')) %>%
  filter( !is.na(plot)) %>%
  mutate( anom = value - station ) %>%
  filter( stat == 'TMEAN', abs(anom) < 20 ) %>%
  mutate( month = month(date), year = year(date)) %>%
  left_join(quads, by = c('plot' = 'QuadName')) %>%
  left_join(seasons, by = 'month') %>%
  select( month, plot, Treatment, type, year, season, date, anom, PrecipGroup)



m1 <- lmer( data = daily_anoms, anom ~ Treatment + type + (1|plot) + (1|date))
summary(m1)
summary( m1 )

gg1 <- emmeans(m1, ~ Treatment ) %>%
  data.frame() %>%
  ggplot( aes(  x = Treatment , y = emmean, ymin = asymp.LCL, ymax = asymp.UCL )) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_errorbar(position = position_dodge(width = 0.2)) +
  ylab( 'Daily Plot Tmean - Station Tmean deg. C') +
  theme( axis.text.x = element_text(size = 14 ), axis.title.x = element_blank())

m1 <- lm( data = daily_anoms, anom ~ Treatment + type )

gg2 <- emmeans(m1, ~ Treatment ) %>%
  data.frame() %>%
  ggplot( aes(  x = Treatment , y = emmean, ymin = lower.CL, ymax = upper.CL )) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_errorbar(position = position_dodge(width = 0.2)) +
  ylab( 'Daily Plot Tmean - Station Tmean deg. C') +
  theme( axis.text.x = element_text(size = 14 ), axis.title.x = element_blank())

myplot <- egg::ggarrange(gg1 + ylim(-1,1) + ggtitle('Random effects: (1|plot) + (1|date)'),
               gg2 +
                 theme(axis.title.y = element_blank()) +
                 ylim(-1, 1) +
                 ggtitle('No random effects'), nrow = 1)

ggsave(myplot, filename = 'figures/shelter_temp_effect.png', width = 8, height = 5)
