rm(list = ls())

library( tidyverse )
library(zoo) # for rollapply

# input ---------------------------------------------------- #

weather <- sheepweather::usses_weather # use sheepweather package 

# output ---------------------------------------------------- #

rainfall_outfile <- 'data/temp_data/daily_station_dat_rainfall.RDS'

# ---------------------------------------------------------------------------------------

weather$date <- as.POSIXct( strptime( weather$date, '%Y-%m-%d', tz = 'MST') )

weather <-
  weather %>%
  spread( ELEMENT, value)

weather <-
  weather %>%
  mutate( TMEAN = ( TMAX + TMIN ) / 2 ) %>%
  select(date, PRCP, TMEAN) %>%
  mutate( rainfall = rollapply(PRCP, 2, sum, fill = 0, na.rm = TRUE, align = 'right') ) %>%
  mutate( rainfall = ifelse( rainfall > 0.0 & TMEAN > 3 & !is.na(rainfall), 'rainy', 'not rainy')) %>%
  mutate( rainfall = ifelse( is.na(rainfall), 'not rainy', rainfall))

# create a factor listing each rainy period, including the day before the rain

weather <-
  weather %>%
  arrange( desc(date) ) %>%
  mutate( prerain = lag( rainfall, 1) ) %>%
  mutate( prerain = ifelse( prerain == 'rainy' & rainfall == 'not rainy', TRUE, FALSE)) %>%
  arrange( date) %>%
  mutate( prcp_event = factor( cumsum ( prerain ) )) %>%
  group_by( prcp_event, prerain) %>%
  mutate( total_rain = cumsum(ifelse(is.na(PRCP), 0, PRCP) ))

weather <-
  weather %>%
  ungroup() %>%
  mutate( simple_date = as.Date( date, tz = 'MST')) %>%
  mutate( year = strftime( simple_date, '%Y', tz = 'MST')) %>%
  group_by( year ) %>%
  arrange( year, simple_date ) %>%
  mutate( ann_cum_PRCP = ifelse(is.na(PRCP), 0, PRCP))

saveRDS(weather, rainfall_outfile )
