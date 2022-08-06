rm(list = ls())

library( tidyverse )
library(lubridate)
library(sheepweather)
# input ---------------------------------------------------- #

decagon <- usses_decagon # comes from sheepweather package 

rainfall <- read_csv('data/temp/daily_station_dat_rainfall.csv')
# comes from the make_rainfall script

# output ---------------------------------------------------- #

decagon_outfile <- 'data/temp/decagon_data_with_station_data.RDS'

# ---------------------------------------------------------------------------------------

decagon <-
  decagon %>%
  left_join(rainfall, by = 'date')

saveRDS(decagon, decagon_outfile )


