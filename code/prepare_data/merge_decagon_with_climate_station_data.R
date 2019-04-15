rm(list = ls())

library( tidyverse )
library(lubridate)
# input ---------------------------------------------------- #

decagon <- sheepweather::usses_decagon # comes from sheepweather package 

rainfall <- readRDS('data/temp_data/daily_station_dat_rainfall.RDS')
  # comes from the make_rainfall script

# output ---------------------------------------------------- #

decagon_outfile <- 'data/temp_data/decagon_data_with_station_data.RDS'

# ---------------------------------------------------------------------------------------

decagon <-
  decagon %>%
  left_join(rainfall, by = 'date')

saveRDS(decagon, decagon_outfile )


