
library(tidyverse)
load('data/temp_data/st_anom.rda')

biomass <- c(472, 308, 367, 379, 437, 310)
biomass <- data.frame( year = 1942:1947, biomass = biomass ) 

biomass_prcp <- 
  st_anom %>% 
  filter( year > 1939, year < 1950 ) %>% 
  filter( end == 7, len == 12, var == 'PRCP') %>% 
  left_join(biomass)

biomass_prcp %>% 
  ggplot( aes( x = year, y = rllm)) + 
  geom_line( aes( y = biomass/10)) + 
  geom_line(color = 'blue') + 
  scale_x_yearmon()

biomass_prcp %>% 
  ggplot( aes( x = rllm, y= biomass/10 )) + 
  geom_point()

biomass_prcp <- biomass_prcp[ complete.cases(biomass_prcp), ] 
cor(biomass_prcp$rllm, biomass_prcp$biomass)

species <- 'PSSP'

intra <- paste0( 'W.' , species )
dat <- readRDS( paste0 ( 'data/temp_data/', paste0( species, '_size.RDS')))
dat <- dat[ dat$year < 1960 , ]



# primitave Area Lat/Lon 44.240660, -112.176595
PA_lat <- 44.240660
PA_lon <- -72.179003
library(MODISTools)
# subset <- mt_subset(product = "MOD13Q1",
#                     lat = PA_lat,
#                     lon = PA_lon,
#                     band = "250m_16_days_NDVI",
#                     start = "2004-01-01",
#                     end = "2019-12-30",
#                     km_lr = 0,
#                     km_ab = 0,
#                     site_name = "testsite",
#                     internal = TRUE,
#                     progress = TRUE)

save(subset, file = 'data/ndvi.rda')

subset %>% 
  mutate( date = lubridate::ymd( calendar_date) ) %>% 
  ggplot( aes( x = date, y = value )) + 
  geom_line()
