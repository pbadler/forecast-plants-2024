
#library(prism)
library(lubridate)
library(imputeTS)
library(pollen)

sheepweather::usses_weather %>% head 
daily %>% head 
sheepweather::usses_soilwat

daily <- read_csv('~/Dropbox/projects/sheepweather/data-raw/daily_station_data_ncdc_download.csv')
daily

load('~/Dropbox/projects/sheepweather/data/usses_soilwat.rda')  

daily <- 
  data.frame( DATE = seq( min( date( daily$DATE ) ), max(date(daily$DATE)), by = 1) ) %>% 
  mutate ( year = year( DATE)) %>% 
  left_join( daily %>% mutate( DATE = ymd(DATE)) , by = 'DATE') %>% 
  mutate( PPT_ts = ts(PRCP) , TMAX_ts = ts( TMAX), TMIN_ts = ts(TMIN )) %>%
  mutate( PPT_ts = na_interpolation(PPT_ts), 
          TMIN_ts = na_interpolation(TMIN_ts), 
          TMAX_ts = na_interpolation(TMAX_ts)) %>% 
  ungroup() %>% 
  select(DATE, PPT_ts, TMAX_ts, TMIN_ts ) 
  
usses_soilwat <- 
  usses_soilwat %>% 
  mutate(d = row_number() -1 ) %>% 
  mutate( DATE = as.Date( d, origin = '1925-01-01')) %>% 
  mutate( Lyr_1 = ts( Lyr_1), Lyr_2 = ts(Lyr_2)) %>% 
  select( DATE, Lyr_1, Lyr_2)

daily <- 
  daily %>% 
  left_join(usses_soilwat) %>% 
  mutate( WEEK = lubridate::week(DATE), YEAR = year( DATE ), MONTH = month(DATE ) ) 

monthly <- 
  daily %>% 
  group_by(YEAR, MONTH) %>% 
  mutate( GDD = pollen::gdd( TMAX_ts, TMIN_ts, tbase = 10, tbase_max = 25)) %>% 
  mutate( growdays1 = TMAX_ts > 10 & TMAX_ts < 25 & Lyr_1 > 0.13, growdays2 = TMAX_ts > 10 & TMAX_ts < 25 & Lyr_2 > 0.13 ) %>% 
  summarise( PPT = sum(PPT_ts), TMAX_ts = mean(TMAX_ts), TMIN_ts = mean(TMIN_ts), wd1 = sum(Lyr_1 > 0.125), wd2 = sum(Lyr_2 > 0.125), GDD = max(GDD), growdays1 = sum(growdays1), growdays2 = sum(growdays2))

weekly <- 
  daily %>% 
  group_by(YEAR, WEEK) %>% 
  mutate( GDD = pollen::gdd( TMAX_ts, TMIN_ts, tbase = 10, tbase_max = 25)) %>% 
  mutate( growdays1 = TMAX_ts > 10 & TMAX_ts < 25 & Lyr_1 > 0.13, growdays2 = TMAX_ts > 10 & TMAX_ts < 25 & Lyr_2 > 0.13 ) %>% 
  summarise( PPT = sum(PPT_ts), TMAX_ts = mean(TMAX_ts), TMIN_ts = mean(TMIN_ts), wd1 = sum(Lyr_1 > 0.125), wd2 = sum(Lyr_2 > 0.125), GDD = max(GDD), growdays1 = sum(growdays1), growdays2 = sum(growdays2))


weekly 


#  mutate( GDD = pollen::gdd( TMAX_ts, TMIN_ts, tbase = 10, tbase_max = 25 )) 


usses_soilwat %>% 
  filter( Date > '1940-01-01' & Date < '1944-01-01') %>% 
  ggplot( aes( x = Date, y = Lyr_2)) + geom_line() + 
  geom_hline( aes( yintercept = 0.125)) + 
  scale_x_date(date_breaks = '6 month')

summer_wet_days <- 
  usses_soilwat %>% 
  group_by( Year ) %>% 
  filter( month( Date) %in% c( 5:9 )) %>% 
  summarise( Lyr1 = sum( Lyr_1 > 0.125), Lyr2 = sum( Lyr_2 > 0.125)) %>% 
  gather( layer, wetdays , starts_with('Lyr')) 





  ggplot( aes( x = Year, y= wetdays, color = layer)) + 
  geom_line() 

ET_dat <- list( Tmax = zoo::zoo(test$TMIN_ts), Tmin = zoo::zoo ( test$TMAX_ts), Date.daily = zoo::as.Date( test$DATE ))

constants$lat
constants$Elev
constants$lat_rad

ET_dat$Date.monthly <- as.yearmon(as.character(test$DATE))
ET_dat$Date.daily <- ET_dat$Date.daily
ET_dat$Tmax <- ET_dat$Tmax
ET_dat$Tmin <- ET_dat$Tmin

data( "processeddata")
constants$Elev <- 1661
constants$lat_rad <- 0.7710218863
constants$lambda
constants$Gsc 

constants[! names( constants ) %in% c('Elev', 'lat_rad', 'lambda', 'Gsc')  ] <- NA 

ET_dat$J <- zoo( lubridate::yday(ET_dat$Date.daily) ) 
ET_dat$i <- zoo( lubridate::month(ET_dat$Date.daily) ) 

ET_dat$Ndays <- 
  zoo( 
    test %>% group_by( year, month( DATE)) %>% 
      summarise( ndays = n() ) %>% .$ndays )

ET_dat$Tmax[ ET_dat$Tmax < 5 ] <- 4 
ET_dat$Tmin[ ET_dat$Tmin < 4 ] <- 2 
out <- ET.HargreavesSamani( ET_dat, constants, ts = 'daily')



out <- ET.HargreavesSamani( test_et  , constants)
constants

test %>% 
  filter( year %in% c( 1925:1940 )) %>% 
  ggplot( aes(x = DATE, y = gdd )) + 
  geom_line() +  
  geom_point( data = test %>% 
                filter( year %in% c(1925:1940)) %>% 
                filter( month(DATE) == 6, day(DATE) == 1) ) 

test %>% 
  filter( year %in% c( 1935:1940 )) %>% 
  ggplot( aes(x = DATE, y = PPT_ts )) + 
  geom_point() 
  

hist( ( test$PRCP)  )

USSES_ppt <- test %>% 
  mutate( date = ymd (DATE )) %>% 
  mutate( year = year( date), month = month(date)) %>% 
  group_by( year, month ) %>% 
  summarise( USSES = mean ( PRCP ))

prism_extract <- read_csv('~/Downloads/PRISM_ppt_tmin_tmean_tmax_tdmean_vpdmin_vpdmax_stable_4km_190101_201701_44.2410_-112.1759.csv', skip = 10) 

prism_ppt <- 
  prism_extract %>% 
  separate(Date, c('year', 'month')) %>% 
  mutate( year = as.numeric( year ) , month = as.numeric(month)) %>% 
  mutate( date = ymd ( paste( year, month, 1 )) ) %>% 
  mutate( year = year( date) , month = month( date)) %>% 
  mutate( PRISM = `ppt (mm)`) %>% 
  select( year, month, PRISM)


prism_ppt  %>% 
  left_join( USSES_ppt ) %>% 
  mutate( USSES = USSES*25.4 ) %>% 
  ggplot( aes( x = PRISM, y= USSES )) + 
  geom_point() + 
  geom_abline(aes( intercept = 0 , slope = 1))


#  ggplot( aes( x = date, y = `ppt (mm)`)) + 
#  geom_line() 


