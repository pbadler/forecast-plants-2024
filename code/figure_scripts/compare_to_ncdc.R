
# compare to ncdc 
weather %>% filter( ELEMENT == 'TMAX', is.na(value))

daily_ncdc <- read_csv('~/Desktop/1701628.csv')
daily_ncdc$date <- ymd( daily_ncdc$DATE )

weather %>% 
  mutate( year = year(date)) %>% 
  group_by( year, ELEMENT) %>% 
  summarise( n() )

monthly_ncdc_from_daily <- daily_ncdc %>% 
  mutate( year = year(date), month = month(date)) %>%
  group_by( year, month) %>% 
  mutate( TAVG = (TMIN + TMAX)/2 ) %>% 
  summarise( TAVG = mean(TAVG, na.rm  = T),  PRCP = sum(PRCP, na.rm = T))

joined_daily <- 
  weather %>% left_join(
    daily_ncdc %>% 
      gather( ELEMENT, value, PRCP:TMIN), by = c('date', 'ELEMENT')) %>% 
  select( date, ELEMENT, value.x, value.y) 

joined_daily %>% 
  ggplot(aes( x = value.x, y = value.y)) + 
  geom_point()

joined_daily %>% 
  filter( ELEMENT == "PRCP") %>% 
  filter( value.x == 0 ) 

missing_prcp <- 
  joined_daily %>% 
  filter( ELEMENT == 'PRCP' ) %>% 
  group_by( month = month(date), year  = year(date) ) %>% 
  summarise( missing = sum(is.na(value.x))) %>% 
  filter( missing > 0) %>% 
  filter( year < 2019) 


ncdc_monthly <- read_csv('~/Desktop/1701639.csv')

ncdc_monthly <- 
  ncdc_monthly %>% select( DATE, PRCP, TAVG) %>% 
  mutate( year = str_extract(DATE, '[0-9]{4}'), 
          month = str_extract(DATE, '[0-9]{2}$')) %>% 
  mutate( YEAR = as.numeric(year)) %>% 
  mutate( MONTH = as.numeric(month)) %>% 
  select( YEAR, MONTH, PRCP, TAVG) %>% 
  gather( ELEMENT, value, PRCP:TAVG)

monthly_ncdc_from_daily  
ncdc_monthly

joined_test <- 
  monthly %>% 
  gather( ELEMENT, value, PRCP:TAVG) %>% 
  left_join(ncdc_monthly, by = c('YEAR', 'MONTH', 'ELEMENT')) 

joined_test %>% 
  filter( is.na(value.y))

joined_test %>% 
  filter( ELEMENT == 'PRCP') %>% 
  mutate( diff = value.x - value.y) 

